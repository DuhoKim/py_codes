"""
Functions for personal use
written by https://github.com/DuhoKim
"""

from astropy import coordinates as coords
import astropy.units as u
import numpy as np
from datetime import datetime
from astropy.cosmology import Planck18
from astropy.coordinates import Distance, SkyCoord
from astropy import wcs
import math
from astropy.io import ascii
import glob
import os
from os import path
import shutil

def closest(cat1, pos1):
    best_dist = 360.0 * 60.0 * u.arcminute
    best_index = 0
    for i in range(i, len(cat1)):
        dist = pos1.separation(cat1[i]).arcminute
        if dist < best_dist:
            best_dist = dist
            best_index = i
    return best_index, best_dist

def find_nearest(array, value):
    array = np.array(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx

# When there is a increasing or decreasing
# arrays of x corresponding to y,
# find a x0 that corresponds to y0 which is in the middle
def find_x0(x, y, y0):
    grt_y = np.where(y - y0 > 0)[0]  # index where y > y0
    if y[0] < y[-1]:                # increasing
        if len(grt_y) == 0:                     # y0 is smaller than all y
            x0 = x[0]
        elif grt_y[0] == len(x) - 1:       # y0 is greater than all y
            x0 = x[-1]
        else:
            a = y[grt_y][0] - y0
            b = y0 - y[grt_y][-1]
            e = x[grt_y][0] - x[grt_y][-1]
            dx = e * b / (a + b)
            x0 = x[grt_y][-1] + dx
    else:                                   # decreasing
        if grt_y[-1] == len(x)-1:                     # y0 is smaller than all y
            x0 = x[-1]
        elif len(grt_y) == 0:                 # y0 is greater than all y
            x0 = x[0]
        else:
            a = y[grt_y][-1] - y0
            b = y0 - y[grt_y[-1]+1]
            e = x[grt_y[-1]+1] - x[grt_y[-1]]
            dx = e * a / (a + b)
            x0 = x[grt_y[-1]] + dx
    return x0

def print_time():
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("Current Time =", current_time)

####################### Jarrett+03 ########################################################
#  Modified logarithmic visualization method P' = sqrt( log(1+P/(n*sig))  ), where P is the pixel
# intensity value, sig is the image rms "noise",
#  and n is a threshold throttle (w/ satisfactory values between 5 and 10
def jarrett(pix_values, sigma, n):
    return np.sqrt(np.log10(1+pix_values/(n*sigma)))
############################################################################################

def clcendist(z1, z2, coord1, coord2):
    dz1 = Distance(z=z1, cosmology=Planck18)
    dz2 = Distance(z=z2, cosmology=Planck18)
    c1 = SkyCoord(coord1, distance=dz1)
    c2 = SkyCoord(coord2, distance=dz2)
    sep = c1.separation_3d(c2)
    return sep.Mpc

def get_gala_sky(cluster, band, hdu_stack, coord):
    cel_coord = [[coord.ra.value, coord.dec.value], [0, 0]]
    tile_no = 0
    tile_no_last = 8 if cluster == 'A754' else 9
    for jj in range(1, len(hdu_stack)):
        # read WCS
        w = wcs.WCS(hdu_stack[jj].header)
        pixcrd = w.wcs_world2pix(cel_coord, 1)
        if (pixcrd[0][0] > 0) & (pixcrd[0][0] < hdu_stack[jj].shape[1])  &   \
                (pixcrd[0][1] > 0) & (pixcrd[0][1] < hdu_stack[jj].shape[0]):
            tile_no = jj
    if tile_no:
        dir_tile = f'/Users/duhokim/work/abell/galapagos/{cluster}_{band}/t{tile_no}/'
        if not os.path.isdir(dir_tile):
            while not os.path.isdir(dir_tile):
                tile_no = tile_no + 1 if tile_no < tile_no_last else 1
                dir_tile = f'/Users/duhokim/work/abell/galapagos/{cluster}_{band}/t{tile_no}/'
        rs_out = sorted(glob.glob(dir_tile + 'galfit/*_outsky'), key=os.path.getmtime, reverse=True)
        if len(rs_out) == 0:
            while len(rs_out) == 0:
                tile_no = tile_no + 1 if tile_no < tile_no_last else 1
                dir_tile = f'/Users/duhokim/work/abell/galapagos/{cluster}_{band}/t{tile_no}/'
                rs_out = sorted(glob.glob(dir_tile + 'galfit/*_outsky'), key=os.path.getmtime, reverse=True)
        rs_id = []
        for rs in rs_out:
            fn = rs.split('/')[-1]
            tid = fn.split('_')[0]
            rs_id.append(int(tid[2:])-1)
        rs_id_sorted = sorted(rs_id)
        outcat = ascii.read(dir_tile+f't{tile_no}r.outcat')
        outcat = outcat[rs_id_sorted]
        coords_out = SkyCoord(outcat['col13'], outcat['col14'], unit='deg')
        d2d = coord.separation(coords_out)
        dist = d2d.arcsec
        closest_id = outcat['col1'][np.argmin(dist)]
        res = open(dir_tile + f'galfit/t{tile_no}{closest_id}_r_outsky')
        first_line = res.readline()
        line_spl = first_line.split(' ')
        return float(next(s for s in line_spl if s))
    else:
        return math.nan

def copy_gala_result(cluster, hdu_stack, coord, copy_dir, cat_id):
    cel_coord = [[coord.ra.value, coord.dec.value], [0, 0]]
    tile_no = 0
    tile_no_last = 8 if cluster == 'A754' else 9
    for jj in range(1, len(hdu_stack)):
        # read WCS
        w = wcs.WCS(hdu_stack[jj].header)
        pixcrd = w.wcs_world2pix(cel_coord, 1)
        if (pixcrd[0][0] > 0) & (pixcrd[0][0] < hdu_stack[jj].shape[1])  &   \
                (pixcrd[0][1] > 0) & (pixcrd[0][1] < hdu_stack[jj].shape[0]):
            tile_no = jj
    if tile_no:
        ### check GALAPAGOS result of red sequence galaxies
        dir_tile = f'/Users/duhokim/work/abell/galapagos/{cluster}_r/t{tile_no}/'
        out_cat = ascii.read(dir_tile + f't{tile_no}r.outcat')
        coords_out = SkyCoord(out_cat['col13'], out_cat['col14'], unit='deg')
        d2d = coord.separation(coords_out)
        matched_sex = (d2d.arcsec < 1)
        if sum(matched_sex):
            galfit_id = out_cat['col1'][matched_sex][0]
            fn_img_gala = dir_tile + f'/galfit/png/t{tile_no}{galfit_id}_gala_vir_10000.png'
            if path.exists(fn_img_gala):
                shutil.copyfile(fn_img_gala, copy_dir + f'{cat_id}_gala.png')
                return
        else:
            return

        ### then check result of spectroscopic member
        dir_tile = f'/Users/duhokim/work/abell/galapagos/{cluster}_test/t{tile_no}/'
        out_cat = ascii.read(dir_tile + f't{tile_no}r.outcat')
        coords_out = SkyCoord(out_cat['col13'], out_cat['col14'], unit='deg')
        d2d = coord.separation(coords_out)
        matched_sex = (d2d.arcsec < 1)
        if sum(matched_sex):
            galfit_id = out_cat['col1'][matched_sex][0]
            fn_img_gala = dir_tile + f'/galfit/png/t{tile_no}{galfit_id}_gala_vir.png'
            if path.exists(fn_img_gala):
                shutil.copyfile(fn_img_gala, copy_dir + f'{cat_id}_gala.png')


