import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
# import matplotlib
# matplotlib.use('TkAgg')     # https://stackoverflow.com/questions/28757348/how-to-clear-memory-completely-of-all-matplotlib-plots
from reproject.mosaicking import find_optimal_celestial_wcs, reproject_and_coadd
from reproject import reproject_interp
from astropy.visualization.wcsaxes import SphericalCircle
from astropy import units as u
import abell_cluster_module as ab
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.visualization import make_lupton_rgb
from astropy.io import ascii
from astropy import wcs
import my_module as mm
import img_scale
import importlib
importlib.reload(img_scale)
importlib.reload(mm)
import pylab as py
import sys
from astropy.convolution import Gaussian2DKernel
from scipy.signal import convolve as scipy_convolve
from astropy.convolution import convolve
from astropy.cosmology import Planck18 as Cosmo
from scipy import constants as const
import os
import calc_kcor as kcor

qs = [8 + exp * 4 for exp in range(3, 4)]
# qs = [10]
stretches = [0.08 + 0.02 * exp for exp in range(-2, -1)]
# stretches = [0.08]
x_stddev = [0.5 * exp for exp in range(1, 2)]
n = 3

prog = ['SDSS', 'DECaLS']
pixscale = [0.4, 0.262]

# save_dir = ab.work_dir + 'gui/A3558/'
cmap = 'PuOr'

def find_and_cut(cluster, program, coord, hthumb):
    # read the celestial coordinate
    cel_coord = [[coord.ra.value,
                  coord.dec.value], [0, 0]]
    for jj in range(1, 26):
        with fits.open(ab.work_dir + f'fits/{prog[program]}_extracted/{ab.clusters[k]}_rsi_{jj}.fits') as hdu_r:
            # read WCS
            w = wcs.WCS(hdu_r[0].header)
            pixcrd = w.wcs_world2pix(cel_coord, 1)
            # if w.footprint_contains(sky_coord):
            yy = int(pixcrd[0][0])
            xx = int(pixcrd[0][1])
            xb = hdu_r[0].shape[0]
            yb = hdu_r[0].shape[1]
            ht = hthumb

            if (xx > 0) & (xx < xb) & (yy > 0) & (yy < yb):
                y1 = yy - ht if yy - ht > 0 else 0
                y2 = yy + ht if yy + ht < yb else yb
                x1 = xx - ht if xx - ht > 0 else 0
                x2 = xx + ht if xx + ht < xb else xb

                if program:     # DECaLS
                    with fits.open(ab.work_dir + f'fits/{prog[program]}_extracted/{ab.clusters[k]}_gsi_{jj}.fits') as hdu_g:
                        g_array, footprint = reproject_interp(hdu_g[0], hdu_r[0].header)
                        cutout_u = g_array[x1:x2, y1:y2] * 0        # no data in u
                        cutout_g = g_array[x1:x2, y1:y2]
                        cutout_r = hdu_r[0].data[x1:x2, y1:y2]
                        # print(f'{ht} {x1} {x2} {y1} {y2} {xb} {yb}')
                        return cutout_u, cutout_g, cutout_r, jj, xx, yy   # FITS file [y, x]
                else:           # SDSS
                    with fits.open(ab.work_dir + f'fits/{prog[program]}_extracted/{ab.clusters[k]}_gsi_{jj}.fits') as hdu_g, \
                            fits.open(ab.work_dir + f'fits/{prog[program]}_extracted/{ab.clusters[k]}_usi_{jj}.fits') as hdu_u:
                        g_array, footprint = reproject_interp(hdu_g[0], hdu_r[0].header)
                        u_array, footprint = reproject_interp(hdu_u[0], hdu_r[0].header)
                        cutout_u = np.rot90(hdu_r[0].data[x1:x2, y1:y2], k=3)
                        cutout_g = np.rot90(g_array[x1:x2, y1:y2], k=3)
                        cutout_r = np.rot90(u_array[x1:x2, y1:y2], k=3)
                        # print(f'{ht} {x1} {x2} {y1} {y2} {xb} {yb}')
                        return cutout_u, cutout_g, cutout_r, jj, xx, yy   # FITS file [y, x]
    return 0, 0, 0, 0, 0, 0

for k in range(2, 3):   # 2399 and 2670 only which are covered by SDSS and DECaLS
    # save_dir = ab.work_dir + f'pics/{ab.clusters[k]}_new_rs/'
    save_dir = ab.work_dir + f'pics/vs_sdss_decals/'

    with open(ab.work_dir + f'spec/{ab.clusters[k]}_spec_ind.npy', 'rb') as spec_ind, \
        open(ab.work_dir+f'spec/{ab.clusters[k]}_rs_each_ind.npy', 'rb') as rs_ind:  \

        sex_cat = ascii.read(ab.sex_dir + f'DECam_merged_SEx_cat_{ab.clusters[k]}_Gal_ext_corrected_20rmag_psf.txt')
        ind_spec = np.load(spec_ind)
        ind_rs = np.load(rs_ind)
        ind_xor = ind_spec ^ ind_rs
        ind_spec_only = ind_xor & ind_spec
        ind_spec_or_rs = ind_spec | ind_rs
        cat_rs = sex_cat[ind_spec_or_rs]

        nan_id = []
        # for i in range(0, 1):
        # for i in range(0, len(cat_rs)):
        i = np.where(cat_rs['NUMBER'] == 402)[0][0]

        c = SkyCoord(f"{cat_rs['ALPHA_J2000'][i]} {cat_rs['DELTA_J2000'][i]}", unit='deg')
        z = ab.redshifts[k]  # redshift
        kpc_arcmin = Cosmo.kpc_proper_per_arcmin(z)  # xxx kpc / arcmin
        arcmin_5kpc = 5.0 / kpc_arcmin  # xxx arcmin / 5 kpc

        for j in range(0, 2):   # 0: SDSS, 1: DECaLS
            half_thumb = int(10 * cat_rs['A_WORLD'][i] / pixscale[j])  # 4 x major axis as a box size for SDSS
            cut_u, cut_g, cut_r, tile_no, cut_x, cut_y = find_and_cut(ab.clusters[k], j, c, half_thumb)

            if tile_no:
                fig = plt.figure(figsize=(5, 5))
                img_sqrt_0_3 = np.zeros((cut_r.shape[0], cut_r.shape[1], 3), dtype=float)
                img_jarret = np.zeros((cut_r.shape[0], cut_r.shape[1], 3), dtype=float)
                minval = 0
                maxval = 1

                img_sqrt_0_3[:,:,0] = img_scale.sqrt(cut_r, scale_min=0, scale_max=0.1)
                img_sqrt_0_3[:,:,1] = img_scale.sqrt(cut_g, scale_min=0, scale_max=0.1)
                img_sqrt_0_3[:,:,2] = img_scale.sqrt(cut_u, scale_min=0, scale_max=0.1)

                img_jarret[:,:,0] = mm.jarrett(cut_r, np.nanstd(cut_r), n)
                img_jarret[:,:,1] = mm.jarrett(cut_g, np.nanstd(cut_g), n)
                img_jarret[:,:,2] = mm.jarrett(cut_u, np.nanstd(cut_u), n)

                plt.imshow(img_sqrt_0_3, aspect='equal', origin='lower')
                plt.savefig(save_dir + f'{cat_rs["NUMBER"][i]}_rgb_sqrt_0_01_{prog[j]}.png')
                plt.close(fig)
                plt.imshow(img_jarret, aspect='equal', origin='lower')
                plt.savefig(save_dir + f'{cat_rs["NUMBER"][i]}_rgb_jar_{n}_{prog[j]}.png')
                plt.close(fig)

                for q in qs:
                    for stretch in stretches:
                        fig = plt.figure(figsize=(5, 5))
                        array_rgb = make_lupton_rgb(cut_r, cut_g, cut_u, Q=q, stretch=stretch)
                        plt.imshow(array_rgb, origin='lower')
                        plt.savefig(save_dir + f'{cat_rs["NUMBER"][i]}_rgb_Q{q}_str{stretch}_{prog[j]}.png')
                        plt.close(fig)
