from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table, vstack
import numpy as np
import os
import shutil
import glob
from pandas import DataFrame, read_csv
import pandas as pd
from astropy.cosmology import Planck15 as Cosmo
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.backends.backend_pdf as pdf
from matplotlib.colors import LogNorm
import math
from astropy import constants as const
import my_module as mm
import math
import copy
import img_scale

work_dir = '/Users/duhokim/work/abell/galapagos/A3558_rsr_all/t1/galfit/'
cat = ascii.read("/Users/duhokim/work/abell/spec/Shapley/catalog.dat")

n = 7
cmap_perc_uni = ['viridis', 'plasma', 'inferno', 'magma', 'cividis']
cmap_seq = ['Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds', 'YlOrBr', 'YlOrRd', 'OrRd',
            'PuRd', 'BuPu', 'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn']
cmap_div = ['PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic']
cmap_mis = ['twilight', 'twilight_shifted', 'hsv', 'ocean', 'gist_earth', 'terrain', 'gist_stern', 'gnuplot', 'gnuplot2',
            'CMRmap', 'cubehelix', 'brg', 'gist_rainbow', 'rainbow', 'jet', 'turbo', 'nipy_spectral', 'gist_ncar']
cmap_fin = ['BrBG', 'cividis', 'coolwarm', 'jet', 'PiYG', 'PRGn', 'PuOr', 'RdBu', 'RdGy', 'RdYlBu', 'Spectral']

cmap = ['PuOr']

kron_fac = 2.2
rot_box = 1.5

def rotate(fits_f, cx, cy, pa_f, boa_f, rad_max):
    ylen, xlen = fits_f.shape
    pa_rad = math.radians(pa_f)

    yy, xx = np.ogrid[:ylen, :xlen]

    fits_f_rotated = copy.deepcopy(fits_f)

    ind1 = 1 <= ((yy - cy) * math.cos(pa_rad) - (xx - cx) * math.sin(pa_rad)) ** 2 / rad_max ** 2 + \
           ((yy - cy) * math.sin(pa_rad) + (xx - cx) * math.cos(pa_rad)) ** 2 / (rad_max * boa_f) ** 2

    res_sum = 0
    sub_sum = 0
    for xi in xx[0]:
        for yi in yy:
            if not ind1[yi[0], xi]:
                x_rot = int(2 * cx - xi)
                y_rot = int(2 * cy - yi[0])
                fits_f_rotated[y_rot, x_rot] = fits_f[yi[0], xi]
                res_sum += np.abs(fits_f[yi[0], xi])
                sub_sum += np.abs(fits_f[yi[0], xi] - fits_f[y_rot, x_rot])

    return fits_f_rotated, sub_sum / res_sum



for j in range(5, 6):

    fn = sorted(glob.glob(f'/Users/duhokim/work/abell/galapagos/A2670_test/t{j}/galfit/*_gf.fits'), key=os.path.getmtime, reverse=True)

    for k in cmap:
        for i in range(0, len(fn)):
        # for i in range(188, 189):
            fits_name = fn[i].split('/')[-1]
            name = fits_name.split('_')[0]
            hdu = fits.open(fn[i])
            img_orig = hdu[0].data
            img_model = hdu[1].data
            img_sub = hdu[2].data
            # img_sigma = hdu[3].data

            # img_orig[np.where(img_orig < 0)] = 0
            # img_model[np.where(img_model < 0)] = 0
            # img_sub[np.where(img_sub < 0)] = 0

            # z = cat['col18'][int(name)]/const.c*1e3                     # redshift
            # kpc_arcmin = Cosmo.kpc_proper_per_arcmin(z)         # xxx kpc / arcmin
            # arcmin_5kpc = 5.0 / kpc_arcmin                      # xxx arcmin / 5 kpc
            # frac_5kpc = arcmin_5kpc * 100.0 / img_orig.shape[0]       # calculate fraction of a length of 5 kpc in fig size

            fig = plt.figure(figsize=(12, 4))
            # fig = plt.figure(figsize=(15, 3))

            # h1 = fig.add_axes([0, 0, 0.2, 1])
            # h2 = fig.add_axes([0.2, 0, 0.2, 1])
            # h3 = fig.add_axes([0.4, 0, 0.2, 1])
            # h4 = fig.add_axes([0.6, 0, 0.2, 1])
            # h5 = fig.add_axes([0.8, 0, 0.2, 1])

            h1 = fig.add_axes([0, 0, 0.33, 1])
            h2 = fig.add_axes([0.33, 0, 0.33, 1])
            h3 = fig.add_axes([0.66, 0, 0.33, 1])

            sky = np.nanmedian(img_orig)
            ser_idx = -1
            # mv_m = np.nanmedian(img_model)
            # mv_s = np.nanmedian(img_sub)
            # mv_g = np.nanmedian(img_sigma)

            res = open(fn[i][:-4] + 'galfit.01')

            for line in res:
                if 'Sky background' in line:
                    sky = float(line.split(' ')[2])
                if 'Position x' in line:
                    pos_x = float(line.split(' ')[2])-1
                if 'Position y' in line:
                    pos_y = float(line.split(' ')[2])-1
                if 'Integrated magnitude' in line:
                    mag = float(line.split(' ')[2])
                if 'R_e' in line:
                    reff = float(line.split(' ')[2])
                if 'Sersic index' in line:
                    ser_idx = float(line.split(' ')[2])
                if 'Axis ratio' in line:
                    boa = float(line.split(' ')[2])
                if 'Position angle' in line:
                    pa = float(line.split(' ')[1])
                if 'Component number: 3' in line:
                    break


            # img_sub_rotated, ylb_ret, yhb_ret, xlb_ret, xhb_ret = rotate(img_sub, pos_x, pos_y,  -pa, boa, kron_fac * reff)
            # img_sub_rotated, a_res = rotate(img_sub, pos_x, pos_y, pa, boa, kron_fac * reff)
            # img_sub_rotated = rotate(img_sub, pos_x, pos_y, pa, boa,
            #                          img_sub.shape[0] if img_sub.shape[0] > img_sub.shape[1] else img_sub.shape[1])

            # stretch_orig = mm.jarrett(img_orig - mv_o, np.nanstd(img_orig), n)
            # stretch_model = mm.jarrett(img_model - mv_m, np.nanstd(img_model), n)
            # stretch_sub = mm.jarrett(img_sub - mv_s, np.nanstd(img_sub), n)
            # stretch_sigma = mm.jarrett(img_sigma - mv_g, np.nanstd(img_sigma), n)
            # orig_sqrt = img_scale.sqrt(img_orig - sky, scale_min=0, scale_max=0.1)
            # h1.imshow(orig_sqrt, origin='lower')
            h1.imshow(img_orig - sky, norm=LogNorm(vmin=0.1, vmax=10), origin='lower', cmap=k)
            # h1.plot([0.05, 0.05 + frac_5kpc.value], [0.02, 0.02], color='white', transform=h1.transAxes)


            # h2.imshow(img_orig - sky, norm=LogNorm(vmin=0.1, vmax=1000), origin='lower')
            h2.imshow(img_model - sky, norm=LogNorm(vmin=0.1, vmax=100), origin='lower', cmap=k)
            h2.text(0.1, 0.9, f'Sersic index = {ser_idx}', color='white', transform=h2.transAxes)
            h3.imshow(img_sub, norm=LogNorm(vmin=0.01, vmax=3), origin='lower', cmap=k)
            # h4.imshow(img_sub_rotated, norm=LogNorm(vmin=0.1, vmax=10), origin='lower')
            # h5.imshow(img_sub - img_sub_rotated, norm=LogNorm(vmin=0.1, vmax=10), origin='lower')
            # h5.text(0.1, 0.9, f'{kron_fac:.1f} * reff={kron_fac * reff:.1f}', color='black', transform=h5.transAxes)
            # h5.text(0.1, 0.8, f'A_res={a_res:.1f}', color='black', transform=h5.transAxes)
            # h5.text(0.1, 0.8, f'xmax={img_sub.shape[0]:.1f}', color='black', transform=h5.transAxes)
            # h5.text(0.1, 0.7, f'ymax={img_sub.shape[1]:.1f}', color='black', transform=h5.transAxes)

            # ell = patches.Ellipse((pos_x, pos_y),
            #                       2 * kron_fac * reff,
            #                       2 * kron_fac * reff * boa,
            #                       angle=90+pa,
            #                       linewidth=1,
            #                       fill = False
            #                       )

            # h4.add_patch(ell)


            # h1.imshow(img_orig - sky, origin='lower')
            # h2.imshow(img_orig - sky, origin='lower')
            # h3.imshow(img_model - sky, origin='lower')
            # h4.imshow(img_sub, origin='lower')
            # h1.imshow(stretch_orig, origin='lower', cmap=cmap)
            # h2.imshow(stretch_model, origin='lower', cmap=cmap)
            # h3.imshow(stretch_sub, origin='lower', cmap=cmap)
            # h4.imshow(stretch_sigma, origin='lower', cmap=cmap)


            h1.get_xaxis().set_visible(False)
            h1.get_yaxis().set_visible(False)
            h2.get_xaxis().set_visible(False)
            h2.get_yaxis().set_visible(False)
            h3.get_xaxis().set_visible(False)
            h3.get_yaxis().set_visible(False)
            # h4.get_xaxis().set_visible(False)
            # h4.get_yaxis().set_visible(False)
            # h5.get_xaxis().set_visible(False)
            # h5.get_yaxis().set_visible(False)

            results_dir = f'/Users/duhokim/work/abell/galapagos/A2670_test/t{j}/galfit/png/'

            if not os.path.isdir(results_dir):
                os.makedirs(results_dir)

            fig.savefig(results_dir + name + f'_{k}3.png')
            plt.close(fig)
