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
import importlib
importlib.reload(img_scale)
from astropy.coordinates import SkyCoord

work_dir = '/Users/duhokim/work/abell/galapagos/A3558_rsr_all/t1/galfit/'
cat = ascii.read("/Users/duhokim/work/abell/spec/Shapley/catalog.dat")
sex_cat = ascii.read("/Users/duhokim/work/abell/sex/background/A2670_back.cat")

coords_sex = SkyCoord(sex_cat['ALPHA_J2000'], sex_cat['DELTA_J2000'], unit='deg')

n = 7
cmap = 'gray'
kron_fac = 2.2
rot_box = 1.5

with open('/Users/duhokim/work/abell/galapagos/A2670_test/A2670_back_values.txt', 'w') as rec:
    rec.writelines(f'# id ra dec galapagos_sky nanmedian sig_med_clip sig_mean_clip sex \n')

    for j in range(1, 10):

        fn = sorted(glob.glob(f'/Users/duhokim/work/abell/galapagos/A2670_test/t{j}/galfit/*_gf.fits'), key=os.path.getmtime, reverse=True)
        out_cat = ascii.read(f'/Users/duhokim/work/abell/galapagos/A2670_test/t{j}/t{j}r.outcat')

        for i in range(0, len(fn)):
        # for i in range(188, 189):
            fits_name = fn[i].split('/')[-1]
            name = fits_name.split('_')[0]
            ind = np.where(out_cat['col1'] == int(name[2:]))
            coord_gal = SkyCoord(out_cat['col13'][ind], out_cat['col14'][ind], unit='deg')

            hdu = fits.open(fn[i])
            img_orig = hdu[0].data
            img_model = hdu[1].data
            img_sub = hdu[2].data

            fig = plt.figure(figsize=(15, 3))

            h1 = fig.add_axes([0, 0, 0.2, 1])
            h2 = fig.add_axes([0.2, 0, 0.2, 1])
            h3 = fig.add_axes([0.4, 0, 0.2, 1])
            h4 = fig.add_axes([0.6, 0, 0.2, 1])
            h5 = fig.add_axes([0.8, 0, 0.2, 1])

            sky = np.nanmedian(img_orig)
            ser_idx = -1

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

            orig_sqrt = img_scale.sqrt(img_orig - sky, scale_min=0, scale_max=10)
            h1.imshow(orig_sqrt, origin='lower', cmap='rainbow')
            h1.text(0.1, 0.9, f'sky ={sky:.3f}', color='black', transform=h1.transAxes)
            orig_asinh = img_scale.asinh(img_orig - sky, scale_min=0, scale_max=10)
            h2.imshow(orig_asinh, origin='lower', cmap='rainbow')
            h2.text(0.1, 0.9, f'nanmedian ={np.nanmedian(img_orig):.3f}', color='black', transform=h2.transAxes)
            orig_log = img_scale.log(img_orig - sky, scale_min=0, scale_max=10)
            h3.imshow(orig_log, origin='lower', cmap='rainbow')

            sky_med_sig_clip1, num_med_iter1 = img_scale.sky_median_sig_clip(img_orig, 1, 0.0000001, max_iter=100)
            sky_med_sig_clip2, num_med_iter2 = img_scale.sky_median_sig_clip(img_orig, 2, 0.0000001, max_iter=100)
            sky_med_sig_clip3, num_med_iter3 = img_scale.sky_median_sig_clip(img_orig, 3, 0.0000001, max_iter=100)
            h3.text(0.1, 0.9, f'sky_median_sig_clip1 ={sky_med_sig_clip1:.3f}', color='black', transform=h3.transAxes)
            h3.text(0.1, 0.8, f'num_med_iter1 ={num_med_iter1}', color='black', transform=h3.transAxes)
            h3.text(0.1, 0.7, f'sky_median_sig_clip2 ={sky_med_sig_clip2:.3f}', color='black', transform=h3.transAxes)
            h3.text(0.1, 0.6, f'num_med_iter2 ={num_med_iter2}', color='black', transform=h3.transAxes)
            h3.text(0.1, 0.5, f'sky_median_sig_clip3 ={sky_med_sig_clip3:.3f}', color='black', transform=h3.transAxes)
            h3.text(0.1, 0.4, f'num_med_iter3 ={num_med_iter3}', color='black', transform=h3.transAxes)

            orig_linear = img_scale.linear(img_orig - sky, scale_min=0, scale_max=10)
            h4.imshow(orig_linear, origin='lower', cmap='rainbow')

            sky_mean_sig_clip1, num_mean_iter1 = img_scale.sky_mean_sig_clip(img_orig, 1, 0.000001)
            sky_mean_sig_clip2, num_mean_iter2 = img_scale.sky_mean_sig_clip(img_orig, 2, 0.000001)
            sky_mean_sig_clip3, num_mean_iter3 = img_scale.sky_mean_sig_clip(img_orig, 3, 0.000001)
            h4.text(0.1, 0.9, f'sky_mean_sig_clip1 ={sky_med_sig_clip1:.3f}', color='black', transform=h4.transAxes)
            h4.text(0.1, 0.8, f'num_mean_iter1 ={num_mean_iter1}', color='black', transform=h4.transAxes)
            h4.text(0.1, 0.7, f'sky_mean_sig_clip2 ={sky_med_sig_clip2:.3f}', color='black', transform=h4.transAxes)
            h4.text(0.1, 0.6, f'num_mean_iter2 ={num_mean_iter2}', color='black', transform=h4.transAxes)
            h4.text(0.1, 0.5, f'sky_mean_sig_clip3 ={sky_med_sig_clip3:.3f}', color='black', transform=h4.transAxes)
            h4.text(0.1, 0.4, f'num_mean_iter3 ={num_mean_iter3}', color='black', transform=h4.transAxes)

            d2d = coord_gal.separation(coords_sex)
            matched_sex = (d2d.arcsec < 1)
            if sum(matched_sex):
                back_sex = sex_cat['BACKGROUND'][matched_sex][0]
            else:
                back_sex = 0.0
            h5.imshow(orig_sqrt, origin='lower', cmap='winter')
            h5.text(0.1, 0.9, f'sky_sex ={back_sex:.3f}', color='black', transform=h5.transAxes)

            rec.writelines(f'{name} {coord_gal[0].ra.value} {coord_gal[0].dec.value} {sky} {np.nanmedian(img_orig):.3f} '
                           f'{sky_med_sig_clip2:.3f} {sky_med_sig_clip2:.3f} {back_sex:.3f} \n')


            h1.get_xaxis().set_visible(False)
            h1.get_yaxis().set_visible(False)
            h2.get_xaxis().set_visible(False)
            h2.get_yaxis().set_visible(False)
            h3.get_xaxis().set_visible(False)
            h3.get_yaxis().set_visible(False)
            h4.get_xaxis().set_visible(False)
            h4.get_yaxis().set_visible(False)
            h5.get_xaxis().set_visible(False)
            h5.get_yaxis().set_visible(False)

            results_dir = f'/Users/duhokim/work/abell/galapagos/A2670_test/t{j}/galfit/png/'

            if not os.path.isdir(results_dir):
                os.makedirs(results_dir)

            # fig.savefig(results_dir + name + '.png')
            plt.close(fig)
