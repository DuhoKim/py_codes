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

work_dir = '/Users/duhokim/work/abell/galapagos/A3558_rsr_all/t1/galfit/'
cat = ascii.read("/Users/duhokim/work/abell/spec/Shapley/catalog.dat")

n = 7
cmap = 'gray'

for j in range(1, 10):

    fn = sorted(glob.glob(f'/Users/duhokim/work/abell/galapagos/A3558_rsr_cat_all/t{j}/galfit/*_gf.fits'), key=os.path.getmtime, reverse=True)

    for i in range(0, len(fn)):
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

        # h1 = fig.add_axes([0, 0, 0.25, 1])
        # h2 = fig.add_axes([0.25, 0, 0.25, 1])
        # h3 = fig.add_axes([0.5, 0, 0.25, 1])
        # h4 = fig.add_axes([0.75, 0, 0.25, 1])

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
                line_spl = line.split(' ')
                sky = float(line_spl[2])
            if 'Sersic index' in line:
                line_spl = line.split(' ')
                ser_idx = float(line_spl[2])
                break

        # stretch_orig = mm.jarrett(img_orig - mv_o, np.nanstd(img_orig), n)
        # stretch_model = mm.jarrett(img_model - mv_m, np.nanstd(img_model), n)
        # stretch_sub = mm.jarrett(img_sub - mv_s, np.nanstd(img_sub), n)
        # stretch_sigma = mm.jarrett(img_sigma - mv_g, np.nanstd(img_sigma), n)

        h1.imshow(img_orig - sky, norm=LogNorm(vmin=0.1, vmax=100), origin='lower')
        # h1.plot([0.05, 0.05 + frac_5kpc.value], [0.02, 0.02], color='white', transform=h1.transAxes)


        # h2.imshow(img_orig - sky, norm=LogNorm(vmin=0.1, vmax=1000), origin='lower')
        h2.imshow(img_model - sky, norm=LogNorm(vmin=0.1, vmax=100), origin='lower')
        h2.text(0.1, 0.9, f'Sersic index = {ser_idx}', color='white', transform=h2.transAxes)
        h3.imshow(img_sub, norm=LogNorm(vmin=0.1, vmax=10), origin='lower')
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

        results_dir = f'/Users/duhokim/work/abell/galapagos/A3558_rsr_cat_all/t{j}/galfit/png2/'

        if not os.path.isdir(results_dir):
            os.makedirs(results_dir)

        fig.savefig(results_dir + name + '.png')
        plt.close(fig)
