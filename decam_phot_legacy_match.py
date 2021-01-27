from astropy.io import ascii, fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy.coordinates import FK4, Angle, Latitude, Longitude
from astropy import coordinates as coords
import os
from os import listdir
from os.path import isfile, join
import astropy.units as u
from astropy.table import Table, vstack
import pandas as pd
from pandas import DataFrame, read_csv
# import matplotlib
# matplotlib.rcParams['text.usetex'] = True
from matplotlib.patches import Circle
import time
import my_module

plot_dir=("/Users/duhokim/work/abell/plot/")
sex_dir=("/Users/duhokim/work/abell/sex/cat/")
cat_dir=("/Users/duhokim/work/abell/cat/")

ver = 'v1.1'
fwhm_lim = 0.0 * 0.263
class_star_lim = 0.9
# mag_lim = [21, 22, 22]
mag_lim = [99, 99, 99]
mag_sex = 'MAG_AUTO'
mag_sdss = 'petroMag'
mag_lim_faint = 26
mag_lim_bright = 14

clusters = ['A2670']

bands = ['u', 'g', 'r']
date_col = ['red', 'green', 'blue'] # for 19aug, 20aug, 21aug
date_lab = ['19aug', '20aug', '21aug']
# for_gal_ext_u = [0.159, 0.188, 0.157] # irsa.ipac.caltech.edu, S and F (2011)
# for_gal_ext_g = [0.124, 0.146, 0.122]
# for_gal_ext_r = [0.086, 0.101, 0.085]

for_gal_ext_u = [0, 0, 0] # SDSS Galaxy cat didn't correct for Galactic extinction?
for_gal_ext_g = [0, 0, 0]
for_gal_ext_r = [0, 0, 0]

max_sep = 1.0                   # matching radius limit among each individual exposures

params = ['MEAN_MAG_DIFF', 'DATE', 'EXP_TIME', 'SEEING', 'AIRMASS', 'TOTAL_MATCH_NUM']
num_of_param = len(params)

am_cat = pd.read_csv('/Users/duhokim/work/abell/sex/cat/airmass.csv')

fig, axs = plt.subplots(2, 2, tight_layout=True, figsize=(8, 8))

for k in range(0, len(clusters)):
    sex_cat = ascii.read(sex_dir+'DECam_19_21_aug_2014_single_best_exposure_SEx_cat_'+clusters[k]+'_match_rad_1as_'
                                    'Gal_ext_corrected_'+ver+'.txt')
    sex_coords = SkyCoord(sex_cat['ALPHA_J2000'], sex_cat['DELTA_J2000'], unit='deg')

    is_first = True

    for file in os.listdir(cat_dir + 'legacy/a2670'):
        if file.endswith('.fits'):
            data = Table.read(cat_dir + 'legacy/a2670/' + file)
            # hdu1 = fits.open(cat_dir + 'legacy/a2670/' + file)
            # data = hdu1[1].data

            leg_coords = SkyCoord(data['ra'], data['dec'], unit='deg')

            idx, d2d, d3d = sex_coords.match_to_catalog_sky(leg_coords)
            sep_constraint = d2d.arcsec < max_sep
            sex_matches = sex_cat[sep_constraint]
            leg_matches = data[idx[sep_constraint]]

            if is_first:
                sex_total = sex_matches
                leg_total = leg_matches
                is_first = False
            else:
                sex_total = vstack([sex_total, sex_matches])
                leg_total = vstack([leg_total, leg_matches])

    is_good = (sex_total[mag_sex+'_g'] < 900) & (leg_total['flux_g'] > 0)
    axs[0, 0].scatter(sex_total[mag_sex+'_g'][is_good],
                   22.5 - 2.5 * np.log10(leg_total['flux_g'][is_good]),
                   alpha = 0.01,
                   s = 1)
    axs[1, 0].hist(sex_total[mag_sex+'_g'][is_good],
                   bins=np.arange(15, 26, 0.5),
                   label='Our 60s single exposure (SEx)',
                   histtype='step')
    axs[1, 0].hist(22.5 - 2.5 * np.log10(leg_total['flux_g'][is_good]),
                   bins=np.arange(15.2, 26, 0.5),
                   label='DECaLS ~140s exposure(s) (Tractor)',
                   histtype='step')

    is_good = (sex_total[mag_sex+'_r'] < 900) & (leg_total['flux_r'] > 0)
    axs[0, 1].scatter(sex_total[mag_sex+'_r'][is_good],
                   22.5 - 2.5 * np.log10(leg_total['flux_r'][is_good]),
                   alpha = 0.01,
                   s = 1)
    axs[1, 1].hist(sex_total[mag_sex+'_r'][is_good],
                   bins=np.arange(15, 26, 0.5),
                   label='Our 300s single exposure (SEx)',
                   histtype='step')
    axs[1, 1].hist(22.5 - 2.5 * np.log10(leg_total['flux_r'][is_good]),
                   bins=np.arange(15.2, 26, 0.5),
                   label='DECaLS ~100s exposure(s) (Tractor)',
                   histtype='step')

    axs[0, 0].text(15, 13, 'A2670', fontsize=24)
    axs[0, 0].set_title('g band')
    axs[0, 1].set_title('r band')
    axs[0, 0].set_xlim([mag_lim_faint, mag_lim_bright])
    axs[0, 1].set_xlim([mag_lim_faint, mag_lim_bright])
    axs[0, 0].set_ylim([mag_lim_faint, mag_lim_bright])
    axs[0, 1].set_ylim([mag_lim_faint, mag_lim_bright])
    # axs[0].invert_xaxis()
    # axs[1].invert_xaxis()
    axs[0, 0].set_xlabel('SEx MAG_AUTO (60s)')
    axs[0, 1].set_xlabel('SEx MAG_AUTO (300s)')
    axs[0, 0].set_ylabel('DECaLS Tractor (~140s)')
    axs[0, 1].set_ylabel('DECaLS Tractor (~100s)')
    axs[0, 0].plot([mag_lim_faint, mag_lim_bright], [mag_lim_faint, mag_lim_bright], linestyle=':', alpha=0.5)
    axs[0, 1].plot([mag_lim_faint, mag_lim_bright], [mag_lim_faint, mag_lim_bright], linestyle=':', alpha=0.5)

    axs[0, 0].set_aspect('equal', 'box')
    axs[0, 1].set_aspect('equal', 'box')

    axs[1, 0].legend(loc='upper left', fontsize='small')
    axs[1, 1].legend(loc='upper left', fontsize='small')
fig.savefig(plot_dir + 'DECaL_vs_DECam_A2670.png')
