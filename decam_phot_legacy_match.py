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
import abell_cluster_module as ab
import importlib
importlib.reload(ab)

plot_dir=("/Users/duhokim/work/abell/plot/")
sex_dir=("/Users/duhokim/work/abell/sex/cat/")
cat_dir=("/Users/duhokim/work/abell/cat/")

merging_decals = True

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

fig = plt.figure(figsize=(8, 10))

gs = fig.add_gridspec(3, 2, width_ratios = (5, 5), height_ratios=(4, 2 ,4))

for k in range(0, len(clusters)):
    if not merging_decals:
        sex_cat = ascii.read(ab.sex_dir+f'DECam_merged_SEx_cat_{clusters[k]}_Gal_ext_corrected_{ab.ver}.txt')
        sex_coords = SkyCoord(sex_cat['ALPHA_J2000'], sex_cat['DELTA_J2000'], unit='deg')
    else:
        leg_total = Table()

    is_first = True

    for file in os.listdir(cat_dir + 'legacy/a2670'):
        if file.endswith('.fits'):
            data = Table.read(cat_dir + 'legacy/a2670/' + file)
            # hdu1 = fits.open(cat_dir + 'legacy/a2670/' + file)
            # data = hdu1[1].data

            if not merging_decals:
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
            else:
                leg_total = vstack([leg_total, data])
    leg_total.write(f'/Users/duhokim/work/abell/catalogue/a2670_legacy.cat', format='ascii', overwrite=True)

    if not merging_decals:

        is_good = (sex_total[mag_sex+'_g'] < 900) & (leg_total['flux_g'] > 0)
        ax00 = fig.add_subplot(gs[0, 0])
        ax00.scatter(sex_total[mag_sex+'_g'][is_good],
                       22.5 - 2.5 * np.log10(leg_total['flux_g'][is_good]),
                       alpha = 0.01,
                       s = 1)
        ax10 = fig.add_subplot(gs[1, 0], sharex=ax00)
        ax10.scatter(sex_total[mag_sex+'_g'][is_good],
                          sex_total[mag_sex + '_g'][is_good] - (22.5 - 2.5 * np.log10(leg_total['flux_g'][is_good])),
                          alpha=0.01,
                          s=1)
        ax20 = fig.add_subplot(gs[2, 0])
        ax20.hist(sex_total[mag_sex+'_g'][is_good],
                       bins=np.arange(15, 26, 0.5),
                       label='Our catalog',
                       histtype='step')
        ax20.hist(22.5 - 2.5 * np.log10(leg_total['flux_g'][is_good]),
                       bins=np.arange(15.2, 26, 0.5),
                       label='DECaLS ~140s \n exposure(s) (Tractor)',
                       histtype='step')

        is_good = (sex_total[mag_sex+'_r'] < 900) & (leg_total['flux_r'] > 0)
        ax01 = fig.add_subplot(gs[0, 1])
        ax01.scatter(sex_total[mag_sex+'_r'][is_good],
                       22.5 - 2.5 * np.log10(leg_total['flux_r'][is_good]),
                       alpha = 0.01,
                       s = 1)
        ax11 = fig.add_subplot(gs[1, 1], sharex=ax01)
        ax11.scatter(sex_total[mag_sex+'_r'][is_good],
                       sex_total[mag_sex+'_r'][is_good] - (22.5 - 2.5 * np.log10(leg_total['flux_r'][is_good])),
                       alpha = 0.01,
                       s = 1)
        ax21 = fig.add_subplot(gs[2, 1])
        ax21.hist(sex_total[mag_sex+'_r'][is_good],
                       bins=np.arange(15, 26, 0.5),
                       label='Our catalog',
                       histtype='step')
        ax21.hist(22.5 - 2.5 * np.log10(leg_total['flux_r'][is_good]),
                       bins=np.arange(15.1, 26, 0.5),
                       label='DECaLS ~100s \n exposure(s) (Tractor)',
                       histtype='step')

        # axs[0, 0].text(25, 13, 'Match to DECaLS (A2670)', fontsize=24)
        ax00.text(25, 16, 'g', fontsize=24)
        ax01.text(25, 16, 'r', fontsize=24)
        # axs[0, 0].set_title('g band')
        # axs[0, 1].set_title('r band')
        ax00.set_xlim([mag_lim_faint, mag_lim_bright])
        ax01.set_xlim([mag_lim_faint, mag_lim_bright])
        ax00.set_ylim([mag_lim_faint, mag_lim_bright])
        ax01.set_ylim([mag_lim_faint, mag_lim_bright])
        # axs[0].invert_xaxis()
        # axs[1].invert_xaxis()
        ax10.set_xlabel('Our catalog')
        ax11.set_xlabel('Our catalog')
        ax10.set_ylim([-0.5, 0.5])
        ax11.set_ylim([-0.5, 0.5])
        ax10.set_ylabel('Ours - DECaLS (mag)')

        ax00.set_ylabel('DECaLS Tractor (~140s)')
        ax01.set_ylabel('DECaLS Tractor (~100s)')
        ax00.plot([mag_lim_faint, mag_lim_bright], [mag_lim_faint, mag_lim_bright], linestyle=':', alpha=0.5)
        ax01.plot([mag_lim_faint, mag_lim_bright], [mag_lim_faint, mag_lim_bright], linestyle=':', alpha=0.5)

        ax00.set_aspect('equal', 'box')
        ax01.set_aspect('equal', 'box')

        ax20.legend(loc='upper left', fontsize='small')
        ax21.legend(loc='upper left', fontsize='small')

        ax20.set_xlabel('g')
        ax21.set_xlabel('r')
        ax20.set_ylabel('N')
        ax21.set_ylabel('N')
if not merging_decals:
    fig.savefig(plot_dir + 'DECaL_vs_DECam_merged_A2670.png')
