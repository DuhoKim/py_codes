from astropy.io import ascii
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

fig, axs = plt.subplots(2, 3, tight_layout=True, figsize=(12, 8))

for k in range(0, len(clusters)):
    sex_cat = ascii.read(sex_dir+'DECam_19_21_aug_2014_single_best_exposure_SEx_cat_'+clusters[k]+'_match_rad_1as_'
                                    'Gal_ext_corrected_'+ver+'.txt')

    good_g = (sex_cat[mag_sex+'_r'] < mag_lim[2]) & (sex_cat[mag_sex+'_g'] < mag_lim[1])
    good_u = (sex_cat[mag_sex+'_r'] < mag_lim[2]) & (sex_cat[mag_sex+'_u'] < mag_lim[0])

    axs[0, k].scatter(sex_cat[mag_sex+'_r'][good_g],
                      sex_cat[mag_sex+'_g'][good_g] - sex_cat[mag_sex+'_r'][good_g],
                      label='DECam ({:5}) \n g<{} & r<{} & FWHM_r>{:4.2f}\" \n & CLASS_STAR_r<{:4.2f} & '
                            'DQM_r,cen=0,128'.format(len(sex_cat[good_g]), mag_lim[1], mag_lim[2],  fwhm_lim,
                                                     class_star_lim), alpha=0.5, s=1)
    axs[1, k].scatter(sex_cat[mag_sex+'_r'][good_u],
                      sex_cat[mag_sex+'_u'][good_u] - sex_cat[mag_sex+'_r'][good_u],
                      label='DECam ({:4}) \n u<{} & r<{} & FWHM_r>{:4.2f}\" \n & CLASS_STAR_r<{:4.2f} & '
                            'DQM_r,cen=0,128'.format(len(sex_cat[good_u]), mag_lim[0], mag_lim[1], fwhm_lim,
                                                     class_star_lim), alpha=0.5, s=1)

    axs[0, k].set_title(clusters[k])
    axs[0, k].set_xlim([13, 22])
    axs[1, k].set_xlim([13, 22])
    axs[0, k].set_ylim([-1, 2])
    axs[1, k].set_ylim([-2, 5])
    # axs[0, k].invert_xaxis()
    # axs[1, k].invert_xaxis()
    axs[0, k].set_xlabel('r')
    axs[0, k].set_xlabel('r')
    axs[0, k].set_ylabel('g - r')
    axs[1, k].set_ylabel('u - r')
    # axs[2, 0].plot([13, 21], [0, 0], linestyle=':', alpha=0.5)

    if k < 2:
        sdss_galaxy_all = ascii.read(cat_dir + 'sdss_galaxy_' + clusters[k] + '.csv')
        bright_u = (sdss_galaxy_all[mag_sdss+'_u'] < mag_lim[0]) & (sdss_galaxy_all[mag_sdss+'_r'] < mag_lim[2])
        bright_g = (sdss_galaxy_all[mag_sdss+'_g'] < mag_lim[1]) & (sdss_galaxy_all[mag_sdss+'_r'] < mag_lim[2])
        sdss_galaxy_u = sdss_galaxy_all[bright_u]
        sdss_galaxy_g = sdss_galaxy_all[bright_g]
        axs[0, k].scatter(sdss_galaxy_g[mag_sdss+'_r'],
                          sdss_galaxy_g[mag_sdss+'_g'] - sdss_galaxy_g[mag_sdss+'_r'],
                          alpha = 0.1,
                          s=1,
                          label = 'SDSS Galaxy Catalog ({:5}) \n g<{} & r<{}'.format(len(sdss_galaxy_all[bright_g]),
                          mag_lim[1], mag_lim[2]))
        axs[1, k].scatter(sdss_galaxy_u[mag_sdss+'_r'],
                          sdss_galaxy_u[mag_sdss+'_u'] - sdss_galaxy_u[mag_sdss+'_r'],
                          alpha = 0.1,
                          s=1,
                          label = 'SDSS Galaxy Catalog ({:5}) \n u<{} & r<{}'.format(len(sdss_galaxy_all[bright_u]),
                          mag_lim[0], mag_lim[1]))

    axs[0, k].legend(loc='lower left', fontsize='small')
    axs[1, k].legend(loc='lower left', fontsize='x-small')
fig.savefig(plot_dir + 'SDSS_cat_vs_DECam_CMD_standardized_best_single_exposure_'+ver+'.png')
