from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import abell_cluster_module as ab

tag = 'deblend'

size = 1
alpha = 0.05

params = ['MEAN_MAG_DIFF', 'DATE', 'EXP_TIME', 'SEEING', 'AIRMASS', 'TOTAL_MATCH_NUM']
num_of_param = len(params)

am_cat = pd.read_csv('/Users/duhokim/work/abell/sex/cat/airmass.csv')

fig, axs = plt.subplots(2, len(ab.clusters), tight_layout=True, figsize=(20, 6))
# fig2, axs2 = plt.subplots(2, 3, tight_layout=True, figsize=(12, 8))

for k in range(0, len(ab.clusters)):
    sex_cat = ascii.read(ab.sex_dir+'DECam_merged_SEx_cat_'+ab.clusters[k]+'_Gal_ext_corrected_'+ab.ver+
                         '_'+tag+'_merged.txt')
    # sex_cat2 = ascii.read(sex_dir+'DECam_19_21_aug_2014_SEx_cat_'+clusters[k]+'_Gal_ext_corrected_'+ver+
    #                       '_brighter_than_18mag_central_deblended_merged.txt')

    good_g = (sex_cat[ab.mag_sys+'_g'] < 30)
    good_u = (sex_cat[ab.mag_sys+'_u'] < 30)

    axs[0, k].scatter(sex_cat[ab.mag_sys+'_r'][good_g],
                      sex_cat[ab.mag_sys+'_g'][good_g] - sex_cat[ab.mag_sys+'_r'][good_g],
                      label='DECam ({:5}) \n g<{} & r<{} \n & CLASS_STAR_r<{:4.2f} & '
                            'DQM_r,cen=0,128'.format(len(sex_cat[good_g]), ab.mag_lim, ab.mag_lim,
                                                     ab.class_star_lim), alpha=alpha, s=size)
    axs[1, k].scatter(sex_cat[ab.mag_sys+'_r'][good_u],
                      sex_cat[ab.mag_sys+'_u'][good_u] - sex_cat[ab.mag_sys+'_r'][good_u],
                      label='DECam ({:4}) \n u<{} & r<{} \n & CLASS_STAR_r<{:4.2f} & '
                            'DQM_r,cen=0,128'.format(len(sex_cat[good_u]), ab.mag_lim, ab.mag_lim,
                                                     ab.class_star_lim), alpha=alpha, s=size)

    # axs[0, k].scatter(sex_cat[mag_sex+'_r'],
    #                   sex_cat[mag_sex+'_g'] - sex_cat[mag_sex+'_r'],
    #                   label='DECam ({:5}) \n g<{} & r<{} \n & CLASS_STAR_r<{:4.2f} & '
    #                         'DQM_r,cen=0,128'.format(len(sex_cat), mag_lim[1], mag_lim[2],
    #                                                  class_star_lim), alpha=alpha, s=size)
    # axs[1, k].scatter(sex_cat[mag_sex+'_r'],
    #                   sex_cat[mag_sex+'_u'] - sex_cat[mag_sex+'_r'],
    #                   label='DECam ({:4}) \n u<{} & r<{} \n & CLASS_STAR_r<{:4.2f} & '
    #                         'DQM_r,cen=0,128'.format(len(sex_cat), mag_lim[0], mag_lim[1],
    #                                                  class_star_lim), alpha=alpha, s=size)

    axs[0, k].set_title(ab.clusters[k])
    axs[0, k].set_xlim([13, 20])
    axs[1, k].set_xlim([13, 20])
    axs[0, k].set_ylim([-1, 3])
    axs[1, k].set_ylim([0, 6])
    # axs[0, k].invert_xaxis()
    # axs[1, k].invert_xaxis()
    axs[0, k].set_xlabel('r')
    axs[0, k].set_xlabel('r')
    axs[0, k].set_ylabel('g - r')
    axs[1, k].set_ylabel('u - r')

    # if k == 1 or k == 2:
    #     sdss_galaxy_all = ascii.read(cat_dir + 'sdss_galaxy_' + ab.clusters[k] + '.csv')
    #     bright_u = (sdss_galaxy_all[mag_sdss+'_u'] < mag_lim[0]) & (sdss_galaxy_all[mag_sdss+'_r'] < mag_lim[2])
    #     bright_g = (sdss_galaxy_all[mag_sdss+'_g'] < mag_lim[1]) & (sdss_galaxy_all[mag_sdss+'_r'] < mag_lim[2])
    #     sdss_galaxy_u = sdss_galaxy_all[bright_u]
    #     sdss_galaxy_g = sdss_galaxy_all[bright_g]
    #     axs[0, k].scatter(sdss_galaxy_g[mag_sdss+'_r'],
    #                       sdss_galaxy_g[mag_sdss+'_g'] - sdss_galaxy_g[mag_sdss+'_r'],
    #                       alpha = alpha,
    #                       s=size,
    #                       label = 'SDSS Galaxy Catalog ({:5}) \n g<{} & r<{}'.format(len(sdss_galaxy_all[bright_g]),
    #                       mag_lim[1], mag_lim[2]))
    #     axs[1, k].scatter(sdss_galaxy_u[mag_sdss+'_r'],
    #                       sdss_galaxy_u[mag_sdss+'_u'] - sdss_galaxy_u[mag_sdss+'_r'],
    #                       alpha = alpha,
    #                       s=size,
    #                       label = 'SDSS Galaxy Catalog ({:5}) \n u<{} & r<{}'.format(len(sdss_galaxy_all[bright_u]),
    #                       mag_lim[0], mag_lim[1]))

    # axs[0, k].legend(loc='lower left', fontsize='small')
    # axs[1, k].legend(loc='lower left', fontsize='x-small')

    # axs2[0, k].legend(loc='lower left', fontsize='small')
    # axs2[1, k].legend(loc='lower left', fontsize='x-small')
fig.savefig(ab.plot_dir + 'CMD_merged_'+ab.ver+'_'+tag+'.png')
# fig2.savefig(plot_dir + 'SDSS_cat_vs_DECam_CMD_standardized_merged_'+ver+'.png')
