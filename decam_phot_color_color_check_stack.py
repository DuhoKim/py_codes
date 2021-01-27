from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import abell_cluster_module as ab

size = 1
alpha = 0.05

tag = 'deblend'

params = ['MEAN_MAG_DIFF', 'DATE', 'EXP_TIME', 'SEEING', 'AIRMASS', 'TOTAL_MATCH_NUM']
num_of_param = len(params)

am_cat = pd.read_csv('/Users/duhokim/work/abell/sex/cat/airmass.csv')

fig, axs = plt.subplots(2, len(ab.clusters), tight_layout=True, figsize=(20, 6))
# fig2, axs2 = plt.subplots(2, 3, tight_layout=True, figsize=(12, 8))

for k in range(0, len(ab.clusters)):
    sex_cat = ascii.read(ab.sex_dir+'DECam_merged_SEx_cat_'+ab.clusters[k]+'_'+
                                    'Gal_ext_corrected_'+ab.ver+'_central_deblended_'+tag+'_stack.txt')

    good_g = (sex_cat['col6'] < 30)
    good_u = (sex_cat['col4'] < 30)

    axs[0, k].scatter(sex_cat['col8'][good_g],
                      sex_cat['col6'][good_g] - sex_cat['col8'][good_g],
                      label='DECam ({:5}) \n g<{} & r<{} \n & CLASS_STAR_r<{:4.2f} '
                            .format(len(sex_cat[good_g]), ab.mag_lim, ab.mag_lim,
                                                     ab.class_star_lim), alpha=alpha, s=size)
    axs[1, k].scatter(sex_cat['col8'][good_u],
                      sex_cat['col4'][good_u] - sex_cat['col8'][good_u],
                      label='DECam ({:4}) \n u<{} & r<{} \n & CLASS_STAR_r<{:4.2f} '
                            .format(len(sex_cat[good_u]), ab.mag_lim, ab.mag_lim,
                                                     ab.class_star_lim), alpha=alpha, s=size)

    axs[0, k].set_title(ab.clusters[k])
    axs[0, k].set_xlabel('r')
    axs[0, k].set_xlabel('r')
    axs[0, k].set_ylabel('g - r')
    axs[1, k].set_ylabel('u - r')
    axs[0, k].set_xlim([13, 20])
    axs[1, k].set_xlim([13, 20])
    axs[0, k].set_ylim([-1, 3])
    axs[1, k].set_ylim([0, 6])

fig.savefig(ab.plot_dir + 'CMD_stacked_'+ab.ver+'_'+tag+'.png')

