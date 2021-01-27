from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
import time
import my_module
from astropy.table import Table, vstack

plot_dir=("/Users/duhokim/work/abell/plot/")
sex_dir=("/Users/duhokim/work/abell/sex/cat/")

clusters = ['A2399', 'A2670', 'A3716']

ver = 'deblended_restand'
band = ['u', 'g', 'r']

# for default SEx conf.
# a_init = [
#     [4.05, 6.06, 6.38],
#     [3.96, 6.04, 6.22],
#     [3.96, 6.15, 6.42]
# ]

# for debleding SEx conf.
# a_init = [
#     [4.177, 6.084, 6.398],
#     [3.989, 6.064, 6.234],
#     [4.008, 6.183, 6.429]
# ]

# for debleding SEx conf. and excluding neighboring standard stars
a_init = [
    [4.164, 6.084, 6.401],
    [3.987, 6.092, 6.237],
    [4.002, 6.208, 6.433]
]

max_sep = 1.0

fig, axs = plt.subplots(3, 3, tight_layout=True, figsize=(12, 12))
fig1, axs1 = plt.subplots(3, 3, tight_layout=True, figsize=(12, 12))
fig2, axs2 = plt.subplots(3, 3, tight_layout=True, figsize=(12, 12))

for k in range(0, len(clusters)):
    single_cat = ascii.read(sex_dir + 'DECam_19_21_aug_2014_single_best_exposure_SEx_cat_' + clusters[k] + '_match_rad_1as_'
                                        'Gal_ext_corrected_' + ver + '.txt')
    single_coords = SkyCoord(single_cat['ALPHA_J2000'], single_cat['DELTA_J2000'], unit='deg')

    for i in range(0, len(band)):
        start_time = time.time()
        stacked_cat = ascii.read(sex_dir + 'stack/' + clusters[k] + '_' + band[i] + 'si.cat')
        stacked_coords = SkyCoord(stacked_cat['ALPHA_J2000'], stacked_cat['DELTA_J2000'], unit='deg')

        idx, d2d, d3d = single_coords.match_to_catalog_sky(stacked_coords)
        sep_constraint = d2d.arcsec < max_sep
        single_matches = single_cat[sep_constraint]
        stacked_matches = stacked_cat[idx[sep_constraint]]

        minchi = 1e9
        a_best = 0
        a_guess = a_init[k][i] + np.linspace(-0.05, 0.05, 101)

        for j in a_guess:
            chi = 0
            for l in range(0, len(single_matches)):
                chi += np.abs(single_matches['MAG_AUTO_'+band[i]][l] - stacked_matches['MAG_AUTO'][l] - j) \
                       / np.sqrt(single_matches['MAGERR_AUTO_'+band[i]][l]**2 + stacked_matches['MAGERR_AUTO'][l]**2)
            if chi < minchi:
                minchi = chi
                a_best = j
        # a_best = a_init[k][i]

        axs[k, i].scatter(single_matches['MAG_AUTO_'+band[i]],
                          stacked_matches['MAG_AUTO'] + a_best,
                          alpha = 0.1,
                          s = 0.1)
        axs2[k, i].scatter(single_matches['MAGERR_AUTO_'+band[i]],
                          stacked_matches['MAGERR_AUTO'],
                          alpha = 0.1,
                          s = 0.1)
        axs[k, i].plot([12, 30], [12, 30], '--', alpha=0.1)
        axs2[k, i].plot([0, 0.5], [0, 0.5], '--', alpha=0.1)
        axs[k, i].set_xlim([12, 30])
        axs[k, i].set_ylim([12, 30])
        axs2[k, i].set_xlim([0, 0.5])
        axs2[k, i].set_ylim([0, 0.5])
        axs[k, i].set_title(clusters[k])
        axs1[k, i].set_title(clusters[k])
        axs2[k, i].set_title(clusters[k])
        axs[k, i].set_ylabel('stacked')
        axs[k, i].set_xlabel('single exposure')
        axs2[k, i].set_ylabel('stacked')
        axs2[k, i].set_xlabel('single exposure')
        axs[k, i].invert_xaxis()
        axs[k, i].invert_yaxis()
        axs[k, i].text(29, 13, band[i]+' = '+band[i]+'_inst (ZP=25) + {:6.4f}'.format(a_best))
        # axs[k, i].text(29, 14, r'$\chi^2$/# = {:7.3f}'.format(minchi/len(single_matches)))
        axs[k, i].text(29, 15, '# of sample = {:4d}'.format(len(single_matches)))
        # axs[k, i].legend(title='Standard Star Field (Airmass)(#)', loc='lower right', prop={'size':7})
        axs2[k, i].text(0.1, 0.4, band[i])

        axs1[k, i].hist(single_cat['MAG_AUTO_'+band[i]],
                       bins=np.arange(12, 30, 0.5),
                       label='single',
                       histtype='step')
        axs1[k, i].hist(stacked_cat['MAG_AUTO'] + a_best,
                       bins=np.arange(12.1, 30.1, 0.5),
                       label='stacked',
                       histtype='step')
        axs1[k, i].legend(loc='upper left', fontsize='small')

        print("%s %s %5.2f" % (clusters[k], band[i], a_best))
        # print(a_guess)
        print("--- %s minutes ---" % ((time.time() - start_time)/60.0))

fig.savefig(plot_dir+'DECam_stacked_vs_single_best_'+ver+'.png')
fig1.savefig(plot_dir+'DECam_stacked_vs_single_best_hist_'+ver+'.png')
fig2.savefig(plot_dir+'DECam_stacked_vs_single_best_err_'+ver+'.png')
# plt.show(block=False)

