from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
import time
import my_module
from astropy.table import Table, vstack
import statsmodels.api as sm
from statsmodels.sandbox.regression.predstd import wls_prediction_std
import pandas as pd
from scipy.stats import sigmaclip

plot_dir=("/Users/duhokim/work/abell/plot/")
sex_dir=("/Users/duhokim/work/abell/sex/cat/")

clusters = ['A754', 'A2399', 'A2670', 'A3558', 'A3574', 'A3659', 'A3716']

ver = '2013-2014'
band = ['u', 'g', 'r']

max_sep = 1.0

sig = 3.0

fig, axs = plt.subplots(len(clusters), 3, tight_layout=True, figsize=(10, 14))
fig1, axs1 = plt.subplots(len(clusters), 3, tight_layout=True, figsize=(10, 14))
fig2, axs2 = plt.subplots(len(clusters), 3, tight_layout=True, figsize=(10, 14))

for k in range(0, len(clusters)):
    single_cat = ascii.read(sex_dir + 'DECam_short_exposure_SEx_cat_' + clusters[k] + '_match_rad_1as_'
                                        'Gal_ext_corrected_' + ver + '.txt')
    single_coords = SkyCoord(single_cat['ALPHA_J2000'], single_cat['DELTA_J2000'], unit='deg')

    for i in range(0, len(band)):
        start_time = time.time()
        stacked_all_cat = ascii.read(sex_dir + 'stack/' + clusters[k] + '_' + band[i] + 'si.cat')
        stacked_cat = stacked_all_cat[stacked_all_cat['MAG_AUTO'] < 90]
        stacked_coords = SkyCoord(stacked_cat['ALPHA_J2000'], stacked_cat['DELTA_J2000'], unit='deg')

        idx, d2d, d3d = single_coords.match_to_catalog_sky(stacked_coords)
        sep_constraint = (d2d.arcsec < max_sep) & (single_cat['MAG_AUTO_'+band[i]] < 90)
        single_matches = single_cat[sep_constraint]
        stacked_matches = stacked_cat[idx[sep_constraint]]

        c, low, upp = sigmaclip(single_matches['MAG_AUTO_'+band[i]] - stacked_matches['MAG_AUTO'], sig, sig)

        ind = np.where((single_matches['MAG_AUTO_'+band[i]] - stacked_matches['MAG_AUTO'] > low) &
                       (single_matches['MAG_AUTO_'+band[i]] - stacked_matches['MAG_AUTO'] < upp))
        ind2 = np.where((single_matches['MAG_AUTO_'+band[i]] - stacked_matches['MAG_AUTO'] < low) |
                        (single_matches['MAG_AUTO_'+band[i]] - stacked_matches['MAG_AUTO'] > upp))

        Data_Frame = {'delm': single_matches['MAG_AUTO_'+band[i]][ind] - stacked_matches['MAG_AUTO'][ind],
                      'am': np.zeros(len(ind[0])),
                      'err': single_matches['MAGERR_AUTO_'+band[i]][ind]**2 + stacked_matches['MAGERR_AUTO'][ind]**2
                      }

        df = pd.DataFrame(Data_Frame, columns=['delm', 'am', 'err'])

        X = df[['am']]
        Y = df['delm']

        # with statsmodels
        X = sm.add_constant(X)

        model = sm.WLS(Y, X, weights=1 / df['err']).fit()
        # results = model.fit()
        prstd, iv_l, iv_u = wls_prediction_std(model)

        print("{} {} {} {} {}".format(clusters[k], band[i], model.params[0], model.bse[0], len(ind[0])))


        axs[k, i].scatter(single_matches['MAG_AUTO_'+band[i]][ind],
                          single_matches['MAG_AUTO_'+band[i]][ind] - stacked_matches['MAG_AUTO'][ind],
                          alpha = 0.1,
                          s = 0.1)
        axs[k, i].scatter(single_matches['MAG_AUTO_'+band[i]][ind2],
                          single_matches['MAG_AUTO_'+band[i]][ind2] - stacked_matches['MAG_AUTO'][ind2],
                          alpha = 0.1,
                          s = 0.1,
                          color = 'gray')
        axs2[k, i].scatter(single_matches['MAGERR_AUTO_'+band[i]],
                           stacked_matches['MAGERR_AUTO'],
                           alpha = 0.1,
                           s = 0.1)
        axs[k, i].plot([12, 30], [low, low], '--', alpha=0.1, color='salmon')
        axs[k, i].plot([12, 30], [upp, upp], '--', alpha=0.1, color='salmon')
        axs[k, i].plot([12, 30], [model.params[0], model.params[0]], '--', alpha=0.2)
        axs[k, i].fill_between([12, 30],
                         [model.params[0] + model.bse[0], model.params[0] + model.bse[0]],
                         [model.params[0] - model.bse[0], model.params[0] - model.bse[0]],
                         alpha=0.1,
                         color='blue'
                         )

        axs2[k, i].plot([0, 0.5], [0, 0.5], '--', alpha=0.1)
        axs[k, i].set_xlim([12, 30])
        # axs[k, i].set_ylim([12, 30])
        axs2[k, i].set_xlim([0, 0.5])
        axs2[k, i].set_ylim([0, 0.5])
        axs[k, i].set_title(clusters[k])
        axs1[k, i].set_title(clusters[k])
        axs2[k, i].set_title(clusters[k])
        if i == 0:
            axs[k, i].set_ylabel(r'$\Delta m$ ($m_{single} - m_{stacked}$)')
        if k == 2:
            axs[k, i].set_xlabel(r'$m_{single}$')
        axs2[k, i].set_ylabel('stacked')
        axs2[k, i].set_xlabel('single exposure')
        # axs[k, i].invert_xaxis()
        # axs[k, i].invert_yaxis()
        axs[k, i].text(0.1, 0.9,
                       band[i]+' = '+band[i]+'_inst (ZP=25) + {:6.4f}'.format(model.params[0]),
                       transform = axs[k, i].transAxes)
        # axs[k, i].text(29, 14, r'$\chi^2$/# = {:7.3f}'.format(minchi/len(single_matches)))
        axs[k, i].text(0.1, 0.8,
                       '# of sample = {:4d}'.format(len(ind[0])),
                       transform = axs[k, i].transAxes)
        # axs[k, i].legend(title='Standard Star Field (Airmass)(#)', loc='lower right', prop={'size':7})
        axs2[k, i].text(0.1, 0.4, band[i])

        axs1[k, i].hist(single_cat['MAG_AUTO_'+band[i]],
                       bins=np.arange(12, 30, 0.5),
                       label='single',
                       histtype='step')
        axs1[k, i].hist(stacked_cat['MAG_AUTO'] + model.params[0],
                       bins=np.arange(12.1, 30.1, 0.5),
                       label='stacked',
                       histtype='step')
        axs1[k, i].legend(loc='upper left', fontsize='small')

        print("%s %s %5.2f" % (clusters[k], band[i], model.params[0]))
        # print(a_guess)
        print("--- %s minutes ---" % ((time.time() - start_time)/60.0))

fig.savefig(plot_dir+'DECam_stacked_vs_single_best_'+ver+'.png')
fig1.savefig(plot_dir+'DECam_stacked_vs_single_best_hist_'+ver+'.png')
fig2.savefig(plot_dir+'DECam_stacked_vs_single_best_err_'+ver+'.png')
# plt.show(block=False)

