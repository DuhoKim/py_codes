from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
import time
import my_module
from astropy.io import fits
from astropy.table import Table, vstack
import statsmodels.api as sm
from statsmodels.sandbox.regression.predstd import wls_prediction_std
import pandas as pd
from scipy.stats import sigmaclip
import abell_cluster_module as ab
import glob
from multiprocessing import Pool
from astropy import wcs
import importlib
importlib.reload(ab)

work_dir=("/Users/duhokim/work/abell/")
plot_dir=("/Users/duhokim/work/abell/plot/")
sex_dir=("/Users/duhokim/work/abell/sex/cat/")

use_saved = True

fs = 20     # fontsize

max_sep = 1.0
sig = 3.0

ccd_x = 2046
ccd_y = 4094

margin = 200

def dqm_check(n):
    prev_time = start_time
    is_dqm_zero = np.ones(ni, dtype=bool)
    hdu_r_single = fits.open(ab.work_dir + 'fits/best_single/' + ab.short_dqm_fn[k][i])
    hdu_r_stack = fits.open(ab.work_dir + 'fits/stacked/' + ab.stack_dqm_fn[k][i])

    for ii in range(0, ni):
        # read the celestial coordinate
        cel_coord = [[stacked_matches['ALPHA_J2000'][ii+n*ni],
                      stacked_matches['DELTA_J2000'][ii+n*ni]], [0, 0]]

        for jj in range(1, len(hdu_r_single)):
            # read WCS
            w = wcs.WCS(hdu_r_single[jj].header)
            pixcrd = w.wcs_world2pix(cel_coord, 1)
            # if w.footprint_contains(sky_coord):
            if (pixcrd[0][0] > ab.dqm_boundary) & (pixcrd[0][0] < hdu_r_single[jj].shape[1] - ab.dqm_boundary) & \
                    (pixcrd[0][1] > ab.dqm_boundary) & (pixcrd[0][1] < hdu_r_single[jj].shape[0] - ab.dqm_boundary):
                dqm_fits = hdu_r_single[jj].data
                # check the value of DQM
                if k == 0:      # Somehow A754 has value 7 (cosmic ray) in galaxies
                    if dqm_fits[int(pixcrd[0][1])][int(pixcrd[0][0])] % 7:
                        # toggle on saturated bool array
                        is_dqm_zero[ii] = False
                else:           # somehow other clusters have value 128 (?) in galaxies
                    if dqm_fits[int(pixcrd[0][1])][int(pixcrd[0][0])] % 128:
                        # toggle on saturated bool array
                        is_dqm_zero[ii] = False

        for jj in range(1, len(hdu_r_stack)):
            # read WCS
            w = wcs.WCS(hdu_r_stack[jj].header)
            pixcrd = w.wcs_world2pix(cel_coord, 1)
            # if w.footprint_contains(sky_coord):
            if (pixcrd[0][0] > 0) & (pixcrd[0][0] < hdu_r_stack[jj].shape[1]) & \
                    (pixcrd[0][1] > 0) & (pixcrd[0][1] < hdu_r_stack[jj].shape[0]):
                dqm_fits = hdu_r_stack[jj].data
                # check the value of DQM
                if dqm_fits[int(pixcrd[0][1])][int(pixcrd[0][0])]:
                    # toggle on saturated bool array
                    is_dqm_zero[ii] = False

        if ii % 1000 == 0:
            print("--- %06.2f minutes ---" % ((time.time() - prev_time) / 60.0))
            print("--- %d / %d - %d (DQM check) ---" % (ii, n*ni, (n+1)*ni))
            prev_time = time.time()

    hdu_r_single.close()
    hdu_r_stack.close()
    return is_dqm_zero

fig, axs = plt.subplots(len(ab.clusters), 3, tight_layout=True, figsize=(10, 14))
# fig11, axs11 = plt.subplots(1, 3, tight_layout=True, figsize=(10, 3))
fig11, axs11 = plt.subplots(1, 3, figsize=(10, 3))
fig1, axs1 = plt.subplots(len(ab.clusters), 3, tight_layout=True, figsize=(10, 14))
fig2, axs2 = plt.subplots(len(ab.clusters), 3, tight_layout=True, figsize=(10, 14))
fig3, axs3 = plt.subplots(len(ab.clusters), 3, tight_layout=True, figsize=(10, 14))

for k in range(0, len(ab.clusters)):
# for k in range(0, 1):
    for i in range(0, len(ab.bands)):
        if not use_saved:
            single_cat_fns = glob.glob(f'{work_dir}sex/run_short/{ab.clusters[k]}_{ab.bands[i]}si_*_psf.cat')
            single_all_cat = Table()
            for single_cat_fn in single_cat_fns:
                single_cat = ascii.read(single_cat_fn)
                single_all_cat = vstack([single_all_cat, single_cat])
            # flag out source w/ flag other than 0 and mag err > 0.198, which is S/N < 5
            single_all_cat = single_all_cat[(single_all_cat['FLAGS'] == 0) & (single_all_cat['MAGERR_BEST'] < 0.198)]
            single_all_cat[ab.mag_sys] = single_all_cat[ab.mag_sys] + ab.short_a[k][i] + ab.short_b[k][i]
            single_coords = SkyCoord(single_all_cat['ALPHA_J2000'], single_all_cat['DELTA_J2000'], unit='deg')

            stacked_cat_fns = glob.glob(f'{ab.work_dir}sex/run_stack/{ab.clusters[k]}_{ab.bands[i]}si_*_psf_nan.cat')
            stacked_all_cat = Table()
            for stacked_cat_fn in stacked_cat_fns:
                stacked_cat = ascii.read(stacked_cat_fn)
                stacked_all_cat = vstack([stacked_all_cat, stacked_cat])
            stacked_all_cat = stacked_all_cat[(stacked_all_cat['FLAGS'] == 0) & (stacked_all_cat['MAGERR_BEST'] < 0.198)]
            stacked_coords = SkyCoord(stacked_all_cat['ALPHA_J2000'], stacked_all_cat['DELTA_J2000'], unit='deg')

            idx, d2d, d3d = single_coords.match_to_catalog_sky(stacked_coords)
            sep_constraint = (d2d.arcsec < max_sep)
            single_matches = single_all_cat[sep_constraint]
            stacked_matches = stacked_all_cat[idx[sep_constraint]]

            num_process = 16
            num_gal_r = len(stacked_matches)
            ni = int(num_gal_r / num_process)
            start_time = time.time()
            with Pool(num_process) as p:
                is_dqm_r_zero = p.map(dqm_check, range(num_process))
            is_dqm_r_zero = np.array(is_dqm_r_zero).reshape(-1)

            single_matches = single_matches[0:len(is_dqm_r_zero)][is_dqm_r_zero]
            stacked_matches = stacked_matches[0:len(is_dqm_r_zero)][is_dqm_r_zero]

            single_matches.write(f'{ab.work_dir}catalogue/matched_flag0_5sig_dqm/{ab.clusters[k]}_{ab.bands[i]}_single_nopsf.cat',
                                 format='ascii', overwrite=True)
            stacked_matches.write(f'{ab.work_dir}catalogue/matched_flag0_5sig_dqm/{ab.clusters[k]}_{ab.bands[i]}_stack_nopsf.cat',
                                  format='ascii', overwrite=True)
        else:
            start_time = time.time()
            single_matches = ascii.read(f'{ab.work_dir}catalogue/matched_flag0_5sig_dqm/{ab.clusters[k]}_{ab.bands[i]}_single.cat')
            stacked_matches = ascii.read(f'{ab.work_dir}catalogue/matched_flag0_5sig_dqm/{ab.clusters[k]}_{ab.bands[i]}_stack.cat')

        is_inside = (single_matches['X_IMAGE'] > margin) & (single_matches['X_IMAGE'] < ccd_x - margin) & \
                    (single_matches['Y_IMAGE'] > margin) & (single_matches['Y_IMAGE'] < ccd_y - margin)

        single_matches = single_matches[is_inside]
        stacked_matches = stacked_matches[is_inside]

        is_complete = (single_matches['MAG_BEST'] < ab.short_m90[k][i]) & (stacked_matches['MAG_BEST'] < ab.stack_m90[k][i])

        single_matches = single_matches[is_complete]
        stacked_matches = stacked_matches[is_complete]

        is_galaxy = (single_matches['CLASS_STAR'] < 0.05) & (stacked_matches['CLASS_STAR'] < 0.05)

        single_matches = single_matches[is_galaxy]
        stacked_matches = stacked_matches[is_galaxy]

        c, low, upp = sigmaclip(single_matches[ab.mag_sys] - stacked_matches[ab.mag_sys], sig, sig)

        ind = np.where((single_matches[ab.mag_sys] - stacked_matches[ab.mag_sys] > low) &
                       (single_matches[ab.mag_sys] - stacked_matches[ab.mag_sys] < upp))
        ind2 = np.where((single_matches[ab.mag_sys] - stacked_matches[ab.mag_sys] < low) |
                        (single_matches[ab.mag_sys] - stacked_matches[ab.mag_sys] > upp))

        Data_Frame = {'delm': single_matches[ab.mag_sys][ind] - stacked_matches[ab.mag_sys][ind],
                      'am': np.zeros(len(ind[0])),
                      'err': single_matches[ab.magerr_sys][ind]**2 + stacked_matches[ab.magerr_sys][ind]**2
                      }

        df = pd.DataFrame(Data_Frame, columns=['delm', 'am', 'err'])

        df['err'][df['err'] == 0] = np.max(df['err'])

        X = df[['am']]
        Y = df['delm']

        # with statsmodels
        X = sm.add_constant(X)

        model = sm.WLS(Y, X, weights=1 / df['err']).fit()
        # results = model.fit()
        prstd, iv_l, iv_u = wls_prediction_std(model)

        print("{} {} {} {} {}".format(ab.clusters[k], ab.bands[i], model.params[0], model.bse[0], len(ind[0])))


        axs[k, i].scatter(single_matches[ab.mag_sys][ind],
                          single_matches[ab.mag_sys][ind] - stacked_matches[ab.mag_sys][ind],
                          alpha = 0.1,
                          s = 0.1)
        axs[k, i].scatter(single_matches[ab.mag_sys][ind2],
                          single_matches[ab.mag_sys][ind2] - stacked_matches[ab.mag_sys][ind2],
                          alpha = 0.1,
                          s = 0.1,
                          color = 'gray')
        if k == 0:
            axs11[i].scatter(single_matches[ab.mag_sys][ind],
                              single_matches[ab.mag_sys][ind] - stacked_matches[ab.mag_sys][ind],
                              alpha = 0.1,
                              s = 0.1)
            axs11[i].scatter(single_matches[ab.mag_sys][ind2],
                              single_matches[ab.mag_sys][ind2] - stacked_matches[ab.mag_sys][ind2],
                              alpha = 0.1,
                              s = 0.1,
                              color = 'gray')
            axs11[i].plot([12, 30], [low, low], '--', alpha=0.5, color='salmon')
            axs11[i].plot([12, 30], [upp, upp], '--', alpha=0.5, color='salmon')
            axs11[i].plot([12, 30], [model.params[0], model.params[0]], '--', alpha=0.7)
            axs11[i].fill_between([12, 30],
                                   [model.params[0] + model.bse[0], model.params[0] + model.bse[0]],
                                   [model.params[0] - model.bse[0], model.params[0] - model.bse[0]],
                                   alpha=0.3,
                                   color='blue'
                                   )
            axs11[i].set_xlim([12, 30])
            # axs11[i].set_ylabel(r'$\Delta m$ ($m_{single} - m_{combined}$)', fontsize=fs)
            if i == 0:
                axs11[i].set_ylabel(r'$a\prime$', fontsize=fs)
            axs11[i].set_xlabel(rf'$m_{{ {ab.bands[i]} \prime}}$', fontsize=fs)
            # axs11[i].text(0.1, 0.8, rf'${ab.bands[i]}\prime$', transform=axs[k, i].transAxes, size=30)
            axs11[i].set_ylim([model.params[0] - 1, model.params[0] + 1])
            axs11[i].tick_params(direction='in', top=True, right=True)
            if i == 2:
                axs11[i].yaxis.set_label_position('right')
                axs11[i].set_ylabel(f'{ab.clusters[k]}', fontsize=fs, rotation=270, labelpad=20)
            box = axs11[i].get_position()
            box.y0 = box.y0 + 0.08
            box.y1 = box.y1 + 0.08
            axs11[i].set_position(box)



        axs2[k, i].scatter(single_matches[ab.magerr_sys],
                           stacked_matches[ab.magerr_sys],
                           alpha = 0.1,
                           s = 0.1)
        axs3[k, i].scatter(single_matches[ab.mag_sys],
                           stacked_matches[ab.mag_sys] + model.params[0],
                           alpha = 0.1,
                           s = 0.1)
        axs[k, i].plot([12, 30], [low, low], '--', alpha=0.5, color='salmon')
        axs[k, i].plot([12, 30], [upp, upp], '--', alpha=0.5, color='salmon')
        axs[k, i].plot([12, 30], [model.params[0], model.params[0]], '--', alpha=0.7)
        axs[k, i].fill_between([12, 30],
                         [model.params[0] + model.bse[0], model.params[0] + model.bse[0]],
                         [model.params[0] - model.bse[0], model.params[0] - model.bse[0]],
                         alpha=0.3,
                         color='blue'
                         )

        axs2[k, i].plot([0, 0.5], [0, 0.5], '--', alpha=0.1)
        axs3[k, i].plot([12, 30], [12, 30], '--', alpha=0.1)
        axs[k, i].set_xlim([12, 30])
        axs3[k, i].set_xlim([12, 30])
        # axs[k, i].set_ylim([12, 30])
        axs2[k, i].set_xlim([0, 0.5])
        axs2[k, i].set_ylim([0, 0.5])
        #axs[k, i].set_title(ab.clusters[k])
        axs1[k, i].set_title(ab.clusters[k])
        axs2[k, i].set_title(ab.clusters[k])
        if i == 0 and k == 3:
            # axs[k, i].set_ylabel(r'$\Delta m$ ($m_{single} - m_{combined}$)', fontsize=fs)
            axs[k, i].set_ylabel(r'$a\prime$', fontsize=fs)
        if i == 1 and k == 6:
            axs[k, i].set_xlabel(r'$m_{single}$', fontsize=fs)
        if k < 6:
            plt.setp(axs[k, i].get_xticklabels(), visible=False)
        if i == 2:
            axs[k, i].yaxis.set_label_position('right')
            axs[k, i].set_ylabel(f'{ab.clusters[k]}', fontsize=fs, rotation=270, labelpad=20)
        if k == 0:
            axs[k, i].text(0.1, 0.8, rf'${ab.bands[i]}\prime$', transform = axs[k, i].transAxes, size=30)
        axs2[k, i].set_ylabel('stacked')
        axs2[k, i].set_xlabel('single exposure')
        # axs[k, i].invert_xaxis()
        # axs[k, i].invert_yaxis()
        # axs[k, i].text(0.1, 0.85,
        #                ab.bands[i]+' = '+ab.bands[i]+'_inst (ZP=25) + {:6.4f}'.format(model.params[0]),
        #                transform = axs[k, i].transAxes)
        # # axs[k, i].text(29, 14, r'$\chi^2$/# = {:7.3f}'.format(minchi/len(single_matches)))
        # axs[k, i].text(0.1, 0.7,
        #                '# of sample = {:4d}'.format(len(ind[0])),
        #                transform = axs[k, i].transAxes)
        # axs[k, i].legend(title='Standard Star Field (Airmass)(#)', loc='lower right', prop={'size':7})
        axs[k, i].set_ylim([model.params[0]-1, model.params[0]+1])
        axs[k, i].tick_params(direction='in', top=True, right=True)
        axs2[k, i].text(0.1, 0.4, ab.bands[i])

        axs1[k, i].hist(single_matches[ab.mag_sys],
                       bins=np.arange(12, 30, 0.5),
                       label='single',
                       histtype='step')
        axs1[k, i].hist(stacked_matches[ab.mag_sys] + model.params[0],
                       bins=np.arange(12.1, 30.1, 0.5),
                       label='stacked',
                       histtype='step')
        axs1[k, i].legend(loc='upper left', fontsize='small')

        print("%s %s %5.2f" % (ab.clusters[k], ab.bands[i], model.params[0]))
        # print(a_guess)
        print("--- %s minutes ---" % ((time.time() - start_time)/60.0))



fig11.savefig(plot_dir+'DECam_stacked_vs_single_best_dqm_5sig_flag_inside_complete_005galaxy_A754.png')
fig.savefig(plot_dir+'DECam_stacked_vs_single_best_dqm_5sig_flag_inside_complete_005galaxy.png')
fig1.savefig(plot_dir+'DECam_stacked_vs_single_best_hist_dqm_5sig_flag_inside_complete_005galaxy.png')
fig2.savefig(plot_dir+'DECam_stacked_vs_single_best_err_dqm_5sig_flag_inside_complete_005galaxy.png')
fig3.savefig(plot_dir+'DECam_stacked_vs_single_best_1on1_dqm_5sig_flag_inside_complete_005galaxy.png')
# plt.show(block=False)

