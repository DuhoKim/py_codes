from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astropy.coordinates import SkyCoord
import astropy.units as u
import time
import my_module
from astropy.table import Table, vstack
import statsmodels.api as sm
from statsmodels.sandbox.regression.predstd import wls_prediction_std
import pandas as pd
from scipy.stats import sigmaclip
import math

work_dir=("/Users/duhokim/work/abell/")
ss_dir=("/Users/duhokim/old_macbook_documents/work/cat/SouthernStandardStars/remove_asterisk/")
check_dir=("/Volumes/APPLE SSD/Server/work/sex/ss_deblend/")

dates = ['10apr', '11apr', '19aug', '20aug', '21aug']
title = ['2013-04-10', '2013-04-11', '2014-08-19', '2014-08-20', '2014-08-21']
bands = ['u', 'g', 'r']
exp_t = [100, 10, 10]
mag_col = ['col4', 'col7', 'col10']
err_col = ['col5', 'col8', 'col11']

# column name for u, g, r in SDSS Std catalogue
std_mag_col = ['psfMag_u', 'psfMag_g', 'psfMag_r']     # the median magnitude
std_err_col = ['psfMagErr_u', 'psfMagErr_g', 'psfMagErr_r']      # the standard error for mean magnitudes which is x 1/1.25 of median

sig = 3.0
fs = 20     # fontsize
fst = 15    # fontsize for title

max_sep = 1.0
class_star_lim = 0.9
flag_lim = 0
near_nei_dis_lim = 1.0  # in aperture radii which is 7.43"

fig = plt.figure(figsize=(10, 12))

gs0 = gridspec.GridSpec(2, 1)

gs00 = gridspec.GridSpecFromSubplotSpec(3, 5, subplot_spec=gs0[0], hspace=0, wspace=0)

gs01 = gs0[1].subgridspec(3, 5, hspace=0, wspace=0)

for n in range(0, len(bands)):
    for m in range(0, len(dates)):
        ax1 = fig.add_subplot(gs00[n, m])
        ax2 = fig.add_subplot(gs01[n, m])

        if n == 0:
            ax1.set_ylim([-2.7, -1.5])
            ax1.set_yticks(np.arange(-2.5, -1.5, 0.5))
            ax2.set_ylim([-2, -1])
            ax2.set_yticks(np.arange(-1.8, -1, 0.4))
        elif n == 1:
            ax1.set_ylim([-0.2, 0.5])
            ax1.set_yticks(np.arange(-0.1, 0.4, 0.2))
            ax2.set_ylim([0, 1])
            ax2.set_yticks(np.arange(0.1, 0.9, 0.3))
        else:
            ax1.set_ylim([0.1, 0.8])
            ax2.set_ylim([0, 1])

        ax1.set_xlim([1, 2])
        ax2.set_xlim([11, 22])
        ax1.set_xticks(np.arange(0, 1.8, 0.6))



        if n == 0:
            if m == 0:
                ax1.axis('off')
                ax2.axis('off')
                # ax1.set_ylabel(r'$\Delta m$', color='white')
                # ax2.set_ylabel(r'$\Delta m$ (Airmass term corrected)', color='white')
            ax1.set_title(title[m], fontsize=fst)

        if m == 0:
            if n == 0:
                continue
            ss_name = ['E4-A_0040', 'Std4_0208', 'E4-A_0211', 'Std4_0337', 'Std5_0340', 'Std6_0507',
                       'E6-A_0510', 'Std5_0636', 'E6-A_0639', 'E6-A_0732', 'E5-A_0735', 'Std6_0925']
            ss_cat_name = ['E4_a', 'Std4', 'E4_a', 'Std4', 'Std5', 'Std6',
                           'E6_a', 'Std5', 'E6_a', 'E6_a', 'E5_a', 'Std6']
            ss_airmass = [1.04, 1.16, 1.08, 1.24, 1.16, 1.2, 1.06, 1.53, 1.04, 1.08, 1.54, 1.78]
            ss_seeing = [0.8] * len(ss_name)
        elif m == 1:
            if n == 0:
                ss_name = ['E4-A_0425', 'Std5_0430', 'E4-A_0602', 'Std6_0714',
                           'Std6_0814', 'Std6_0841', 'E4-A_0907', 'Std6_0912']
                ss_cat_name = ['E4_a', 'Std5', 'E4_a', 'Std6', 'Std6', 'Std6', 'E4_a', 'Std6']
                ss_airmass = [1.1, 1.18, 1.04, 1.21, 1.37, 1.5, 1.27, 1.7]
            else:
                ss_name = ['E4-A_0046', 'Std4_0214', 'E4-A_0217', 'Std4_0328', 'Std5_0332', 'E4-A_0425', 'Std5_0430',
                           'E4-A_0602', 'Std6_0714', 'Std6_0814', 'Std6_0841', 'E4-A_0907', 'Std6_0912']
                ss_cat_name = ['E4_a', 'Std4', 'E4_a', 'Std4', 'Std5', 'E4_a', 'Std5',
                               'E4_a', 'Std6', 'Std6', 'Std6', 'E4_a', 'Std6']
                ss_airmass = [1.04, 1.16, 1.09, 1.23, 1.16, 1.1, 1.18, 1.04, 1.21, 1.37, 1.5, 1.27, 1.7]
            ss_seeing = [0.8] * len(ss_name)
        elif m == 2:
            ss_name = ['220100', 'E8-A', 'LSE_259', '190000']
            ss_legend_name = ['220100 (1.96)', 'E8-A (1.21)', 'LSE_259 (1.75)', '190000 (1.23)']
            ss_cat_name = ['220100-300000', 'E8_a', 'LSE_259', '190000-295600']
            ss_airmass = [1.96, 1.21, 1.75, 1.23]
            ss_seeing = [1.04, 0.95, 1.03, 0.92]

        elif m == 3:
            ss_name = ['E5-A', 'E6-A', 'LSE_259_1_11', '220100_2_44', '220100_1_01', 'LSE_259_1_71']
            ss_legend_name = ['E5-A (1.7)', 'E6-A (1.12)', 'LSE_259 (1.11)', '220100 (2.44)', '220100 (1.01)', 'LSE_259 (1.71)']
            ss_cat_name = ['E5_a', 'E6_a', 'LSE_259', '220100-300000', '220100-300000', 'LSE_259']
            ss_airmass = [1.7, 1.12, 1.11, 2.44, 1.01, 1.71]
            ss_seeing = [1.75, 1.49, 1.31, 1.77, 3.5, 3.8]

        elif m == 4:
            ss_name = ['E5-A', 'E6-A', 'LSE_259_1_11', 'LSE_259_1_8']
            ss_legend_name = ['E5-A (1.65)', 'E6-A (1.11)', 'LSE_259 (1.11)', 'LSE_259 (1.8)']
            ss_cat_name = ['E5_a', 'E6_a', 'LSE_259', 'LSE_259']
            ss_airmass = [1.65, 1.11, 1.11, 1.8]
            ss_seeing = [1.46, 1.03, 1.07, 1.58]

        delm = []
        am = []
        merr = []
        mag = []
        for i in range(0, len(ss_name)):
            # sex_result = ascii.read(work_dir+'sex/cat/ss_deblend/'+dates[m]+'/'+ss_name[i] + '_' + bands[n] + 'a.cat')
            if m < 2:
                sex_result = ascii.read('/Volumes/APPLE SSD/Server/work/sex/ss2/' + dates[m] + '/' + ss_name[i] + '_' +
                                        bands[n] + 'a.cat')
            else:
                sex_result = ascii.read('/Users/duhokim/old_macbook_documents/work/sex/'+dates[m]+'/'+ss_name[i] + '_' +
                                        bands[n] + 'a.cat')
            # sex_result_good = (sex_result['MAG_ISO'] < 90) & (sex_result['FLAGS'] <= flag_lim) & \
            #                   (sex_result['CLASS_STAR'] > class_star_lim) & ~np.isnan(sex_result['FLUX_APER_5']) & \
            #                   ~np.isnan(sex_result['FLUX_APER_6']) & ~np.isnan(sex_result['FLUXERR_APER_5']) & \
            #                   ~np.isnan(sex_result['FLUXERR_APER_6'])
            sex_result_good = (sex_result['MAG_ISO'] < 90) & (sex_result['FLAGS'] <= flag_lim) & \
                              (sex_result['CLASS_STAR'] > class_star_lim) & (sex_result['FLUX_APER_5'] > 0) & \
                              (sex_result['FLUX_APER_6'] > 0) & (sex_result['FLUXERR_APER_5'] > 0) & \
                              (sex_result['FLUXERR_APER_6'] > 0)
            sex_coords = SkyCoord(sex_result['ALPHA_J2000'][sex_result_good],
                                  sex_result['DELTA_J2000'][sex_result_good],
                                  unit='deg')

            if ss_cat_name[i][0:3] == 'Std':
                ss_cat = ascii.read(work_dir + 'cat/' + ss_cat_name[i] + '.csv')
                cat_coords = SkyCoord(ss_cat['ra'], ss_cat['dec'], unit='deg')
            else:
                ss_cat = ascii.read(ss_dir + ss_cat_name[i] + '.dat.trimmed', format='no_header')
                cat_coords = SkyCoord(ss_cat['col2'], ss_cat['col3'], unit=(u.hourangle, u.deg))

            idx, d2d, d3d = sex_coords.match_to_catalog_sky(cat_coords)

            if ss_cat_name[i][0:3] == 'Std':
                sep_constraint = (d2d.arcsec < max_sep)
            else:
                sep_constraint = (d2d.arcsec < max_sep) & (ss_cat[idx][err_col[n]] > 0) & (
                        ss_cat[idx]['col20'] > near_nei_dis_lim)
            sex_matches = sex_result[sex_result_good][sep_constraint]
            cat_matches = ss_cat[idx[sep_constraint]]

            if ss_cat_name[i][0:3] == 'Std':
                mag_col_txt = std_mag_col[n]
                err_col_txt = std_err_col[n]
            else:
                mag_col_txt = mag_col[n]
                err_col_txt = err_col[n]

            mag5 = 25 - 2.5 * np.log10(sex_matches['FLUX_APER_5'] / exp_t[n])
            mag6 = 25 - 2.5 * np.log10(sex_matches['FLUX_APER_6'] / exp_t[n])

            # temporarily save mag after aperture correction
            sex_matches['MAG_ISO'] = mag5 * 0.186 + mag6 * 0.814 - 0.05 * ss_seeing[i] + 0.04
            not_nan_ind = np.where( (sex_matches['FLUX_APER_6'] - sex_matches['FLUXERR_APER_6']) > 0 )

            sex_matches['MAGERR_ISO'][not_nan_ind] = - 2.5 * np.log10(
                (sex_matches['FLUX_APER_6'][not_nan_ind] - sex_matches['FLUXERR_APER_6'][not_nan_ind]) / exp_t[n]) \
                                        + 2.5 * np.log10(
                (sex_matches['FLUX_APER_6'][not_nan_ind] + sex_matches['FLUXERR_APER_6'][not_nan_ind]) / exp_t[n])

            delm.extend((cat_matches[mag_col_txt] - sex_matches['MAG_ISO']).data.tolist())
            merr.extend((cat_matches[err_col_txt] ** 2 + sex_matches['MAGERR_ISO'] ** 2).data.tolist())
            mag.extend(cat_matches[mag_col_txt].data.tolist())

            delm_print = (cat_matches[mag_col_txt] - sex_matches['MAG_ISO']).data.tolist()

            am.extend([ss_airmass[i]] * len(sex_matches))

            print('{} {} {} {}'.format(ss_name[i], ss_airmass[i], len(sex_matches), delm_print))

        c, low, upp = sigmaclip(delm, sig, sig)

        ind = np.where((delm > low) & (delm < upp))
        ind2 = np.where((delm < low) | (delm > upp))

        delm = np.array(delm)
        am = np.array(am)
        merr = np.array(merr)
        mag = np.array(mag)

        Data_Frame = {'delm': delm[ind],
                      'am': am[ind],
                      'err': merr[ind]
                      }

        df = pd.DataFrame(Data_Frame, columns=['delm', 'am', 'err'])

        X = df[['am']]
        Y = df['delm']

        # with statsmodels
        X = sm.add_constant(X)

        model = sm.WLS(Y, X, weights=1/df['err']).fit()
        # results = model.fit()
        prstd, iv_l, iv_u = wls_prediction_std(model)

        print("{} {} {} {} {} {} {}".format(dates[m], bands[n], model.params[0], model.bse[0],
                                            model.params[1], model.bse[1], len(ind[0])))

        ax1_xran = [min(am), max(am)]
        ax2_xran = [min(mag), max(mag)]

        ax1.plot(ax1_xran, [f * model.params[1] + model.params[0] for f in ax1_xran], '--', alpha=0.2)
        ax2.plot(ax2_xran, [model.params[0], model.params[0]], '--', alpha=0.2)
        ax1.plot(ax1_xran, [low, low], '--', color='salmon', alpha=0.2)
        ax1.plot(ax1_xran, [upp, upp], '--', color='salmon', alpha=0.2)
        ax1.fill_between(ax1_xran,
                            [f * (model.params[1]+model.bse[1]) + model.params[0] for f in ax1_xran],
                            [f * (model.params[1]-model.bse[1]) + model.params[0] for f in ax1_xran],
                            alpha=0.1,
                            color='blue'
                            )
        ax2.fill_between(ax2_xran,
                            [model.params[0]+model.bse[0]]*2,
                            [model.params[0]-model.bse[0]]*2,
                            alpha=0.1,
                            color='blue'
                            )
        # ax2.set_xlim([12, 20])

        ax1.scatter(am[ind], delm[ind], alpha=0.3, s=5)
        ax2.scatter(mag[ind], delm[ind] - model.params[1] * am[ind], alpha=0.3, s=5)
        ax1.scatter(am[ind2], delm[ind2], alpha=0.3, s=5, color='gray')
        ax2.scatter(mag[ind2], delm[ind2] - model.params[1] * am[ind2], alpha=0.3, s=5, color='gray')

        # if m == 0 and n == 2:
        #     a = 1

        if m == 0 and n==1:
            ax1.set_ylabel(r'$\Delta m$', fontsize=fs)
            # ax2.set_ylabel(r'$\Delta m$ (Airmass term corrected)', fontsize=fs)
            ax2.set_ylabel('a', fontsize=fs)
        # elif m == 1 and n == 0:
        #     ax1.set_ylabel(r'$\Delta m$')
        #     ax2.set_ylabel(r'$\Delta m$ (Airmass term corrected)')
        elif m == 2 and n==2:
            ax2.set_xlabel(r'$m$', fontsize=fs)
            ax1.set_xlabel('Airmass', fontsize=fs)
        elif m == 4:
            ax1.text(0.7, 0.7, r'$' + bands[n] + '\prime$', transform=ax1.transAxes, size=fs)
            ax2.text(0.7, 0.7, r'$' + bands[n] + '\prime$', transform=ax2.transAxes, size=fs)

        ax1.tick_params(direction='in', top=True, right=True)
        ax2.tick_params(direction='in', top=True, right=True)

        if n < 2:
            plt.setp(ax1.get_xticklabels(), visible=False)
            plt.setp(ax2.get_xticklabels(), visible=False)
        if (m > 0):
            if (m != 1) | (n != 0):
                plt.setp(ax1.get_yticklabels(), visible=False)
                plt.setp(ax2.get_yticklabels(), visible=False)


    # if n == 0:
    #     gs.tight_layout(fig, rect=[0.047, 1 - (1 + n) / 3, 1, 1 - n / 3])
    # else:
    #     gs.tight_layout(fig, rect=[0, 1-(1+n)/3, 1, 1-n/3])

fig.savefig(work_dir+'plot/'+'DECam_standardized_vs_trimmed_aper_corr_flag0_{}sigclip_default_2013_new_fig.png'.format(sig))

