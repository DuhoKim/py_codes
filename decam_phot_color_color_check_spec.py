from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import abell_cluster_module as ab
import importlib
importlib.reload(ab)
import statsmodels.api as sm
from astropy.cosmology import Planck18 as Cosmo
from astropy.coordinates import SkyCoord
from scipy import constants as const
import calc_kcor as kcor
from matplotlib import gridspec

best_a = -0.039476
best_b = 1.3851
rs_wid = 0.1

size = 15
alpha = 0.5

size2 = 30
alpha2 = 0.15

g_yran = [0, 1.5]
u_yran = [1, 4]

cols = ['black', 'red', 'purple', 'darkorange', 'olive', 'darkgreen', 'teal']

params = ['MEAN_MAG_DIFF', 'DATE', 'EXP_TIME', 'SEEING', 'AIRMASS', 'TOTAL_MATCH_NUM']
num_of_param = len(params)

am_cat = pd.read_csv('/Users/duhokim/work/abell/sex/cat/airmass.csv')

fig, axs = plt.subplots(2, len(ab.clusters), tight_layout=True, figsize=(20, 6))
fig2, axs2 = plt.subplots(2, 1, tight_layout=True, figsize=(8, 12))
# fig3, axs3 = plt.subplots(2, 4, tight_layout=True, figsize=(12, 6))
fig3 = plt.figure(tight_layout=True, figsize=(14, 6))
fig4, axs4 = plt.subplots(1, len(ab.clusters), tight_layout=True, figsize=(20, 3))
fig5, axs5 = plt.subplots(4, 2, tight_layout=True, figsize=(10, 20))
# fig6, axs6 = plt.subplots(2, 4, tight_layout=True, figsize=(12, 6))
fig6 = plt.figure(tight_layout=True, figsize=(12, 6))

gs3 = fig3.add_gridspec(2, 4, wspace=0.2, hspace=0)
gs = fig6.add_gridspec(2, 4, wspace=0, hspace=0)

for k in range(0, len(ab.clusters)):
# for k in range(0, 1):
    with open(ab.work_dir+f'spec/{ab.clusters[k]}_rs.radec', 'w') as radec:
        xran = [11, 21] - ab.distmod[k]

        sex_cat = ascii.read(ab.sex_dir+f'DECam_merged_SEx_cat_{ab.clusters[k]}_Gal_ext_corrected_{ab.ver}.txt')

        # good_g = (sex_cat['CLASS_STAR'] < 0.2) & (sex_cat[ab.mag_sys+'_g'] < 30) & (sex_cat[ab.mag_sys+'_r'] < xran[1]) & (sex_cat[ab.magerr_sys+'_r'] > 0)
        # good_u = (sex_cat['CLASS_STAR'] < 0.2) & (sex_cat[ab.mag_sys+'_u'] < 30) & (sex_cat[ab.mag_sys+'_r'] < xran[1]) & (sex_cat[ab.magerr_sys+'_r'] > 0)

        sex_cat[ab.mag_sys+'_u'] -= kcor.calc_kcor('u', ab.redshifts[k], 'u - r',
                                                    sex_cat[ab.mag_sys+'_u'] - sex_cat[ab.mag_sys+'_r'])
        sex_cat[ab.mag_sys + '_g'] -= kcor.calc_kcor('g', ab.redshifts[k], 'g - r',
                                                     sex_cat[ab.mag_sys + '_g'] - sex_cat[ab.mag_sys + '_r'])
        sex_cat[ab.mag_sys + '_r'] -= kcor.calc_kcor('r', ab.redshifts[k], 'g - r',
                                                     sex_cat[ab.mag_sys + '_g'] - sex_cat[ab.mag_sys + '_r'])

        ### cross match to spec
        spec = ascii.read(ab.work_dir+f'spec/{ab.clusters[k]}_spec.txt')

        coords_cat = SkyCoord(sex_cat['ALPHA_J2000'], sex_cat['DELTA_J2000'], unit='deg')
        coords_spec = SkyCoord(spec['col2'], spec['col3'], unit='deg')
        idx_spec, d2d, d3d = coords_cat.match_to_catalog_sky(coords_spec)
        matched_cat = (d2d.arcsec < ab.max_sep) & (sex_cat[ab.mag_sys+'_u'] < 900) & \
                      (sex_cat[ab.mag_sys+'_g'] < 900) & (sex_cat[ab.mag_sys+'_r'] < 900)

        u_spec = sex_cat[ab.mag_sys+'_u'][matched_cat]
        g_spec = sex_cat[ab.mag_sys+'_g'][matched_cat]
        r_spec = sex_cat[ab.mag_sys+'_r'][matched_cat]
        r_spec_abs = r_spec - ab.distmod[k]

        axs[0, k].scatter(r_spec_abs, g_spec - r_spec, alpha=alpha, s=size)
        axs[1, k].scatter(r_spec_abs, u_spec - r_spec, alpha=alpha, s=size)

        ### Projected Phase Space diagram
        kpc_arcmin = Cosmo.kpc_proper_per_arcmin(ab.redshifts[k])
        sep = coords_cat[matched_cat].separation(ab.coords_cl[k])   # in degrees
        sig_z = np.std(spec['col4'])
        xx = (sep * 6e1 * kpc_arcmin / 1e3 / ab.r200[k]).value
        yy = np.abs(spec['col4'][idx_spec[matched_cat]] - ab.redshifts[k])/sig_z/(1+ab.redshifts[k])
        in_vir = yy < (-1.5 / 1.2 * xx + 1.5)
        out_vir = yy > (-1.5 / 1.2 * xx + 1.5)

        axs4[k].scatter(xx, yy, alpha=alpha, s=size)
        axs4[k].plot([0, 1.2], [1.5, 0], linestyle='--')

        axs4[k].set_title(ab.clusters[k])
        axs4[k].set_xlim([0, 3])
        axs4[k].set_ylim([0, 3])
        axs4[k].set_xlabel('r/R200')
        axs4[k].set_ylabel(r'$\Delta$z/(1+z)$\sigma_z$')
        axs4[k].text(1.7, 2.5, f'N_tot={len(sep)}')
        axs4[k].text(1.7, 2.2, f'N_in={sum(in_vir)}')
        axs4[k].text(1.7, 1.9, f'N_out={sum(out_vir)}')

        ### radial vel histogram
        row = int(k / 2)
        col = k % 2
        axs5[row, col].hist(spec['col4'] * const.c / 1e3,
                            histtype='step')
        axs5[row, col].set_title(ab.clusters[k])
        axs5[row, col].annotate(r'$\sigma_{rv}$='+f"{np.std(spec['col4'] * const.c / 1e3):.2f}",
                                xy=(0.1, 0.9),
                                xycoords='axes fraction',
                                fontsize=16)


        ### histogram
        row = int(k / 4)
        col = k % 4
        axs3 = fig3.add_subplot(gs3[row, col])
        r_abs_kcor = sex_cat[ab.mag_sys+'_r'] - ab.distmod[k]
        num_sam = np.sum((r_abs_kcor > -25) & (r_abs_kcor < -20))
        num_spec_sam = np.sum((r_abs_kcor[matched_cat] > -25) & (r_abs_kcor[matched_cat] < -20))


        in_rs = (sex_cat[ab.mag_sys + '_g'] - sex_cat[ab.mag_sys + '_r'] < best_a * sex_cat[ab.mag_sys + '_r'] + best_b + rs_wid) & \
                  (sex_cat[ab.mag_sys + '_g'] - sex_cat[ab.mag_sys + '_r'] > best_a * sex_cat[ab.mag_sys + '_r'] + best_b - rs_wid) & \
                (r_abs_kcor < -20)

        num_sam_rs = np.sum(in_rs)

        for i in range(0, num_sam_rs):
            radec.writelines(f"{sex_cat['ALPHA_J2000'][in_rs][i]} {sex_cat['DELTA_J2000'][in_rs][i]} \n")

        axs3.hist(r_abs_kcor, range=(-25, -20), histtype='step', color='blue')
        axs3.hist(r_abs_kcor[in_rs], range=(-25, -20), histtype='step', color='red')
        axs3.hist(r_abs_kcor[matched_cat], range=(-25, -20), histtype='step', color='green')
        # axs3.set_ylabel('N')#, color = 'blue')
        axs3.annotate(r'N$_{phot}$='+f'{num_sam}', xy=(0.05, 0.75), xycoords='axes fraction', fontsize=10, color='blue')
        axs3.annotate(r'N$_{rs}$='+f'{num_sam_rs}', xy=(0.05, 0.65), xycoords='axes fraction', fontsize=10, color='red')
        axs3.annotate(r'N$_{spec}$=' + f'{num_spec_sam}', xy=(0.05, 0.55), xycoords='axes fraction', fontsize=10, color='green')
        # axs3.set_xlabel(r'M$_r$')
        axs3.set_xlim([-25, -20])
        # axs3.set_title(ab.clusters[k])
        axs3.annotate(f'{ab.clusters[k]}', xy=(0.05, 0.85), xycoords='axes fraction', fontsize=20)

        if row == 0 and col < 3:
            plt.setp(axs3.get_xticklabels(), visible=False, fontsize=20)
        else:
            axs3.set_xlabel("r'", fontsize=20)
        if col == 0:
            # plt.setp(axs3.get_yticklabels(), visible=False, fontsize=20)
        # else:
            axs3.set_ylabel("N", fontsize=20)

        ### CMD
        merr_good_g = np.array(np.sqrt(sex_cat[ab.magerr_sys + '_g'][matched_cat] ** 2 +
                                       sex_cat[ab.magerr_sys + '_r'][matched_cat] ** 2 ))
        merr_good_u = np.array(np.sqrt( sex_cat[ab.magerr_sys + '_u'][matched_cat] ** 2 +
                                        sex_cat[ab.magerr_sys + '_r'][matched_cat] ** 2))
        merr_good_g[np.where(merr_good_g < 0.1)] = 0.1
        merr_good_u[np.where(merr_good_u < 0.1)] = 0.1
        # merr_good_g = np.zeros(len(r_spec))
        # merr_good_u = np.zeros(len(r_spec))
        # merr_good_g[:] = 1
        # merr_good_u[:] = 1

        Data_Frame = {'g_r': np.array(g_spec - r_spec), 'u_r': np.array(u_spec - r_spec), 'r_abs': np.array(r_spec_abs),
                      'r_app': np.array(r_spec), 'g_r_err': merr_good_g, 'u_r_err': merr_good_u}

        df = pd.DataFrame(Data_Frame, columns=['g_r', 'u_r', 'r_abs', 'r_app', 'g_r_err', 'u_r_err'])

        # with statsmodels
        X = sm.add_constant(df[['r_abs']])
        X2 = sm.add_constant(df[['r_app']])

        model = sm.WLS(df['g_r'], X, weights=1/df['g_r_err']).fit()
        model2 = sm.WLS(df['u_r'], X, weights=1 / df['u_r_err']).fit()
        model6 = sm.WLS(df['g_r'], X2, weights=1 / df['g_r_err']).fit()

        print(model.summary())
        # print(model2.summary())

        print("{} g-r {} {} {} {}".format(ab.clusters[k], model6.params[0], model6.bse[0], model6.params[1], model6.bse[1]))
        # print("{} u-r {} {} {} {}".format(ab.clusters[k], model2.params[0], model2.bse[0], model2.params[1], model2.bse[1]))

        axs2[0].scatter(r_spec_abs, g_spec - r_spec, label=f'{ab.clusters[k]}', alpha=alpha2, s=size2, color=cols[k])
        axs2[0].plot(xran, [f * model.params[1] + model.params[0] for f in xran], '--', alpha=0.8, color=cols[k], label=f'{ab.clusters[k]}')
        axs2[0].fill_between(xran,
                            [f * (model.params[1]) + model.params[0] - model.bse[0] for f in xran],
                            [f * (model.params[1]) + model.params[0] + model.bse[0] for f in xran],
                            alpha=0.1, color=cols[k])
        axs2[1].scatter(r_spec_abs, u_spec - r_spec, alpha=alpha2, s=size2, color=cols[k])
        axs2[1].plot(xran, [f * model2.params[1] + model2.params[0] for f in xran], '--', alpha=0.8, color=cols[k], label=f'{ab.clusters[k]}')
        axs2[1].fill_between(xran,
                             [f * (model2.params[1]) + model2.params[0] - model2.bse[0] for f in xran],
                             [f * (model2.params[1]) + model2.params[0] + model2.bse[0] for f in xran],
                             alpha=0.1, color=cols[k])

        xran6 = [13, 20]
        size2 = 10
        alpha2 = 0.1

        # axs6 = plt.subplot(gs[row, col])
        axs6 = fig6.add_subplot(gs[row, col])

        axs6.scatter(r_spec, g_spec - r_spec, label=f'{ab.clusters[k]}', alpha=0.4, s=5, color='teal',
                               facecolors='none')
        axs6.scatter(sex_cat[ab.mag_sys + '_r'], sex_cat[ab.mag_sys + '_g'] - sex_cat[ab.mag_sys + '_r'],
                               label=f'{ab.clusters[k]}', alpha=0.2, s=1, color='grey', facecolors='none')
        # axs6[row, col].plot(xran6, [f * model6.params[1] + model6.params[0] for f in xran6], '--', alpha=0.8, color='grey',
        #              label=f'{ab.clusters[k]}')
        # axs6[row, col].fill_between(xran6,
        #                      [f * (model6.params[1]) + model6.params[0] - model6.bse[0] for f in xran6],
        #                      [f * (model6.params[1]) + model6.params[0] + model6.bse[0] for f in xran6],
        #                      alpha=0.1, color='grey')

        axs6.plot(xran6, [f * best_a + best_b for f in xran6], '--', alpha=0.8, color='grey',
                            label=f'{ab.clusters[k]}')
        # axs6[row, col].fill_between(xran6,
        #                             [f * best_a + best_b - 0.1 for f in xran6],
        #                             [f * best_a + best_b + 0.1 for f in xran6],
        #                             alpha=0.1, color='grey')
        axs6.plot(xran6, [f * best_a + best_b - 0.1 for f in xran6], '--', alpha=0.2, color='grey',
                            label=f'{ab.clusters[k]}')
        axs6.plot(xran6, [f * best_a + best_b + 0.1 for f in xran6], '--', alpha=0.2, color='grey',
                            label=f'{ab.clusters[k]}')
        axs6.axvline(x=-20 + ab.distmod[k], linestyle='--', color='grey')

        if row == 0 and col < 3:
            plt.setp(axs6.get_xticklabels(), visible=False, fontsize=20)
        else:
            axs6.set_xlabel("r'", fontsize=20)
        if col:
            plt.setp(axs6.get_yticklabels(), visible=False, fontsize=20)
        else:
            axs6.set_ylabel("g' - r'", fontsize=20)

        axs6.set_xlim([13, 20])
        axs6.set_ylim([0, 1.1])
        axs6.annotate(f'{ab.clusters[k]}', xy=(0.05, 0.05), xycoords='axes fraction', fontsize=20)

        plt.subplots_adjust(hspace=.0, wspace=.0)



        # for i in range(0, 5):
        #     color = cols[k] if i == 4 else 'grey'
        #     in_2sig = (g_spec - r_spec < model6.params[1] * r_spec + model6.params[0] + 2 * model6.bse[0]) & \
        #               (g_spec - r_spec > model6.params[1] * r_spec + model6.params[0] - 2 * model6.bse[0])
        #
        #     New_Data_Frame = {'g_r': np.array(g_spec[in_2sig] - r_spec[in_2sig]), 'r_app': np.array(r_spec[in_2sig]),
        #                       'g_r_err': merr_good_g[in_2sig]}
        #     df = pd.DataFrame(New_Data_Frame, columns=['g_r', 'r_app', 'g_r_err'])
        #     X2 = sm.add_constant(df[['r_app']])
        #     model6 = sm.WLS(df['g_r'], X2, weights=1 / df['g_r_err']).fit()
        #     print("{} g-r {} {} {} {} {}".format(ab.clusters[k], i, model6.params[0], model6.bse[0], model6.params[1],
        #                                          model.bse[1]))
        #
        #     axs6[row, col].plot(xran6, [f * model6.params[1] + model6.params[0] for f in xran6], '-', alpha=0.8, color=color,
        #                         label=f'{ab.clusters[k]}')
        #     axs6[row, col].fill_between(xran6,
        #                                 [f * (model6.params[1]) + model6.params[0] - 2 * model6.bse[0] for f in xran6],
        #                                 [f * (model6.params[1]) + model6.params[0] + 2 * model6.bse[0] for f in xran6],
        #                                 alpha=0.1, color=color)

        # axs6[row, col].set_title(ab.clusters[k])



        axs[0, k].set_title(ab.clusters[k])
        axs[0, k].set_xlim(xran)
        axs[1, k].set_xlim(xran)
        axs[0, k].set_ylim([0, 2])
        axs[1, k].set_ylim([0, 4])
        axs[0, k].set_xlabel('r - DM(z)')
        axs[0, k].set_xlabel('r - DM(z)')
        axs[0, k].set_ylabel('g - r')
        axs[1, k].set_ylabel('u - r')
        axs[0, k].text(-22, 1.75, f'N={len(r_spec)}')


axs2[0].set_xlim(xran)
axs2[1].set_xlim(xran)
axs2[0].set_ylim(g_yran)
axs2[1].set_ylim(u_yran)
axs2[0].set_xlabel('r - DM(z)')
axs2[1].set_xlabel('r - DM(z)')
axs2[0].set_ylabel('g - r')
axs2[1].set_ylabel('u - r')
axs2[0].legend(loc='lower left', fontsize='large')


# fig3.delaxes(axs3[1, 3])


fig.savefig(ab.plot_dir + f'CMD_merged_spec_{ab.ver}.png')
fig2.savefig(ab.plot_dir + f'CMD_allinone_spec_{ab.ver}.png')
fig3.savefig(ab.plot_dir + f'hist_spec_{ab.ver}.png')
fig4.savefig(ab.plot_dir + f'pps_spec_{ab.ver}.png')
fig5.savefig(ab.plot_dir + f'vel_hist_spec_{ab.ver}.png')
fig6.savefig(ab.plot_dir + f'CMD_each_{ab.ver}_same_err.png')

plt.close(fig)
plt.close(fig2)
# plt.close(fig3)
plt.close(fig4)
plt.close(fig5)
plt.close(fig6)

