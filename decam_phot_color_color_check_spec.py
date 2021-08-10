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
fig3, axs3 = plt.subplots(1, len(ab.clusters), tight_layout=True, figsize=(20, 3))
fig4, axs4 = plt.subplots(1, len(ab.clusters), tight_layout=True, figsize=(20, 3))
fig5, axs5 = plt.subplots(4, 2, tight_layout=True, figsize=(10, 20))

for k in range(0, len(ab.clusters)):
# for k in range(0, 1):
    xran = [11, 21] - ab.distmod[k]

    sex_cat = ascii.read(ab.sex_dir+f'DECam_merged_SEx_cat_{ab.clusters[k]}_Gal_ext_corrected_{ab.ver}.txt')

    # good_g = (sex_cat['CLASS_STAR'] < 0.2) & (sex_cat[ab.mag_sys+'_g'] < 30) & (sex_cat[ab.mag_sys+'_r'] < xran[1]) & (sex_cat[ab.magerr_sys+'_r'] > 0)
    # good_u = (sex_cat['CLASS_STAR'] < 0.2) & (sex_cat[ab.mag_sys+'_u'] < 30) & (sex_cat[ab.mag_sys+'_r'] < xran[1]) & (sex_cat[ab.magerr_sys+'_r'] > 0)

    ### cross match to spec
    spec = ascii.read(ab.work_dir+f'spec/{ab.clusters[k]}_spec.txt')

    coords_cat = SkyCoord(sex_cat['ALPHA_J2000'], sex_cat['DELTA_J2000'], unit='deg')
    coords_spec = SkyCoord(spec['col2'], spec['col3'], unit='deg')
    idx_spec, d2d, d3d = coords_cat.match_to_catalog_sky(coords_spec)
    matched_cat = (d2d.arcsec < ab.max_sep) & (sex_cat[ab.mag_sys+'_u'] < 900) & \
                  (sex_cat[ab.mag_sys+'_g'] < 900) & (sex_cat[ab.mag_sys+'_r'] < 900)

    axs[0, k].scatter(sex_cat[ab.mag_sys+'_r'][matched_cat] - ab.distmod[k],
                      sex_cat[ab.mag_sys+'_g'][matched_cat] - sex_cat[ab.mag_sys+'_r'][matched_cat],
                      alpha=alpha, s=size)
    axs[1, k].scatter(sex_cat[ab.mag_sys+'_r'][matched_cat] - ab.distmod[k],
                      sex_cat[ab.mag_sys+'_u'][matched_cat] - sex_cat[ab.mag_sys+'_r'][matched_cat],
                      alpha=alpha, s=size)

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
    axs3[k].hist(sex_cat[ab.mag_sys+'_r'], bins=range(11, 25),
                    label='DECam SEx Catalog', histtype='step', color='blue')
    axs3[k].set_ylabel('N_all', color = 'blue')
    # ax2 = axs3[k].twinx()
    # ax2.hist(sex_cat[ab.mag_sys + '_r'][matched_cat], bins=range(11, 31),
    #              label='SDSS Galaxy Catalog',
    #              histtype='step', color='orange')
    # ax2.set_ylabel('N_spec', color='orange')
    axs3[k].hist(sex_cat[ab.mag_sys + '_r'][matched_cat], bins=range(11, 31),
                 label='SDSS Galaxy Catalog',
                 histtype='step', color='orange')
    # ax2.set_ylabel('N_spec', color='orange')
    axs3[k].set_xlabel('r')


    axs3[k].set_xlim([11, 20])
    axs3[k].set_title(ab.clusters[k])


    g_r = np.array(sex_cat[ab.mag_sys + '_g'][matched_cat] - sex_cat[ab.mag_sys + '_r'][matched_cat])
    r_good_g = np.array(sex_cat[ab.mag_sys + '_r'][matched_cat] - ab.distmod[k])
    merr_good_g = np.array(np.sqrt(sex_cat[ab.magerr_sys + '_g'][matched_cat] ** 2 +
                                   sex_cat[ab.magerr_sys + '_r'][matched_cat] ** 2 ))

    u_r = np.array(sex_cat[ab.mag_sys + '_u'][matched_cat] - sex_cat[ab.mag_sys + '_r'][matched_cat])
    r_good_u = np.array(sex_cat[ab.mag_sys + '_r'][matched_cat] - ab.distmod[k])
    merr_good_u = np.array(np.sqrt( sex_cat[ab.magerr_sys + '_u'][matched_cat] ** 2 +
                                    sex_cat[ab.magerr_sys + '_r'][matched_cat] ** 2))

    merr_good_g[np.where(merr_good_g < 0.001)] = 0.001
    merr_good_u[np.where(merr_good_u < 0.001)] = 0.001

    Data_Frame = {'g_r': g_r,
                  'am': r_good_g,
                  'err': merr_good_g
                  }

    Data_Frame2 = {'u_r': u_r,
                  'am': r_good_u,
                  'err': merr_good_u
                  }

    df = pd.DataFrame(Data_Frame, columns=['g_r', 'am', 'err'])
    df2 = pd.DataFrame(Data_Frame2, columns=['u_r', 'am', 'err'])

    X = df[['am']]
    X2 = df2[['am']]
    Y = df['g_r']
    Y2 = df2['u_r']

    # with statsmodels
    X = sm.add_constant(X)
    X2 = sm.add_constant(X2)

    model = sm.WLS(Y, X, weights=1/df['err']).fit()
    model2 = sm.WLS(Y2, X2, weights=1 / df2['err']).fit()

    print(model.summary())
    print(model2.summary())

    print("{} g-r {} {} {} {}".format(ab.clusters[k], model.params[0], model.bse[0],
                                                model.params[1], model.bse[1]))
    print("{} u-r {} {} {} {}".format(ab.clusters[k], model2.params[0], model2.bse[0],
                                                    model2.params[1], model2.bse[1]))

    axs2[0].scatter(sex_cat[ab.mag_sys + '_r'][matched_cat] - ab.distmod[k],
                    sex_cat[ab.mag_sys + '_g'][matched_cat] - sex_cat[ab.mag_sys + '_r'][matched_cat],
                    label=f'{ab.clusters[k]}', alpha=alpha2, s=size2, color=cols[k])
    axs2[0].plot(xran, [f * model.params[1] + model.params[0] for f in xran], '--', alpha=0.8, color=cols[k], label=f'{ab.clusters[k]}')
    axs2[0].fill_between(xran,
                        [f * (model.params[1]) + model.params[0] - model.bse[0] for f in xran],
                        [f * (model.params[1]) + model.params[0] + model.bse[0] for f in xran],
                        alpha=0.1, color=cols[k]
                        )
    axs2[1].scatter(sex_cat[ab.mag_sys + '_r'][matched_cat] - ab.distmod[k],
                    sex_cat[ab.mag_sys + '_u'][matched_cat] - sex_cat[ab.mag_sys + '_r'][matched_cat],
                    alpha=alpha2, s=size2, color=cols[k])
    axs2[1].plot(xran, [f * model2.params[1] + model2.params[0] for f in xran], '--', alpha=0.8, color=cols[k], label=f'{ab.clusters[k]}')
    axs2[1].fill_between(xran,
                         [f * (model2.params[1]) + model2.params[0] - model2.bse[0] for f in xran],
                         [f * (model2.params[1]) + model2.params[0] + model2.bse[0] for f in xran],
                         alpha=0.1, color=cols[k]
                         )

    axs[0, k].set_title(ab.clusters[k])
    axs[0, k].set_xlim(xran)
    axs[1, k].set_xlim(xran)
    axs[0, k].set_ylim([0, 2])
    axs[1, k].set_ylim([0, 4])
    axs[0, k].set_xlabel('r - DM(z)')
    axs[0, k].set_xlabel('r - DM(z)')
    axs[0, k].set_ylabel('g - r')
    axs[1, k].set_ylabel('u - r')
    axs[0, k].text(-22, 1.75, f'N={len(g_r)}')

axs2[0].set_xlim(xran)
axs2[1].set_xlim(xran)
axs2[0].set_ylim(g_yran)
axs2[1].set_ylim(u_yran)
axs2[0].set_xlabel('r - DM(z)')
axs2[1].set_xlabel('r - DM(z)')
axs2[0].set_ylabel('g - r')
axs2[1].set_ylabel('u - r')
axs2[0].legend(loc='lower left', fontsize='large')

fig.savefig(ab.plot_dir + f'CMD_merged_spec_{ab.ver}.png')
fig2.savefig(ab.plot_dir + f'CMD_allinone_spec_{ab.ver}.png')
fig3.savefig(ab.plot_dir + f'hist_spec_{ab.ver}.png')
fig4.savefig(ab.plot_dir + f'pps_spec_{ab.ver}.png')
fig5.savefig(ab.plot_dir + f'vel_hist_spec_{ab.ver}.png')
