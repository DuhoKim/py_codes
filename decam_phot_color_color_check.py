from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import abell_cluster_module as ab
import statsmodels.api as sm
import importlib
importlib.reload(ab)

size = 1
alpha = 0.2

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

for k in range(0, len(ab.clusters)):
# for k in range(0, 2):
    xran = [11, 21] #- ab.distmod[k]

    sex_cat = ascii.read(ab.sex_dir+f'DECam_merged_SEx_cat_{ab.clusters[k]}_Gal_ext_corrected_{ab.ver}.txt')

    good_g = (sex_cat[ab.mag_sys+'_g'] < 30) #& (sex_cat[ab.mag_sys+'_r'] < xran[1])
    good_u = (sex_cat[ab.mag_sys+'_u'] < 30) #& (sex_cat[ab.mag_sys+'_r'] < xran[1])

    axs[0, k].scatter(sex_cat[ab.mag_sys+'_r'][good_g],# - ab.distmod[k],
                      sex_cat[ab.mag_sys+'_g'][good_g] - sex_cat[ab.mag_sys+'_r'][good_g],
                      alpha=alpha, s=size)
    axs[1, k].scatter(sex_cat[ab.mag_sys+'_r'][good_u],#- ab.distmod[k],
                      sex_cat[ab.mag_sys+'_u'][good_u] - sex_cat[ab.mag_sys+'_r'][good_u],
                      alpha=alpha, s=size)

    g_r = np.array(sex_cat[ab.mag_sys + '_g'][good_g] - sex_cat[ab.mag_sys + '_r'][good_g])
    r_good_g = np.array(sex_cat[ab.mag_sys + '_r'][good_g] - ab.distmod[k])
    merr_good_g = np.array(sex_cat[ab.magerr_sys + '_r'][good_g])

    u_r = np.array(sex_cat[ab.mag_sys + '_u'][good_u] - sex_cat[ab.mag_sys + '_r'][good_u])
    r_good_u = np.array(sex_cat[ab.mag_sys + '_r'][good_u] - ab.distmod[k])
    merr_good_u = np.array(sex_cat[ab.magerr_sys + '_r'][good_u])

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

    print("{} g-r {} {} {} {}".format(ab.clusters[k], model.params[0], model.bse[0],
                                                model.params[1], model.bse[1]))
    print("{} u-r {} {} {} {}".format(ab.clusters[k], model2.params[0], model2.bse[0],
                                                    model2.params[1], model2.bse[1]))

    if k == 3 or k == 4:
        axs2[0].scatter(sex_cat[ab.mag_sys + '_r'][good_g] - ab.distmod[k],
                           sex_cat[ab.mag_sys + '_g'][good_g] - sex_cat[ab.mag_sys + '_r'][good_g],
                           label=f'{ab.clusters[k]}', alpha=alpha, s=size, color=cols[k])
        axs2[0].plot(xran, [f * model.params[1] + model.params[0] for f in xran], '--', alpha=0.2, color=cols[k], label=f'{ab.clusters[k]}')
        # axs2[0].fill_between(xran,
        #                     [f * (model.params[1]+model.bse[1]) + model.params[0] for f in xran],
        #                     [f * (model.params[1]-model.bse[1]) + model.params[0] for f in xran],
        #                     alpha=0.1, color=cols[k]
        #                     )
        axs2[1].scatter(sex_cat[ab.mag_sys + '_r'][good_u] - ab.distmod[k],
                           sex_cat[ab.mag_sys + '_u'][good_u] - sex_cat[ab.mag_sys + '_r'][good_u],
                           alpha=alpha, s=size, color=cols[k])
        axs2[1].plot(xran, [f * model2.params[1] + model2.params[0] for f in xran], '--', alpha=0.2, color=cols[k], label=f'{ab.clusters[k]}')
        # axs2[1].fill_between(xran,
        #                     [f * (model2.params[1]+model2.bse[1]) + model2.params[0] for f in xran],
        #                     [f * (model2.params[1]-model2.bse[1]) + model2.params[0] for f in xran],
        #                     alpha=0.1, color=cols[k]
        #                     )
    else:
        axs2[0].scatter(sex_cat[ab.mag_sys + '_r'][good_g] - ab.distmod[k],
                        sex_cat[ab.mag_sys + '_g'][good_g] - sex_cat[ab.mag_sys + '_r'][good_g],
                        label=f'{ab.clusters[k]}', alpha=alpha2, s=size2, color=cols[k])
        axs2[0].plot(xran, [f * model.params[1] + model.params[0] for f in xran], '--', alpha=0.8, color=cols[k], label=f'{ab.clusters[k]}')
        # axs2[0].fill_between(xran,
        #                     [f * (model.params[1]+model.bse[1]) + model.params[0] for f in xran],
        #                     [f * (model.params[1]-model.bse[1]) + model.params[0] for f in xran],
        #                     alpha=0.1, color=cols[k]
        #                     )
        axs2[1].scatter(sex_cat[ab.mag_sys + '_r'][good_u] - ab.distmod[k],
                        sex_cat[ab.mag_sys + '_u'][good_u] - sex_cat[ab.mag_sys + '_r'][good_u],
                        alpha=alpha2, s=size2, color=cols[k])
        axs2[1].plot(xran, [f * model2.params[1] + model2.params[0] for f in xran], '--', alpha=0.8, color=cols[k], label=f'{ab.clusters[k]}')
        # axs2[1].fill_between(xran,
        #                      [f * (model2.params[1] + model2.bse[1]) + model2.params[0] for f in xran],
        #                      [f * (model2.params[1] - model2.bse[1]) + model2.params[0] for f in xran],
        #                      alpha=0.1, color=cols[k]
        #                      )


    # axs2[0].contour([sex_cat[ab.mag_sys + '_r'][galaxy],
    #                    sex_cat[ab.mag_sys + '_g'][galaxy] - sex_cat[ab.mag_sys + '_r'][galaxy]],
    #                    label=f'{ab.clusters[k]}')
    # axs2[1].contour([sex_cat[ab.mag_sys + '_r'][galaxy],
    #                    sex_cat[ab.mag_sys + '_u'][galaxy] - sex_cat[ab.mag_sys + '_r'][galaxy]])


    axs[0, k].set_title(ab.clusters[k])
    axs[0, k].set_xlim(xran)
    axs[1, k].set_xlim(xran)
    axs[0, k].set_ylim([0, 2])
    axs[1, k].set_ylim([0, 4])
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

axs2[0].set_xlim(xran)
axs2[1].set_xlim(xran)
axs2[0].set_ylim(g_yran)
axs2[1].set_ylim(u_yran)
axs2[0].set_xlabel('r - DM(z)')
axs2[1].set_xlabel('r - DM(z)')
axs2[0].set_ylabel('g - r')
axs2[1].set_ylabel('u - r')
axs2[0].legend(loc='lower left', fontsize='large')

fig.savefig(ab.plot_dir + 'CMD_merged_'+ab.ver+'.png')
fig2.savefig(ab.plot_dir + 'CMD_allinone_'+ab.ver+'.png')
