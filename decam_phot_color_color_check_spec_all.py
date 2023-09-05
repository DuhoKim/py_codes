from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.lines as mlines
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
from astropy.table import vstack
from scipy.stats import norm
import matplotlib.mlab as mlab
from astropy import units as u
import math
from scipy.optimize import curve_fit
import operator
import os
from os import path
import shutil

def gaus(xx, xx0, sigma):
    return (1/(sigma * np.sqrt(2*np.pi))) * np.exp(-(1/2)*(xx-xx0)**2/sigma**2)

# Define the Gaussian function, B = 1/(2 * sigma^2), sigma = sqrt(1/(2*B))
def Gauss(x, A, B, C):
    y = A*np.exp(-1*(x-B)**2/(2*C**2))
    return y

best_a = -0.039476
best_b = 1.3851
rs_wid = 0.1

size = 15
alpha = 0.5
barWidth = 0.25

size2 = 30
alpha2 = 0.15

g_yran = [0, 1.5]
u_yran = [1, 4]

ins = ['Duho', 'Jake', 'Garreth', 'Adarsh']
cols = ['black', 'red', 'purple', 'darkorange', 'olive', 'darkgreen', 'teal']
img_files = ['_rgb_jar_3.png', '_rgb_sqrt_0_1.png', '_rgb_sqrt_0_01_big.png', '_rgb_sqrt_0_001_big.png', '_gala.png']


### for phase space daigram
classifier = ['B', 'C', 'D', 'E', 'F']
star = mlines.Line2D([], [], color='grey', marker='*', linestyle='None',
                          markersize=10, label='Interacting')
plus = mlines.Line2D([], [], color='grey', marker='P', linestyle='None',
                          markersize=10, label='Post-merger')


params = ['MEAN_MAG_DIFF', 'DATE', 'EXP_TIME', 'SEEING', 'AIRMASS', 'TOTAL_MATCH_NUM']
num_of_param = len(params)

am_cat = pd.read_csv(ab.work_dir+f'sex/cat/airmass.csv')

fig, axs = plt.subplots(2, len(ab.clusters), tight_layout=True, figsize=(20, 6))
fig2, axs2 = plt.subplots(2, 1, tight_layout=True, figsize=(8, 12))
# fig3, axs3 = plt.subplots(2, 4, tight_layout=True, figsize=(12, 6))
fig3 = plt.figure(tight_layout=True, figsize=(14, 6))
fig33 = plt.figure(tight_layout=True, figsize=(12, 6))
fig33_2 = plt.figure(tight_layout=True, figsize=(9, 6))
fig333 = plt.figure(tight_layout=True, figsize=(9, 6))
fig5, axs5 = plt.subplots(2, 3, tight_layout=True, figsize=(15, 10))
# fig6, axs6 = plt.subplots(2, 4, tight_layout=True, figsize=(12, 6))
fig444 = plt.figure(tight_layout=True, figsize=(15, 10))
fig6 = plt.figure(tight_layout=True, figsize=(12, 6))
fig7 = plt.figure(tight_layout=True, figsize=(9, 6))
fig77 = plt.figure(tight_layout=True, figsize=(9, 6))
fig8 = plt.figure(tight_layout=True, figsize=(6, 6))
fig9 = plt.figure(tight_layout=True, figsize=(6, 6))
fig10, axs10 = plt.subplots(tight_layout=True, figsize=(6, 6))
fig11, axs11 = plt.subplots(tight_layout=True, figsize=(6, 6))
fig12 = plt.figure(figsize=(6,6))
fig13 = plt.figure(figsize=(6,6))
fig14 = plt.figure(figsize=(6,6))

gs3 = fig3.add_gridspec(2, 4, wspace=0.2, hspace=0)
gs33 = fig33.add_gridspec(2, 3, wspace=0.2, hspace=0)
gs33_2 = fig33_2.add_gridspec(2, 3, wspace=0.2, hspace=0)
gs333 = fig3.add_gridspec(2, 3, wspace=0.2, hspace=0)
gs = fig6.add_gridspec(2, 4, wspace=0, hspace=0)
gs444 = fig444.add_gridspec(2, 3, wspace=0, hspace=0)
gs7 = fig7.add_gridspec(2, 3, wspace=0, hspace=0)
gs77 = fig77.add_gridspec(2, 3, wspace=0, hspace=0)
axs9 = fig9.add_subplot()
axs12 = fig12.add_subplot()
axs13 = fig13.add_subplot()
axs14 = fig14.add_subplot()

f_int_tot = []
d_int_tot = []
f_pm_tot = []
d_pm_tot = []
f_eith_tot = []
d_eith_tot = []

tot_num = 0
tot = []
dk_int_num = 0
dk_int = []
dk_pm_num = 0
dk_pm = []
gm_int_num = 0
gm_int = []
gm_pm_num = 0
gm_pm = []
jc_int_num = 0
jc_int = []
jc_pm_num = 0
jc_pm = []
ar_int_num = 0
ar_int = []
ar_pm_num = 0
ar_pm = []
vote_int_cum = []   # for histogram of agreement rate in I
vote_pm_cum = []   # for histogram of agreement rate in PM


for k in range(0, len(ab.clusters)):
# for k in range(0, 3):
#     if (k == 4) | (k == 5):
#         continue
    with open(ab.work_dir+f'spec/{ab.clusters[k]}_rs_each_{ab.ver}.radec', 'w') as radec, \
        open(ab.work_dir+f'spec/{ab.clusters[k]}_spec_ind_{ab.ver}.npy', 'wb') as spec_ind,  \
        open(ab.work_dir+f'spec/{ab.clusters[k]}_rs_each_ind_{ab.ver}.npy', 'wb') as rs_ind:
        xran = [11, 21] - ab.distmod[k]

        sex_cat = ascii.read(ab.work_dir+f'catalogue/{ab.clusters[k]}_merged_{ab.ver}.txt')
        vis_cat = ascii.read(ab.work_dir + f'sex/cat/DECam_merged_SEx_cat_{ab.clusters[k]}_Gal_ext_corrected_20rmag_psf.txt')

        sex_cat[ab.mag_sys+'_u'] -= ab.mw_ext[k][0]
        sex_cat[ab.mag_sys + '_g'] -= ab.mw_ext[k][1]
        sex_cat[ab.mag_sys + '_r'] -= ab.mw_ext[k][2]

        vis_cat['MAG_AUTO_r'] -= ab.mw_ext[k][2]

        sex_cat[ab.mag_sys+'_u'] -= kcor.calc_kcor('u', ab.redshifts[k], 'u - r',
                                                    sex_cat[ab.mag_sys+'_u'] - sex_cat[ab.mag_sys+'_r'])
        sex_cat[ab.mag_sys + '_g'] -= kcor.calc_kcor('g', ab.redshifts[k], 'g - r',
                                                     sex_cat[ab.mag_sys + '_g'] - sex_cat[ab.mag_sys + '_r'])
        sex_cat[ab.mag_sys + '_r'] -= kcor.calc_kcor('r', ab.redshifts[k], 'g - r',
                                                     sex_cat[ab.mag_sys + '_g'] - sex_cat[ab.mag_sys + '_r'])

        vis_cat['MAG_AUTO_r'] -= kcor.calc_kcor('r', ab.redshifts[k], 'g - r',
                                                     vis_cat['MAG_AUTO_g'] - vis_cat['MAG_AUTO_r'])

        is_sex_volume = (sex_cat[ab.mag_sys + '_r'] - ab.distmod[k]) < -20
        is_vis_volume = (vis_cat['MAG_AUTO_r'] - ab.distmod[k]) < -20

        ### Volume-limited sample
        sex_cat_orig = sex_cat
        sex_cat = sex_cat[is_sex_volume]
        vis_cat = vis_cat[is_vis_volume]


        ### cross match to spec
        spec = ascii.read(ab.work_dir+f'spec/{ab.clusters[k]}_spec_ra_dec_z_zran_0.05_rran_1.5.txt')

        ### assign row and column number for 2x3 plot row and col
        row_spec = 0 if k < 3 else 1
        col_spec = 1 if k == 6 else k % 3

        if k == 4 or k == 5:
            is_spec = False
        else:
            is_spec = True

        ### radial vel histogram
        if is_spec:
            # best fit of data
            (mu, sigma) = norm.fit(spec['col4'] * const.c / 1e7)

            # fit the A3562 seperately
            if k == 3:
                r_spec = ascii.read(ab.work_dir + f'spec/{ab.clusters[k]}_spec_ra_dec_z_zran_0.05_rran_1.5.txt')
                r_cat = ascii.read(ab.work_dir + f'spec/{ab.clusters[k]}_7_6_VVV_class_uncer_zran_0.05_rran_1.5.csv', format='csv')
                # fit the real A3558 first
                is_3558 = (r_cat['classification'] == 2)
                is_3562 = (r_cat['classification'] == 5)
                (mu58, sigma58) = norm.fit(r_spec['col4'][is_3558] * const.c / 1e7)
                (mu62, sigma62) = norm.fit(r_spec['col4'][is_3562] * const.c / 1e7)
            elif k == 6:
                r_spec = ascii.read(ab.work_dir + f'spec/{ab.clusters[k]}_spec_ra_dec_z_zran_0.05_rran_1.5.txt')
                r_cat = ascii.read(ab.work_dir + f'spec/{ab.clusters[k]}_1_2_EEI_class_uncer_zran_0.05_rran_1.5.csv',
                                   format='csv')
                # fit the real A3558 first
                is_3558 = (r_cat['classification'] == 1)
                is_3562 = (r_cat['classification'] == 2)
                (mu58, sigma58) = norm.fit(r_spec['col4'][is_3558] * const.c / 1e7)
                (mu62, sigma62) = norm.fit(r_spec['col4'][is_3562] * const.c / 1e7)

            # the histogram of the data
            vel = spec['col4'] * const.c / 1e7
            binwidth = mu / 1e2
            n, bins, patches = axs5[row_spec, col_spec].hist(vel, histtype='step',
                                                   bins=np.arange(min(vel), max(vel) + binwidth, binwidth),
                                                   color='black')
            axs5[row_spec, col_spec].set_title(ab.clusters[k], fontsize=20)
            axs5[row_spec, col_spec].set_xlabel(r"Radial Velocity (x10$^4$km/s)", fontsize=20)
            axs5[row_spec, col_spec].set_ylabel("N", fontsize=20)
            axs5[row_spec, col_spec].tick_params(axis='both', which='major', labelsize=16)

            # fit the A3562 seperately
            if k == 3:
                vel_r = r_spec['col4'] * const.c / 1e7
                binwidth58 = mu58 / 1e2
                binwidth62 = mu62 / 1e2
                n58, bins58, patches58 = axs5[row_spec, col_spec].hist(vel_r[is_3558], histtype='step',
                                                                 bins=np.arange(min(vel_r[is_3558]),
                                                                                max(vel_r[is_3558]) + binwidth58,
                                                                                binwidth58),
                                                                 color='green')
                n62, bins62, patches62 = axs5[row_spec, col_spec].hist(vel_r[is_3562], histtype='step',
                                                                 bins=np.arange(min(vel_r[is_3562]),
                                                                                max(vel_r[is_3562]) + binwidth62,
                                                                                binwidth62),
                                                                 color='brown')
            elif k == 6:
                vel_r = r_spec['col4'] * const.c / 1e7
                binwidth58 = mu58 / 1e2
                binwidth62 = mu62 / 1e2
                n58, bins58, patches58 = axs5[row_spec, col_spec].hist(vel_r[is_3558], histtype='step',
                                                                       bins=np.arange(min(vel_r[is_3558]),
                                                                                      max(vel_r[is_3558]) + binwidth58,
                                                                                      binwidth58),
                                                                       color='orange')
                n62, bins62, patches62 = axs5[row_spec, col_spec].hist(vel_r[is_3562], histtype='step',
                                                                       bins=np.arange(min(vel_r[is_3562]),
                                                                                      max(vel_r[is_3562]) + binwidth62,
                                                                                      binwidth62),
                                                                       color='green')

            # add a 'best fit' line
            parameters, covariance = curve_fit(Gauss, bins[:-1] + binwidth / 2, n)
            fit_A = parameters[0]
            fit_B = parameters[1]
            fit_C = np.abs(parameters[2])
            fit_err = np.sqrt(np.diag(covariance))

            fit_y = Gauss(bins[:-1] + binwidth / 2, fit_A, fit_B, fit_C)
            axs5[row_spec, col_spec].plot(bins[:-1] + binwidth / 2, fit_y, '-', color='black')

            axs5[row_spec, col_spec].annotate(r'$\sigma_{rv}$=' + f'{math.ceil(fit_C * 1e4)}km/s',
                                    xy=(0.03, 0.92),
                                    xycoords='axes fraction',
                                    fontsize=20)

            # fit the A3562 seperately
            if k == 3:
                parameters58, covariance58 = curve_fit(Gauss, bins58[:-1] + binwidth58 / 2, n58)
                parameters62, covariance62 = curve_fit(Gauss, bins62[:-1] + binwidth62 / 2, n62)
                fit_A58 = parameters58[0]
                fit_B58 = parameters58[1]
                fit_C58 = np.abs(parameters58[2])
                fit_err58 = np.sqrt(np.diag(covariance58))
                fit_A62 = parameters62[0]
                fit_B62 = parameters62[1]
                fit_C62 = np.abs(parameters62[2])
                fit_err62 = np.sqrt(np.diag(covariance62))

                fit_y58 = Gauss(bins58[:-1] + binwidth58 / 2, fit_A58, fit_B58, fit_C58)
                axs5[row_spec, col_spec].plot(bins58[:-1] + binwidth58 / 2, fit_y58, '-', color='green')

                axs5[row_spec, col_spec].annotate(r'$\sigma_{rv,3558}$=' + f'{math.ceil(fit_C58 * 1e4)}km/s',
                                                  xy=(0.03, 0.82),
                                                  xycoords='axes fraction',
                                                  fontsize=20,
                                                  color='green')
                print(f'fit_A58: {fit_A58}+-{fit_err58[0]}, fit_B58: {fit_B58}+-{fit_err58[1]}, fit_C58: {fit_C58} +- {fit_err58[2]}')
                print(f'n_member58: {sum(is_3558)}')

                fit_y62 = Gauss(bins62[:-1] + binwidth62 / 2, fit_A62, fit_B62, fit_C62)
                axs5[row_spec, col_spec].plot(bins62[:-1] + binwidth62 / 2, fit_y62, '-', color='brown')

                axs5[row_spec, col_spec].annotate(r'$\sigma_{rv,3562}$=' + f'{math.ceil(fit_C62 * 1e4)}km/s',
                                                  xy=(0.03, 0.72),
                                                  xycoords='axes fraction',
                                                  fontsize=20,
                                                  color='brown')
                print(f'fit_A62: {fit_A62}+-{fit_err62[0]}, fit_B62: {fit_B62}+-{fit_err62[1]}, fit_C62: {fit_C62} +- {fit_err62[2]}')
                print(f'n_member62: {sum(is_3562)}')
            elif k == 6:
                parameters58, covariance58 = curve_fit(Gauss, bins58[:-1] + binwidth58 / 2, n58)
                parameters62, covariance62 = curve_fit(Gauss, bins62[:-1] + binwidth62 / 2, n62)
                fit_A58 = parameters58[0]
                fit_B58 = parameters58[1]
                fit_C58 = np.abs(parameters58[2])
                fit_err58 = np.sqrt(np.diag(covariance58))
                fit_A62 = parameters62[0]
                fit_B62 = parameters62[1]
                fit_C62 = np.abs(parameters62[2])
                fit_err62 = np.sqrt(np.diag(covariance62))

                fit_y58 = Gauss(bins58[:-1] + binwidth58 / 2, fit_A58, fit_B58, fit_C58)
                axs5[row_spec, col_spec].plot(bins58[:-1] + binwidth58 / 2, fit_y58, '-', color='orange')

                axs5[row_spec, col_spec].annotate(r'$\sigma_{rv,3716S}$=' + f'{math.ceil(fit_C58 * 1e4)}km/s',
                                                  xy=(0.03, 0.72),
                                                  xycoords='axes fraction',
                                                  fontsize=20,
                                                  color='orange')
                print(
                    f'fit_A58: {fit_A58}+-{fit_err58[0]}, fit_B58: {fit_B58}+-{fit_err58[1]}, fit_C58: {fit_C58} +- {fit_err58[2]}')
                print(f'n_member58: {sum(is_3558)}')

                fit_y62 = Gauss(bins62[:-1] + binwidth62 / 2, fit_A62, fit_B62, fit_C62)
                axs5[row_spec, col_spec].plot(bins62[:-1] + binwidth62 / 2, fit_y62, '-', color='green')

                axs5[row_spec, col_spec].annotate(r'$\sigma_{rv,3716N}$=' + f'{math.ceil(fit_C62 * 1e4)}km/s',
                                                  xy=(0.03, 0.82),
                                                  xycoords='axes fraction',
                                                  fontsize=20,
                                                  color='green')
                print(
                    f'fit_A62: {fit_A62}+-{fit_err62[0]}, fit_B62: {fit_B62}+-{fit_err62[1]}, fit_C62: {fit_C62} +- {fit_err62[2]}')
                print(f'n_member62: {sum(is_3562)}')


            axs5[row_spec, col_spec].axvline(x=fit_B - 3 * fit_C, linestyle='--', color='black', alpha=0.5)
            axs5[row_spec, col_spec].axvline(x=fit_B + 3 * fit_C, linestyle='--', color='black', alpha=0.5)
            axs5[row_spec, col_spec].set_xlim([fit_B - 4 * fit_C, fit_B + 4 * fit_C, ])
            print(f'fit_A: {fit_A}+-{fit_err[0]}, fit_B: {fit_B}+-{fit_err[1]}, fit_C: {fit_C} +- {fit_err[2]}')

            is_member = (vel > fit_B - 3 * fit_C) & (vel < fit_B + 3 * fit_C)
            spec = spec[is_member]

            print(f'orig: {len(vel)}, n_member: {len(spec)}')

            if k == 3:
                sha = ascii.read(ab.work_dir + f'spec/Shapley/catalog.dat')
                sha_vel = sha['col18'] / 1e4

                is_member_sha = (sha_vel > fit_B - 3 * fit_C) & (sha_vel < fit_B + 3 * fit_C)
                sha = sha[is_member_sha]

        ### For volume-limited sample
        coords_cat = SkyCoord(sex_cat['ALPHA_J2000'], sex_cat['DELTA_J2000'], unit='deg')
        coords_vis = SkyCoord(vis_cat['ALPHA_J2000'], vis_cat['DELTA_J2000'], unit='deg')
        coords_spec = SkyCoord(spec['col2'], spec['col3'], unit='deg')
        idx_spec, d2d, d3d = coords_cat.match_to_catalog_sky(coords_spec)
        idx_spec_vis, d2d_vis, d3d = coords_vis.match_to_catalog_sky(coords_spec)
        matched_cat = (d2d.arcsec < ab.max_sep) & (sex_cat[ab.mag_sys+'_u'] < 900) & \
                      (sex_cat[ab.mag_sys+'_g'] < 900) & (sex_cat[ab.mag_sys+'_r'] < 900)
        matched_vis = (d2d_vis.arcsec < ab.max_sep)

        ### For all sample
        coords_orig = SkyCoord(sex_cat_orig['ALPHA_J2000'], sex_cat_orig['DELTA_J2000'], unit='deg')
        coords_spec = SkyCoord(spec['col2'], spec['col3'], unit='deg')
        idx_spec_orig, d2d_orig, d3d = coords_orig.match_to_catalog_sky(coords_spec)
        matched_orig = (d2d_orig.arcsec < ab.max_sep)

        if k == 3:
            coords_sha_str = []
            for j in range(len(sha)):
                coords_sha_str.append(f"{sha['col3'][j]}:{sha['col4'][j]}:{sha['col5'][j]} "
                                      f"{sha['col6'][j]}:{sha['col7'][j]}:{sha['col8'][j]}")
            coords_sha = SkyCoord(coords_sha_str, unit=(u.hourangle, u.deg))
            idx_sha, d2d, d3d = coords_cat.match_to_catalog_sky(coords_sha)
            matched_cat2sha = (d2d.arcsec < ab.max_sep)

        u_spec = sex_cat[ab.mag_sys+'_u'][matched_cat]
        g_spec = sex_cat[ab.mag_sys+'_g'][matched_cat]
        r_spec = sex_cat[ab.mag_sys+'_r'][matched_cat]
        r_spec_abs = r_spec - ab.distmod[k]

        # 1.0 > g-r > 0.5, absolute r_mag cut
        # bright_red = ((sex_cat[ab.mag_sys+'_g'] - sex_cat[ab.mag_sys+'_r']) < 1.0) &    \
        #              ((sex_cat[ab.mag_sys+'_g'] - sex_cat[ab.mag_sys+'_r']) > 0.5) &    \
        #              (sex_cat[ab.mag_sys+'_r'] - ab.distmod[k] < -20)

        # 1.0 > g-r > 0.5, apparent r_mag cut and spec matched
        bright_red = ((sex_cat[ab.mag_sys+'_g'] - sex_cat[ab.mag_sys+'_r']) < 1.0) &    \
                     ((sex_cat[ab.mag_sys+'_g'] - sex_cat[ab.mag_sys+'_r']) > 0.5) &    \
                     (sex_cat[ab.mag_sys+'_r'] < 18) & matched_cat

        g_r_red = sex_cat[ab.mag_sys+'_g'][bright_red] - sex_cat[ab.mag_sys+'_r'][bright_red]
        r_red = sex_cat[ab.mag_sys+'_r'][bright_red]

        axs[0, k].scatter(r_spec_abs, g_spec - r_spec, alpha=alpha, s=size)
        axs[1, k].scatter(r_spec_abs, u_spec - r_spec, alpha=alpha, s=size)

        ### Projected Phase Space diagram
        if is_spec:
            kpc_arcmin = Cosmo.kpc_proper_per_arcmin(ab.redshifts[k])
            sep = coords_cat[matched_cat].separation(ab.coords_cl[k])   # in degrees
            sep_vis = coords_vis[matched_vis].separation(ab.coords_cl[k])  # in degrees
            # sep_vis = coords_vis[matched_vis].separation(coords_x[k])  # in degrees
            sig_z = np.std(spec['col4'])
            xx = (sep * 6e1 * kpc_arcmin / 1e3 / ab.r200[k]).value
            xx_vis = (sep_vis * 6e1 * kpc_arcmin / 1e3 / ab.r200[k]).value

            yy = np.abs(spec['col4'][idx_spec[matched_cat]] * const.c / 1e7 - fit_B) / fit_C
            yy_vis = np.abs(spec['col4'][idx_spec_vis[matched_vis]] * const.c / 1e7 - fit_B) / fit_C

            # yy = np.abs(spec['col4'][idx_spec[matched_cat]] - ab.redshifts[k])/sig_z/(1+ab.redshifts[k])
            # yy_vis = np.abs(spec['col4'][idx_spec_vis[matched_vis]] - ab.redshifts[k]) / sig_z / (1 + ab.redshifts[k])
            # yy_vis2 = np.abs(spec['col4'][idx_spec_vis2[matched_vis2]] - ab.redshifts[k]) / sig_z / (1 + ab.redshifts[k])
            in_vir = yy < (-1.5 / 1.2 * xx + 1.5)
            in_vir_vis = yy_vis < (-1.5 / 1.2 * xx_vis + 1.5)
            # out_vir = yy > (-1.5 / 1.2 * xx + 1.5)
            dist = np.sqrt(yy ** 2 + xx ** 2)
            # dist_vis = np.sqrt(yy_vis ** 2 + xx_vis ** 2)
            dist_vis = xx_vis   # use 2D dist instead of 3D

            n12, bins12, patches12 = axs12.hist(dist, histtype='step', alpha=0.5, color=cols[k], linestyle=':')
            axs12.plot(bins12[:-1], n12/(bins12[:-1]**2), color=cols[k], label=f'{ab.clusters[k]}')
            axs12.set_yscale("log")
            axs12.set_xscale("log")

            if k == 3:
                sep_sha = coords_cat[matched_cat2sha].separation(ab.coords_cl[k])  # in degrees
                sha_z = sha['col18'] / const.c * 1e3
                sig_z = np.std(sha_z)
                xx_sha = (sep_sha * 6e1 * kpc_arcmin / 1e3 / ab.r200[k]).value
                yy_sha = np.abs(sha_z[idx_sha[matched_cat2sha]] * const.c / 1e7 - fit_B) / fit_C
                # yy_sha = np.abs(sha_z[idx_sha[matched_cat2sha]] - ab.redshifts[k]) / sig_z / (1 + ab.redshifts[k])
                in_vir_sha = yy_sha < (-1.5 / 1.2 * xx_sha + 1.5)
                dist_sha = np.sqrt(yy_sha ** 2 + xx_sha ** 2)

            n_three = [0 for x in range(6)]  # number counts for three for in and out of I,PM,other
            d_three = [0 for x in range(3)]  # clustocentric distances for three for I,PM,other

            axs444 = fig444.add_subplot(gs444[row_spec, col_spec])

            # combine visual inspections
            with open(ab.work_dir + f'vis/combined/{ab.clusters[k]}.txt', 'w') as res:
                vis_id = vis_cat['NUMBER'][matched_vis]
                vis1 = ascii.read(ab.work_dir + f'vis/{ins[0]}/{ab.clusters[k]}.vis')
                vis2 = ascii.read(ab.work_dir + f'vis/{ins[1]}/{ab.clusters[k]}.vis')
                vis3 = ascii.read(ab.work_dir + f'vis/{ins[2]}/{ab.clusters[k]}.vis')
                vis4 = ascii.read(ab.work_dir + f'vis/{ins[3]}/{ab.clusters[k]}.vis')
                for i in range(len(vis_id)):
                    if vis_id[i] in vis2['col1']:
                        ind, = np.where(vis2['col1'] == vis_id[i])
                        vote_I = 0
                        vote_PM = 0
                        ind_duho, = np.where(vis1['col1'] == vis_id[i])
                        tot_num += 1
                        tot.append(tot_num)
                        # Duho's vote
                        if vis1['col3'][ind_duho] & 2 ** 7 == 2 ** 7:   # interacting
                            vote_I = vote_I + 1
                            dk_int_num += 1
                        dk_int.append(dk_int_num)
                        if vis1['col3'][ind_duho] & 2 ** 6 == 2 ** 6: # post-merger
                            vote_PM = vote_PM + 1
                            dk_pm_num += 1
                        dk_pm.append(dk_pm_num)
                        # Jake's vote
                        if vis2['col3'][ind] & 2 ** 7 == 2 ** 7:
                            vote_I = vote_I + 1
                            jc_int_num += 1
                        jc_int.append(jc_int_num)
                        if vis2['col3'][ind] & 2 ** 6 == 2 ** 6:
                            vote_PM = vote_PM + 1
                            jc_pm_num += 1
                        jc_pm.append(jc_pm_num)
                        # Garreth's vote
                        if vis3['col3'][ind] & 2 ** 7 == 2 ** 7:
                            vote_I = vote_I + 1
                            gm_int_num += 1
                        gm_int.append(gm_int_num)
                        if vis3['col3'][ind] & 2 ** 6 == 2 ** 6:
                            vote_PM = vote_PM + 1
                            gm_pm_num += 1
                        gm_pm.append(gm_pm_num)
                        # Adarsh's vote
                        if vis4['col3'][ind] & 2 ** 7 == 2 ** 7:
                            vote_I = vote_I + 1
                            ar_int_num += 1
                        ar_int.append(ar_int_num)
                        if vis4['col3'][ind] & 2 ** 6 == 2 ** 6:
                            vote_PM = vote_PM + 1
                            ar_pm_num += 1
                        ar_pm.append(ar_pm_num)

                        vote_int_cum.append(vote_I)
                        vote_pm_cum.append(vote_PM)

                        if vote_I >= 2:       # 2 or more Interacting
                            res.writelines(f"{vis_id[i]} 2 {vote_I}\n")
                            d_three[1] += dist_vis[i]
                            this_alpha = 0.8
                            this_marker = "*"
                            this_size = 50
                            this_col = 'red'
                            if in_vir_vis[i]:
                                n_three[2] += 1
                            else:
                                n_three[3] += 1
                        elif vote_PM >= 2:       # 2 or more PM
                            res.writelines(f"{vis_id[i]} 1 {vote_PM}\n")
                            d_three[2] += dist_vis[i]
                            this_alpha = 0.5
                            this_marker = "P"
                            this_size = 50
                            this_col = 'orange'
                            if in_vir_vis[i]:
                                n_three[4] += 1
                            else:
                                n_three[5] += 1
                        else:
                            res.writelines(f"{vis_id[i]} 0 0\n")
                            d_three[0] += dist_vis[i]
                            this_alpha = 0.1
                            this_marker = "."
                            this_size = 30
                            this_col = 'grey'
                            if in_vir_vis[i]:
                                n_three[0] += 1
                            else:
                                n_three[1] += 1

                        axs444.scatter(xx_vis[i], yy_vis[i], alpha=this_alpha, c=this_col, s=this_size,
                                                         marker=this_marker)

            axs444.plot([0, 1.2], [1.5, 0], linestyle='--')
            axs444.annotate(f'{ab.clusters[k]}', xy=(0.7, 0.9), xycoords='axes fraction', fontsize=30)

            if row_spec == 0 and col_spec < 2:
                plt.setp(axs444.get_xticklabels(), visible=False)
            else:
                plt.setp(axs444.get_xticklabels(), fontsize=15)
                if col_spec == 2:
                    axs444.set_xticks([0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
                axs444.set_xlabel(r'R/R$_{200}$', fontsize=20)

            if col_spec:
                plt.setp(axs444.get_yticklabels(), visible=False)
            else:
                plt.setp(axs444.get_yticklabels(), fontsize=15)
                if row_spec:
                    axs444.set_yticks([0.0, 0.5, 1.0, 1.5, 2.0, 2.5])
                axs444.set_ylabel(r'$\left| \Delta v \right|$/$\sigma$', fontsize=20)
            axs444.set_xlim([0, 3])
            axs444.set_ylim([0, 3])
            axs444.tick_params(direction='in', top=True, right=True)

            num_tot_in = n_three[0] + n_three[2] + n_three[4]
            num_tot_out = n_three[1] + n_three[3] + n_three[5]

            f_int_in = n_three[2] / num_tot_in * 100
            f_int_out = n_three[3] / num_tot_out * 100
            f_pm_in = n_three[4] / num_tot_in * 100
            f_pm_out = n_three[5] / num_tot_out * 100
            f_eith_in = (num_tot_in - n_three[0]) / num_tot_in * 100
            f_eith_out = (num_tot_out - n_three[1]) / num_tot_out * 100

            d_x = []
            d_y = []
            f_y = []

            d_tot = sum(d_three) / sum(n_three)

            if (n_three[2] + n_three[3]):
                d_int = d_three[1] / (n_three[2] + n_three[3])
                d_x.append(1)
                d_y.append(d_int / d_tot)
                d_int_tot.append(d_int / d_tot)
                f_y.append(f_int_out - f_int_in)
                f_int_tot.append(f_int_out - f_int_in)

            if (n_three[4] + n_three[5]):
                d_pm = d_three[2] / (n_three[4] + n_three[5])
                d_x.append(2)
                d_y.append(d_pm / d_tot)
                d_pm_tot.append(d_pm / d_tot)
                f_y.append(f_pm_out - f_pm_in)
                f_pm_tot.append(f_pm_out - f_pm_in)

            if sum(n_three[2:]):
                d_eith = sum(d_three[1:]) / sum(n_three[2:])
                d_x.append(3)
                d_y.append(d_eith / d_tot)
                d_eith_tot.append(d_eith / d_tot)
                f_y.append(f_eith_out - f_eith_in)
                f_eith_tot.append(f_eith_out - f_eith_in)

                # axs10.scatter(d_x, f_y, alpha=0.5, s=ab.m200[k]*2)
                # axs11.scatter(d_x, d_y, alpha=0.5, s=ab.m200[k]*2)
                axs10.scatter(d_x, f_y, alpha=0.5, s=ab.m200[k]*4, marker="o", label=f'{ab.clusters[k]}_A{i+1}')
                axs11.scatter(d_x, d_y, alpha=0.5, s=ab.m200[k]*4, marker="o", label=f'{ab.clusters[k]}_A{i+1}')

                print(f'{ab.clusters[k]} {i}th N_sample:{num_tot_in + num_tot_out}, N_in_vir:{num_tot_in}, '
                      f'N_I:{n_three[2]+n_three[3]}, N_I_in_vir:{n_three[2]}, '
                      f'N_PM:{n_three[4] + n_three[5]}, N_PM_in_vir:{n_three[4]}')


        row = int(k / 4)
        col = k % 4


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

        Data_Frame_red = {'g_r_red': g_r_red, 'r_red': r_red}

        df = pd.DataFrame(Data_Frame, columns=['g_r', 'u_r', 'r_abs', 'r_app', 'g_r_err', 'u_r_err'])
        df_red = pd.DataFrame(Data_Frame_red, columns=['g_r_red', 'r_red'])

        # with statsmodels
        X = sm.add_constant(df[['r_abs']])
        X2 = sm.add_constant(df[['r_app']])
        X3 = sm.add_constant(df_red[['r_red']])

        model = sm.WLS(df['g_r'], X, weights=1/df['g_r_err']).fit()
        model2 = sm.WLS(df['u_r'], X, weights=1 / df['u_r_err']).fit()
        model6 = sm.WLS(df['g_r'], X2, weights=1 / df['g_r_err']).fit()
        model_red = sm.WLS(df_red['g_r_red'], X3).fit()
        one_sig = np.mean(np.abs(model_red.params[1] * r_red + model_red.params[0] - g_r_red))

        print(model_red.summary())
        # print(model2.summary())

        # print("{} g-r {} {} {} {}".format(ab.clusters[k], model6.params[0], model6.bse[0], model6.params[1], model6.bse[1]))
        print("{} g-r {} {} {} {} {}".format(ab.clusters[k], model_red.params[0], model_red.bse[0], model_red.params[1],
                                          model_red.bse[1], one_sig))
        # print("{} u-r {} {} {} {}".format(ab.clusters[k], model2.params[0], model2.bse[0], model2.params[1], model2.bse[1]))

        # axs2[0].scatter(r_spec_abs, g_spec - r_spec, label=f'{ab.clusters[k]}', alpha=alpha2, s=size2, color=cols[k])
        # axs2[0].plot(xran, [f * model.params[1] + model.params[0] for f in xran], '--', alpha=0.8, color=cols[k], label=f'{ab.clusters[k]}')
        # axs2[0].fill_between(xran,
        #                     [f * (model.params[1]) + model.params[0] - model.bse[0] for f in xran],
        #                     [f * (model.params[1]) + model.params[0] + model.bse[0] for f in xran],
        #                     alpha=0.1, color=cols[k])
        # axs2[1].scatter(r_spec_abs, u_spec - r_spec, alpha=alpha2, s=size2, color=cols[k])
        # axs2[1].plot(xran, [f * model2.params[1] + model2.params[0] for f in xran], '--', alpha=0.8, color=cols[k], label=f'{ab.clusters[k]}')
        # axs2[1].fill_between(xran,
        #                      [f * (model2.params[1]) + model2.params[0] - model2.bse[0] for f in xran],
        #                      [f * (model2.params[1]) + model2.params[0] + model2.bse[0] for f in xran],
        #                      alpha=0.1, color=cols[k])

        xran6 = [12, 20]
        size2 = 10
        alpha2 = 0.1

        # axs6 = plt.subplot(gs[row, col])
        axs6 = fig6.add_subplot(gs[row, col])

        if k == 3:
            axs8 = fig8.add_subplot()
            axs8.scatter(r_spec, g_spec - r_spec, label=f'{ab.clusters[k]}', alpha=0.4, s=5, color='teal',
                         facecolors='none')
            axs8.scatter(sex_cat[ab.mag_sys + '_r'], sex_cat[ab.mag_sys + '_g'] - sex_cat[ab.mag_sys + '_r'],
                         label=f'{ab.clusters[k]}', alpha=0.2, s=1, color='grey', facecolors='none')

        axs6.scatter(sex_cat[ab.mag_sys + '_r'], sex_cat[ab.mag_sys + '_g'] - sex_cat[ab.mag_sys + '_r'],
                               label=f'{ab.clusters[k]}', alpha=0.2, s=1, color='grey', facecolors='none')


        if is_spec:
            axs33 = fig33.add_subplot(gs33[row_spec, col_spec])
            axs33_2 = fig33_2.add_subplot(gs33_2[row_spec, col_spec])
            axs333 = fig333.add_subplot(gs333[row_spec, col_spec])
            axs7 = fig7.add_subplot(gs7[row_spec, col_spec])
            axs7.scatter(r_spec, g_spec - r_spec, label=f'{ab.clusters[k]}', alpha=0.4, s=5, color='teal',
                         facecolors='none')
            axs7.scatter(sex_cat[ab.mag_sys + '_r'], sex_cat[ab.mag_sys + '_g'] - sex_cat[ab.mag_sys + '_r'],
                         label=f'{ab.clusters[k]}', alpha=0.2, s=1, color='grey', facecolors='none')

            axs77 = fig77.add_subplot(gs77[row_spec, col_spec])
            axs77.scatter(r_spec, g_spec - r_spec, label=f'{ab.clusters[k]}', alpha=0.4, s=5, color='teal',
                      facecolors='none')
            axs77.scatter(sex_cat[ab.mag_sys + '_r'], sex_cat[ab.mag_sys + '_g'] - sex_cat[ab.mag_sys + '_r'],
                          label=f'{ab.clusters[k]}', alpha=0.2, s=1, color='grey', facecolors='none')
            axs77.axvline(x=-20 + ab.distmod[k], linestyle='--', color='grey')

        #### 1 time fit
        axs6.plot(xran6, [f * model_red.params[1] + model_red.params[0] for f in xran6], '--', alpha=0.8, color='grey',
                     label=f'{ab.clusters[k]}')

        axs6.fill_between(xran6,
                             [f * (model_red.params[1]) + model_red.params[0] - one_sig for f in xran6],
                             [f * (model_red.params[1]) + model_red.params[0] + one_sig for f in xran6],
                             alpha=0.1, color='grey')


        #### BEST FIT
        # axs6.plot(xran6, [f * best_a + best_b for f in xran6], '--', alpha=0.8, color='grey',
        #                     label=f'{ab.clusters[k]}')
        # axs6[row, col].fill_between(xran6,
        #                             [f * best_a + best_b - 0.1 for f in xran6],
        #                             [f * best_a + best_b + 0.1 for f in xran6],
        #                             alpha=0.1, color='grey')
        # axs6.plot(xran6, [f * best_a + best_b - 0.1 for f in xran6], '--', alpha=0.2, color='grey',
        #                     label=f'{ab.clusters[k]}')
        # axs6.plot(xran6, [f * best_a + best_b + 0.1 for f in xran6], '--', alpha=0.2, color='grey',
        #                     label=f'{ab.clusters[k]}')
        axs6.axvline(x=18, linestyle='--', color='grey')

        if is_spec:
            # axs7.axvline(x=18, linestyle=':', color='grey', alpha=0.5)
            axs7.axvline(x=-20 + ab.distmod[k], linestyle='--', color='grey')

        # axs6.axvline(x=-20 + ab.distmod[k], linestyle='--', color='grey')


        if (row == 0 and col < 3):
            plt.setp(axs6.get_xticklabels(), visible=False, fontsize=20)
        else:
            axs6.set_xlabel("r'", fontsize=20)
        if col:
            plt.setp(axs6.get_yticklabels(), visible=False, fontsize=20)
        else:
            axs6.set_ylabel("g' - r'", fontsize=20)
        if is_spec:
            if (row_spec == 0 and col_spec < 2):
                plt.setp(axs33.get_xticklabels(), visible=False, fontsize=20)
                plt.setp(axs33_2.get_xticklabels(), visible=False, fontsize=20)
                plt.setp(axs333.get_xticklabels(), visible=False, fontsize=20)
                plt.setp(axs7.get_xticklabels(), visible=False, fontsize=20)
                plt.setp(axs77.get_xticklabels(), visible=False, fontsize=20)
            else:
                axs33.set_xlabel("r' - DM(z)", fontsize=20)
                axs33_2.set_xlabel("r' - DM(z)", fontsize=20)
                axs333.set_xlabel("r' - DM(z)", fontsize=20)
                axs7.set_xlabel("r'", fontsize=20)
                axs77.set_xlabel("r'", fontsize=20)
            if col_spec:
                # plt.setp(axs33.get_yticklabels(), visible=False, fontsize=20)
                plt.setp(axs7.get_yticklabels(), visible=False, fontsize=20)
                plt.setp(axs77.get_yticklabels(), visible=False, fontsize=20)
            else:
                axs33.set_ylabel("N", fontsize=20)
                axs33_2.set_ylabel(r'N$_{spec}$/N$_{rs}$', fontsize=20)
                axs333.set_ylabel("N", fontsize=20)
                axs7.set_ylabel("g' - r'", fontsize=20)
                axs77.set_ylabel("g' - r'", fontsize=20)
            axs33.set_xlim([-25, -20])
            axs33_2.set_xlim([-25, -20])
            axs7.set_xlim([12, 18])
            axs77.set_xlim([13, 20])
            axs77.set_ylim([0, 1.1])
            axs7.set_ylim([0, 1.1])
            axs33.annotate(f'{ab.clusters[k]}', xy=(0.05, 0.85), xycoords='axes fraction', fontsize=20)
            axs33_2.annotate(f'{ab.clusters[k]}', xy=(0.05, 0.85), xycoords='axes fraction', fontsize=20)
            axs333.annotate(f'{ab.clusters[k]}', xy=(0.05, 0.85), xycoords='axes fraction', fontsize=20)
            axs333.tick_params(direction='in', top=True, right=True)
            axs333.grid(which='both', color='k', linestyle='--', linewidth=0.2, alpha=0.6)
            axs333.set_xticks([-24, -22, -20, -18, -16])
            # axs333.axvline(x=-20, color='k', linestyle='--', line)
            axs77.annotate(f'{ab.clusters[k]}', xy=(0.05, 0.05), xycoords='axes fraction', fontsize=20)
            axs7.annotate(f'{ab.clusters[k]}', xy=(0.05, 0.05), xycoords='axes fraction', fontsize=20)

        axs6.set_xlim([13, 20])
        axs6.set_ylim([0, 1.1])
        axs6.annotate(f'{ab.clusters[k]}', xy=(0.05, 0.05), xycoords='axes fraction', fontsize=20)


        plt.subplots_adjust(hspace=.0, wspace=.0)

        ### red sequence fit
        # factor = 1 if k == 2 else 2
        factor = 1
        slope_fit_lim = 0.001

        # for i in range(0, 10):
        i = 0
        old_slope = 0   # dummy value for initial loop
        while np.abs(old_slope - model_red.params[1]) > slope_fit_lim:
            # color = cols[k] if i == 4 else 'grey'
            i += 1
            color =  'grey'
            old_slope = model_red.params[1]
            in_1sig = (g_r_red < model_red.params[1] * r_red + model_red.params[0] + factor * one_sig) & \
                      (g_r_red > model_red.params[1] * r_red + model_red.params[0] - factor * one_sig)

            New_Data_Frame = {'g_r': np.array(g_r_red[in_1sig]), 'r': np.array(r_red[in_1sig])}
            df = pd.DataFrame(New_Data_Frame, columns=['g_r', 'r'])
            X2 = sm.add_constant(df[['r']])
            model_red = sm.WLS(df['g_r'], X2).fit()
            # one_sig = np.mean(np.abs(model_red.params[1] * r_red + model_red.params[0] - g_r_red))
            one_sig = np.sqrt(np.sum((model_red.params[1] * r_red + model_red.params[0] - g_r_red) ** 2) / len(g_r_red))

            print("{} g-r {} {} {} {} {} {}".format(ab.clusters[k], i, model_red.params[0], model_red.bse[0], model_red.params[1],
                                                 model_red.bse[1], one_sig))

            axs6.plot(xran6, [f * model_red.params[1] + model_red.params[0] for f in xran6], '--', alpha=0.5, color=color,
                                label=f'{ab.clusters[k]}')
            axs6.fill_between(xran6,
                                        [f * (model_red.params[1]) + model_red.params[0] - factor * one_sig for f in xran6],
                                        [f * (model_red.params[1]) + model_red.params[0] + factor * one_sig for f in xran6],
                                        alpha=0.1, color=color)

        if is_spec:
            axs7.plot(xran6, [f * model_red.params[1] + model_red.params[0] for f in xran6], '--', alpha=0.5,
                      color=color,
                      label=f'{ab.clusters[k]}')
            axs7.annotate(rf"g'-r' = r' * {'{:.3f}'.format(model_red.params[1])}", xy=(0.05, 0.93), xycoords='axes fraction', fontsize=10,
                          color='blue')
            axs7.annotate(rf"+ {'{:.3f}'.format(model_red.params[0])}", xy=(0.05, 0.88), xycoords='axes fraction', fontsize=10,
                          color='blue')
        if k != 5:
            axs9.scatter((-20 + ab.distmod[k]) * model_red.params[1] + model_red.params[0], model_red.params[1],
                         s=one_sig*1000, c='black')
            axs9.errorbar((-20 + ab.distmod[k]) * model_red.params[1] + model_red.params[0], model_red.params[1],
                          xerr=model_red.bse[0], yerr=model_red.bse[1], color='black')
            axs9.text((-20 + ab.distmod[k]) * model_red.params[1] + model_red.params[0], model_red.params[1],
                      f'{ab.clusters[k]}', fontsize=20)
            # axs9.text((-20 + ab.distmod[k]) * model_red.params[1] + model_red.params[0], model_red.params[1] + 0.003,
            #           f'z={ab.redshifts[k]}', fontsize=10)
            # axs9.text((-20 + ab.distmod[k]) * model_red.params[1] + model_red.params[0], model_red.params[1] + 0.006,
            #           f'N={len(r_red)}', fontsize=10)
        if is_spec:
            axs7.fill_between(xran6,
                              [f * (model_red.params[1]) + model_red.params[0] - factor * one_sig for f in
                               xran6],
                              [f * (model_red.params[1]) + model_red.params[0] + factor * one_sig for f in
                               xran6],
                              alpha=0.1, color=color)
        print(f'one_sig: {one_sig}')

        ### histogram

        axs3 = fig3.add_subplot(gs3[row, col])
        r_abs_kcor = sex_cat[ab.mag_sys + '_r'] - ab.distmod[k]
        r_abs_kcor_orig = sex_cat_orig[ab.mag_sys + '_r'] - ab.distmod[k]
        num_sam = np.sum(r_abs_kcor < -20)
        num_spec_sam = np.sum(r_abs_kcor[matched_cat] < -20)

        in_rs = (sex_cat[ab.mag_sys + '_g'] - sex_cat[ab.mag_sys + '_r'] < model_red.params[1] * sex_cat[
            ab.mag_sys + '_r'] + model_red.params[0] + one_sig) & \
                (sex_cat[ab.mag_sys + '_g'] - sex_cat[ab.mag_sys + '_r'] > model_red.params[1] * sex_cat[
                    ab.mag_sys + '_r'] + model_red.params[0] - one_sig) & \
                (r_abs_kcor < -20)

        in_rs_spec = (sex_cat[ab.mag_sys + '_g'][matched_cat] - sex_cat[ab.mag_sys + '_r'][matched_cat] < model_red.params[1] * sex_cat[
            ab.mag_sys + '_r'][matched_cat] + model_red.params[0] + one_sig) & \
                (sex_cat[ab.mag_sys + '_g'][matched_cat] - sex_cat[ab.mag_sys + '_r'][matched_cat] > model_red.params[1] * sex_cat[
                    ab.mag_sys + '_r'][matched_cat] + model_red.params[0] - one_sig) & \
                (r_abs_kcor[matched_cat] < -20)

        spec_20 = matched_cat & (r_abs_kcor < -20)

        np.save(rs_ind, in_rs)
        np.save(spec_ind, spec_20)

        num_sam_rs = np.sum(in_rs)
        num_sam_rs_spec = np.sum(in_rs_spec)

        for i in range(0, num_sam_rs):
            radec.writelines(f"{sex_cat['ALPHA_J2000'][in_rs][i]} {sex_cat['DELTA_J2000'][in_rs][i]} \n")

        axs3.hist(r_abs_kcor, range=(-25, -20), histtype='step', color='blue')
        axs3.hist(r_abs_kcor[in_rs], range=(-25, -20), histtype='step', color='red')
        axs3.hist(r_abs_kcor[matched_cat], range=(-25, -20), histtype='step', color='green')

        # axs333.hist(r_abs_kcor_orig, range=(-25, -15), histtype='step', color='blue')


        if is_spec:
            axs333.hist(r_abs_kcor_orig[matched_orig], range=(-25, -15), histtype='step', color='green')
            axs333.annotate(r'N$_{spec}$=' + f'{np.sum(matched_orig)}', xy=(0.05, 0.75), xycoords='axes fraction', fontsize=10,
                          color='green')
            axs333.annotate(f'{ab.specs[k]}', xy=(0.05, 0.65), xycoords='axes fraction', fontsize=10, color='green')
            n_rs, bins_rs, patches_rs = axs33.hist(r_abs_kcor[in_rs], range=(-25, -20), histtype='step', color='red')
            n_spec, bins_spec, patches_spec = axs33.hist(r_abs_kcor[matched_cat][in_rs_spec], range=(-25, -20), histtype='step', color='green')
            is_nonzero = np.where(n_rs > 0)
            axs33_right = axs33.twinx()
            axs33_2.plot(bins_rs[is_nonzero] + 0.45, [a/b for a, b in zip(n_spec[is_nonzero], n_rs[is_nonzero])])
            axs33_right.plot(bins_rs[is_nonzero] + 0.45, [a / b for a, b in zip(n_spec[is_nonzero], n_rs[is_nonzero])],
                             color='blue', alpha=0.5, linestyle=':')

            axs33_right.set_ylim([0, 1])

        # axs3.set_ylabel('N')#, color = 'blue')
        axs3.annotate(r'N$_{phot}$=' + f'{num_sam}', xy=(0.05, 0.75), xycoords='axes fraction', fontsize=10,
                      color='blue')
        axs3.annotate(r'N$_{rs}$=' + f'{num_sam_rs}', xy=(0.05, 0.65), xycoords='axes fraction', fontsize=10,
                      color='red')
        axs3.annotate(r'N$_{spec}$=' + f'{num_spec_sam}', xy=(0.05, 0.55), xycoords='axes fraction',
                      fontsize=10, color='green')

        if is_spec:
            axs33.annotate(r'N$_{rs}$=' + f'{num_sam_rs}', xy=(0.05, 0.65), xycoords='axes fraction', fontsize=15, color='red')
            axs33.annotate(r'N$_{spec}$=' + f'{num_sam_rs_spec}', xy=(0.05, 0.55), xycoords='axes fraction',
                          fontsize=15, color='green')
        # axs3.set_xlabel(r'M$_r$')
        axs3.set_xlim([-25, -20])
        # axs3.set_title(ab.clusters[k])
        axs3.annotate(f'{ab.clusters[k]}', xy=(0.05, 0.85), xycoords='axes fraction', fontsize=20)

        if row == 0 and col < 3:
            plt.setp(axs3.get_xticklabels(), visible=False, fontsize=20)
        else:
            axs3.set_xlabel("r' - DM(z)", fontsize=20)
        if col == 0:
            # plt.setp(axs3.get_yticklabels(), visible=False, fontsize=20)
            # else:
            axs3.set_ylabel("N", fontsize=20)
        if row == 0 and col == 2:
            axs33_right.set_ylabel('Completeness', fontsize=15, color='blue', rotation=270, labelpad=20)
        if row == 1 and col == 2:
            axs33_right.set_ylabel('Completeness', fontsize=15, color='blue', rotation=270, labelpad=20)
        axs33_right.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])

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
    if is_spec:
        axs13.text(tot_num, 30, f'{ab.clusters[k]}', rotation=90)

axs2[0].set_xlim(xran)
axs2[1].set_xlim(xran)
axs2[0].set_ylim(g_yran)
axs2[1].set_ylim(u_yran)
axs2[0].set_xlabel('r - DM(z)')
axs2[1].set_xlabel('r - DM(z)')
axs9.set_xlabel(r'RS intercept at M$_{r\prime}$=-20', fontsize=20)
axs9.set_ylabel(r'RS slope', fontsize=20)
axs9.tick_params(direction='in', top=True, right=True)
axs2[0].set_ylabel('g - r')
axs2[1].set_ylabel('u - r')
axs2[0].legend(loc='lower left', fontsize='large')

# red = mlines.Line2D([], [], color='red', marker='o', linestyle='None',
#                           markersize=10, label='E')
# green = mlines.Line2D([], [], color='green', marker='o', linestyle='None',
#                           markersize=10, label='S0')
# blue = mlines.Line2D([], [], color='blue', marker='o', linestyle='None',
#                           markersize=10, label='Spiral')

axs10.set_ylabel(r'% (out vir) $-$ % (in vir)', fontsize=20)
axs11.set_ylabel(r'd$_{norm}$', fontsize=20)
# plt.setp(axs10, xticks=[1, 2, 3], xticklabels=['Interacting(I)', 'Post-merger(PM)', 'I+PM'])
# plt.setp(axs11, xticks=[1, 2, 3], xticklabels=['Interacting(I)', 'Post-merger(PM)', 'I+PM'])
axs10.set_xticklabels(['I', 'PM', 'I+PM'], fontsize=20)
axs11.set_xticklabels(['I', 'PM', 'I+PM'], fontsize=20)
axs10.tick_params(axis='y', labelsize=15)
axs11.tick_params(axis='y', labelsize=15)
# plt.setp(axs10[1], yticks=[r + barWidth for r in range(5)], yticklabels=['INS1', 'INS2', 'INS3', 'INS4', 'INS5'])
# axs10.legend()
# fig3.delaxes(axs3[1, 3])
axs10.axhline(0, linestyle='--', color='grey')
axs11.axhline(1, linestyle='--', color='grey')
axs10.boxplot([f_int_tot, f_pm_tot, f_eith_tot])
axs11.boxplot([d_int_tot, d_pm_tot, d_eith_tot])
# axs10.legend(ncol=2, loc='upper left', prop={'size': 8})
# axs11.legend(ncol=3, loc='lower right', prop={'size': 7})

axs10.legend(ncol=2, loc='upper left', fontsize=8)
axs11.legend(ncol=3, loc='lower right', fontsize=8)


axs12.set_xlabel(r'$d_{3D}$', fontsize=20)
axs12.set_ylabel(r'$\rho_{gal}$', fontsize=20, labelpad=-5)
axs12.tick_params(axis='both', which='both', direction='in', right=True, labelsize=15, labelright=True)
axs12.text(6, 20, r'N$_{gal}$', fontsize=20, rotation=270)
# axs12.yaxis.set_label_position('right')
# axs12.set_ylabel(f'{ab.clusters[l]}', fontsize=fs, rotation=270, labelpad=20)
xxx = np.arange(0.07, 4, 0.01)
axs12.plot(xxx, 300 / (xxx * (1 + xxx) ** 2), linestyle='--', alpha=0.5, color='grey', label='NFW')
# axs12.text(0.2, 4000 ,r'~NFW', fontsize=20, alpha=0.5)
axs12.plot(xxx, 30 / (xxx * (1 + xxx) ** 2), linestyle='--', alpha=0.5, color='grey')
# axs12.text(4.5, 20, 'N', fontsize=20)
axs12.legend(fontsize=15)

# axs12.plot(xxx, 50 / (xxx * (1 + xxx) ** 2), linestyle='-.', alpha=0.5, color='grey')

axs13.plot(tot, dk_pm, label='A, PM', color='b', linestyle=':')
axs13.plot(tot, dk_int, label='A, M', color='b')
axs13.plot(tot, gm_pm, label='B, PM', color='r', linestyle=':')
axs13.plot(tot, gm_int, label='B, M', color='r')
axs13.plot(tot, jc_pm, label='C, PM', color='g', linestyle=':')
axs13.plot(tot, jc_int, label='C, M', color='g')
axs13.plot(tot, ar_pm, label='D, PM', color='c', linestyle=':')
axs13.plot(tot, ar_int, label='D, M', color='c')
axs13.set_xlabel('Incremental # of gals', fontsize=15)
axs13.set_ylabel('Cumulative # of feature-type gals', fontsize=15)
axs13.tick_params(axis='both', which='both', direction='in', right=True, labelsize=15)
axs13.legend(fontsize=12)

axs14.hist(vote_int_cum, color='b', label='I', bins=3, range=(1, 4), histtype='step')
axs14.hist(vote_pm_cum, color='orange', label='PM', bins=3, range=(1, 4), histtype='step')
axs14.set_xlabel('Vote number', fontsize=20)
axs14.set_ylabel('Counts', fontsize=20)
axs14.set_xticks([1.5, 2.5, 3.5],['1', '2', '3'])
axs14.tick_params(axis='both', which='both', direction='in', labelsize=20)
axs14.legend(fontsize=20)

fig444.tight_layout()

fig5.delaxes(axs5[1][2])

fig.savefig(ab.plot_dir + f'CMD_merged_spec_{ab.ver}_psf.png')
fig2.savefig(ab.plot_dir + f'CMD_allinone_spec_{ab.ver}_psf.png')
fig3.savefig(ab.plot_dir + f'hist_spec_{ab.ver}_psf_each_app_cut.png')
fig33.savefig(ab.plot_dir + f'in_red_sequence_completeness_hist.png')
fig33_2.savefig(ab.plot_dir + f'in_red_sequence_completeness_ratio.png')
fig333.savefig(ab.plot_dir + f'hist_spec_{ab.ver}_no_limit.png')
fig444.savefig(ab.plot_dir + f'pps_spec_{ab.ver}_all.png')
fig5.savefig(ab.plot_dir + f'vel_hist_spec_{ab.ver}_psf.png')
fig6.savefig(ab.plot_dir + f'CMD_each_{ab.ver}_no_err_fit_each_til_001_1sig_rej_psf_proc_app_cut.png')
fig7.savefig(ab.plot_dir + f'CMD_each_{ab.ver}_no_err_fit_each_til_001_1sig_rej_psf_final_with_spec_app_cut.png')
fig77.savefig(ab.plot_dir + f'CMD_each_{ab.ver}_psf_final_with_spec.png')
fig9.savefig(ab.plot_dir + f'rs.png')
fig10.savefig(ab.plot_dir + f'frac_in_vs_out.png')
fig11.savefig(ab.plot_dir + f'clcendist_comp.png')
fig12.savefig(ab.plot_dir + f'radial_number.png')
fig13.savefig(ab.plot_dir + f'class_agree.png')
fig14.savefig(ab.plot_dir + f'hist_agree.png')

plt.close('all')

