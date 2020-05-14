from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy.coordinates import FK4, Angle, Latitude, Longitude
from astropy import coordinates as coords
import os
from os import listdir
from os.path import isfile, join
import astropy.units as u
from astropy.table import Table, vstack
import pandas as pd
from pandas import DataFrame, read_csv
# import matplotlib
# matplotlib.rcParams['text.usetex'] = True
from matplotlib.patches import Circle
import time
import Funcs

plot_dir=("/Users/dkim108/Documents/work/plot/")
sex_dir=("/Users/dkim108/Documents/work/sex/gals/")
cat_dir=("/Users/dkim108/Documents/work/cat/")

# clusters = ['A2399', 'A2670', 'A3716']
clusters = ['A2399', 'A2670']
bands = ['u', 'g', 'r']
date_col = ['red', 'green', 'blue'] # for 19aug, 20aug, 21aug
date_lab = ['19aug', '20aug', '21aug']
# for_gal_ext_u = [0.159, 0.188, 0.157] # irsa.ipac.caltech.edu, S and F (2011)
# for_gal_ext_g = [0.124, 0.146, 0.122]
# for_gal_ext_r = [0.086, 0.101, 0.085]

for_gal_ext_u = [0, 0, 0] # SDSS Galaxy cat didn't correct for Galactic extinction?
for_gal_ext_g = [0, 0, 0]
for_gal_ext_r = [0, 0, 0]

class_star_lim = 0.5
max_sep = 1.0                   # matching radius limit among each individual exposures
mag_lower_lim = [20, 21, 21]          # limit magnitude for analysis
# mag_upper_lim = [14, 14, 13]          # limit magnitude for analysis
mag_upper_lim = [0, 0, 0]          # limit magnitude for analysis
class_lim = [0.5, 0.5, 0.5]     # CLASS_STAR limit for galaxy sample

params = ['MEAN_MAG_DIFF', 'DATE', 'EXP_TIME', 'SEEING', 'AIRMASS', 'TOTAL_MATCH_NUM']
num_of_param = len(params)

am_cat = pd.read_csv('/Users/dkim108/Documents/work/cat/airmass.csv')


for k in range(0, len(clusters)):
    start_time = time.time()

    fig, axs = plt.subplots(3, 3, tight_layout=True, figsize=(12, 8))

    sdss_galaxy_all = ascii.read(cat_dir+'sdss_galaxy_'+clusters[k]+'.csv')


    f = listdir(sex_dir+clusters[k])

    u_f, g_f, r_f = [s for s in f if '_ua_' in s], [s for s in f if '_ga_' in s], [s for s in f if '_ra_' in s]

    u_d = np.empty((num_of_param, len(u_f)))
    g_d = np.empty((num_of_param, len(g_f)))
    r_d = np.empty((num_of_param, len(r_f)))

    u_d[:] = np.nan
    g_d[:] = np.nan
    r_d[:] = np.nan

    for i in range(0, len(u_f)):
        time_str = u_f[i].split('_')[3]
        time_with_colon = time_str[:2] + ':' + time_str[2:]
        if '19aug' in u_f[i]:
            # a, b = 3.349 - 5.0, -0.318   # -1.651
            a, b = -1.642, -0.303
            am_match = am_cat.loc[(am_cat['date'] == 19) & (am_cat['time'] == time_with_colon)]
            u_d[1, i] = 19
            col = 0
        elif '20aug' in u_f[i]:
            # a, b = 3.127 - 5.0, -0.054    # -1.873
            a, b = -1.257, -0.592
            am_match = am_cat.loc[(am_cat['date'] == 20) & (am_cat['time'] == time_with_colon)]
            u_d[1, i] = 20
            col = 1
        else:
            # a, b = 3.200 - 5.0, -0.141    # -1.8
            a, b = -1.796, -0.132
            am_match = am_cat.loc[(am_cat['date'] == 21) & (am_cat['time'] == time_with_colon)]
            u_d[1, i] = 21
            col = 2
        X = am_match.iloc[0, 4]
        u_d[4, i] = X                   # airmass
        u_d[3, i] = am_match.iloc[0, 8] # seeing
        t = am_match.iloc[0, 3]
        u_d[2, i] = t

        sdss_galaxy_flag = (sdss_galaxy_all['petroMag_u'] < mag_lower_lim[0]) & \
                           (sdss_galaxy_all['petroMag_u'] > mag_upper_lim[0])
        sdss_galaxy_cat = sdss_galaxy_all[sdss_galaxy_flag]
        cat_coords = coords.SkyCoord(sdss_galaxy_cat['ra'], sdss_galaxy_cat['dec'], unit=(u.deg, u.deg))

        sex_all = ascii.read(sex_dir+clusters[k]+'/'+u_f[i])
        sex_flag = (sex_all['FLAGS'] == 0) & (sex_all['MAG_PETRO'] + 2.5 * np.log10(t) + a + b * X < mag_lower_lim[0]) \
                   & (sex_all['MAG_PETRO'] + 2.5 * np.log10(t) + a + b * X > mag_upper_lim[0]) \
                   & (sex_all['CLASS_STAR'] < class_lim[0])
        sex_indi = sex_all[sex_flag]
        indi_coords = SkyCoord(sex_indi['ALPHA_J2000'], sex_indi['DELTA_J2000'], unit='deg')

        idx, d2d, d3d = indi_coords.match_to_catalog_sky(cat_coords)
        sep_constraint = (d2d.arcsec < max_sep)
        indi_matches = sex_indi[sep_constraint]
        stk_matches = sdss_galaxy_cat[idx[sep_constraint]]

        u_d[5, i] = len(indi_matches)
        # u_d[0, i] = np.sum(np.abs(indi_matches['MAG_ISO']+ 2.5 * np.log10(t) + a + b * X - stk_matches['cModelMag_u']))/len(indi_matches)
        std_mags = indi_matches['MAG_PETRO'] + 2.5 * np.log10(t) + a + b * X
        u_d[0, i] = np.mean(stk_matches['petroMag_u'] - std_mags)
        axs[2, 0].errorbar(stk_matches['petroMag_u'],
                           stk_matches['petroMag_u'] - std_mags,
                           yerr = np.sqrt(indi_matches['MAGERR_PETRO']**2 + stk_matches['petroMagErr_u']**2),
                           xerr = stk_matches['petroMagErr_u'],
                           fmt = 'o',
                           alpha=0.5,
                           ms=1,
                           color=date_col[col])
    for i in range(0, len(g_f)):
        time_str = g_f[i].split('_')[3]
        time_with_colon = time_str[:2] + ':' + time_str[2:]
        if '19aug' in g_f[i]:
            # a, b = 2.857 - 2.5, -0.183    # 0.357
            a, b = 0.381, -0.186
            am_match = am_cat.loc[(am_cat['date'] == 19) & (am_cat['time'] == time_with_colon)]
            g_d[1, i] = 19
            col = 0
        elif '20aug' in g_f[i]:
            # a, b = 2.835 - 2.5, -0.148    # 0.335
            a, b = 0.361, -0.163
            am_match = am_cat.loc[(am_cat['date'] == 20) & (am_cat['time'] == time_with_colon)]
            g_d[1, i] = 19
            col = 1
        else:
            # a, b = 2.756 - 2.5, -0.095  # 0.256
            a, b = 0.179, -0.036
            am_match = am_cat.loc[(am_cat['date'] == 21) & (am_cat['time'] == time_with_colon)]
            g_d[1, i] = 19
            col = 2
        X = am_match.iloc[0, 4]
        g_d[4, i] = X
        g_d[3, i] = am_match.iloc[0, 8]
        t = am_match.iloc[0, 3]
        g_d[2, i] = t

        sdss_galaxy_flag = (sdss_galaxy_all['petroMag_g'] < mag_lower_lim[1]) & \
                           (sdss_galaxy_all['petroMag_g'] > mag_upper_lim[1])
        sdss_galaxy_cat = sdss_galaxy_all[sdss_galaxy_flag]
        cat_coords = coords.SkyCoord(sdss_galaxy_cat['ra'], sdss_galaxy_cat['dec'], unit=(u.deg, u.deg))

        sex_all = ascii.read(sex_dir+clusters[k]+'/'+g_f[i])
        sex_flag = (sex_all['FLAGS'] == 0) & (sex_all['MAG_PETRO'] + 2.5 * np.log10(t) + a + b * X < mag_lower_lim[1]) \
                   & (sex_all['MAG_PETRO'] + 2.5 * np.log10(t) + a + b * X > mag_upper_lim[1]) \
                   & (sex_all['CLASS_STAR'] < class_lim[1])
        sex_indi = sex_all[sex_flag]
        indi_coords = SkyCoord(sex_indi['ALPHA_J2000'], sex_indi['DELTA_J2000'], unit='deg')

        idx, d2d, d3d = indi_coords.match_to_catalog_sky(cat_coords)
        sep_constraint = (d2d.arcsec < max_sep)
        indi_matches = sex_indi[sep_constraint]
        stk_matches = sdss_galaxy_cat[idx[sep_constraint]]

        g_d[5, i] = len(indi_matches)
        # g_d[0, i] = np.sum(np.abs(indi_matches['MAG_ISO']+ 2.5 * np.log10(t) + a + b * X - stk_matches['cModelMag_g']))/len(indi_matches)
        std_mags = indi_matches['MAG_PETRO'] + 2.5 * np.log10(t) + a + b * X
        g_d[0, i] = np.mean(stk_matches['petroMag_g'] - std_mags)
        axs[2, 1].errorbar(stk_matches['petroMag_g'],
                           stk_matches['petroMag_g'] - std_mags,
                           yerr=np.sqrt(indi_matches['MAGERR_PETRO'] ** 2 + stk_matches['petroMagErr_g'] ** 2),
                           xerr=stk_matches['petroMagErr_g'],
                           fmt='o',
                           alpha=0.5,
                           ms=1,
                           color = date_col[col])
    for i in range(0, len(r_f)):
        time_str = r_f[i].split('_')[3]
        time_with_colon = time_str[:2] + ':' + time_str[2:]
        if '19aug' in r_f[i]:
            # a, b = 3.024 - 2.5, -0.117  # 0.524
            a, b = 0.567, -0.136
            am_match = am_cat.loc[(am_cat['date'] == 19) & (am_cat['time'] == time_with_colon)]
            r_d[1, i] = 19
            col = 0
        elif '20aug' in r_f[i]:
            # a, b = 2.949 - 2.5, -0.065  # 0.449
            a, b = 0.411, -0.047
            am_match = am_cat.loc[(am_cat['date'] == 20) & (am_cat['time'] == time_with_colon)]
            r_d[1, i] = 20
            col = 1
        else:
            # a, b = 3.006 - 2.5, -0.095  # 0.506
            a, b = 0.492, -0.078
            am_match = am_cat.loc[(am_cat['date'] == 21) & (am_cat['time'] == time_with_colon)]
            r_d[1, i] = 21
            col = 2
        X = am_match.iloc[0, 4]
        r_d[4, i] = X
        r_d[3, i] = am_match.iloc[0, 8]
        t = am_match.iloc[0, 3]
        r_d[2, i] = t

        sdss_galaxy_flag = (sdss_galaxy_all['petroMag_r'] < mag_lower_lim[2]) & \
                           (sdss_galaxy_all['petroMag_r'] > mag_upper_lim[2])
        sdss_galaxy_cat = sdss_galaxy_all[sdss_galaxy_flag]
        cat_coords = coords.SkyCoord(sdss_galaxy_cat['ra'], sdss_galaxy_cat['dec'], unit=(u.deg, u.deg))

        sex_all = ascii.read(sex_dir+clusters[k]+'/'+r_f[i])
        sex_flag = (sex_all['FLAGS'] == 0) & (sex_all['MAG_PETRO'] + 2.5 * np.log10(t) + a + b * X < mag_lower_lim[2]) \
                   & (sex_all['MAG_PETRO'] + 2.5 * np.log10(t) + a + b * X > mag_upper_lim[2]) \
                   & (sex_all['CLASS_STAR'] < class_lim[2])
        sex_indi = sex_all[sex_flag]
        indi_coords = SkyCoord(sex_indi['ALPHA_J2000'], sex_indi['DELTA_J2000'], unit='deg')

        idx, d2d, d3d = indi_coords.match_to_catalog_sky(cat_coords)
        sep_constraint = (d2d.arcsec < max_sep)
        indi_matches = sex_indi[sep_constraint]
        stk_matches = sdss_galaxy_cat[idx[sep_constraint]]

        r_d[5, i] = len(indi_matches)
        # r_d[0, i] = np.sum(np.abs(indi_matches['MAG_ISO'] + 2.5 * np.log10(t) + a + b * X  - stk_matches['cModelMag_r']))/len(indi_matches)
        std_mags = indi_matches['MAG_PETRO'] + 2.5 * np.log10(t) + a + b * X
        r_d[0, i] = np.mean(stk_matches['petroMag_r'] - std_mags)
        axs[2, 2].errorbar(stk_matches['petroMag_r'],
                           stk_matches['petroMag_r'] - std_mags,
                           yerr=np.sqrt(indi_matches['MAGERR_PETRO'] ** 2 + stk_matches['petroMagErr_r'] ** 2),
                           xerr=stk_matches['petroMagErr_r'],
                           fmt='o',
                           alpha=0.5,
                           ms=1,
                           color = date_col[col])
    axs[2, 0].invert_xaxis()
    axs[2, 0].set_xlabel('petroMag_u')
    axs[2, 0].set_ylabel('petroMag_u - MAG_PETRO')
    axs[2, 0].plot([13, 21], [0, 0], linestyle=':', alpha=0.5)
    axs[2, 1].invert_xaxis()
    axs[2, 1].set_xlabel('petroMag_g')
    axs[2, 1].set_ylabel('petroMag_g - MAG_PETRO')
    axs[2, 1].plot([13, 21], [0, 0], linestyle=':', alpha=0.5)
    axs[2, 2].invert_xaxis()
    axs[2, 2].set_xlabel('petroMag_r')
    axs[2, 2].set_ylabel('petroMag_r - MAG_PETRO')
    axs[2, 2].plot([13, 21], [0, 0], linestyle=':', alpha=0.5)
    for j in range(0, len(bands)):
        if j == 0:
            plot_data = u_d
        if j == 1:
            plot_data = g_d
            axs[0, j].set_title(clusters[k])
        if j == 2:
            plot_data = r_d
        size = plot_data[2, :]/30
        color = np.empty(len(plot_data[0, :]), int)
        label = ["" for x in range(len(plot_data[0, :]))]
        for i in range(0, len(plot_data[0, :])):
            color[i] = np.int(plot_data[1, i]-19)
            label[i] = date_lab[np.int(plot_data[1, i]-19)]
        scatter = axs[0, j].scatter(plot_data[0, :], plot_data[3, :], alpha=0.5, s=size.astype(int), c=color,
                                    label=label)
        # axs1[i, j].plot([0, 30], [0, 30], '--', alpha=0.1)
        # axs1[i, j].set_xlim([30, 10])
        # axs1[i, j].set_ylim([30, 10])
        axs[0, j].set_xlabel('<SDSS Galaxy petroMag_' + bands[j] + '-MAG_PETRO>')
        axs[0, j].set_ylabel('seeing')

        # handles, labels = scatter.legend_elements(prop="colors", alpha=0.6)
        # kw = dict(prop="colors", num=3)
        # legend1 = axs[0, j].legend(*scatter.legend_elements(prop="colors"), title='date (0:19aug, 1:20, 2:21)', loc='lower right')
        # axs[0, j].add_artist(legend1)

        # legend1 = axs[0, j].legend(title='date (0:19aug, 1:20, 2:21)',
        #                            loc='lower right')
        # axs[0, j].add_artist(legend1)

        # handles, labels = scatter.legend_elements(prop="sizes", alpha=0.6)
        # legend2 = axs[0, j].legend(loc="upper right", title="Exp t(x30s)")

        # kw = dict(prop="sizes", num=3, color=)
            # custom_date = [Circle([0], [0], color='red', radius=5),
            #                   Circle([0], [0], color='green',  radius=5),
            #                   Circle([0], [0],  color='blue', radius=5)]
            # axs[0, j].legend(custom_date, ['19aug', '20aug', '21aug'])
        # elif j == 1:
            # custom_time = [Circle([0], [0], color='green', radius=1),
            #                Circle([0], [0], color='green', radius=5)]
            # axs[0, j].legend(custom_time, ['60s', '300s'])

        axs[1, j].scatter(plot_data[0, :], plot_data[4, :], alpha=0.5, s=size.astype(int), c=color)
        # axs[1, j].set_xlabel(bands[j])
        axs[1, j].set_ylabel('Airmass')


    fig.savefig(plot_dir + clusters[k]+'_SDSS_cat_vs_DECam_standardized_find_best_single_exposure_aper_corr_mean_lt_20mag_petro_no_upper.png')

    print("--- %s minutes ---" % ((time.time() - start_time)/60.0))

    Funcs.print_time()