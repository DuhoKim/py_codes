from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
from astroquery.sdss import SDSS
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

plot_dir=("/Users/dkim108/Documents/work/plot/")
sex_dir=("/Users/dkim108/Documents/work/sex/gals/")

clusters = ['A2399', 'A2670', 'A3716']
for_gal_ext_u = [0.159, 0.188, 0.157] # irsa.ipac.caltech.edu, S and F (2011)
for_gal_ext_g = [0.124, 0.146, 0.122]
for_gal_ext_r = [0.086, 0.101, 0.085]

class_star_lim = 0.5
max_sep = 1.0   # matching radius limit among each individual exposures

params = ['MAG_ISO', 'MAGERR_ISO', 'FLUX_RADIUS', 'A_WORLD', 'B_WORLD', 'CLASS_STAR', 'TOTAL_EXP_TIME']
num_of_param = len(params)

am_cat = pd.read_csv('/Users/dkim108/Documents/work/cat/airmass.csv')

for k in range(0, len(clusters)):
    fh = open(sex_dir+'DECam_19_21_aug_2014_single_best_exposure_SEx_cat_'+clusters[k]+'_match_rad_1as_'
                                    'Gal_ext_corrected.txt', 'w')

    fh.writelines('#   1 NUMBER                 Running object number \n')
    fh.writelines('#   2 ALPHA_J2000            Right ascension of barycenter (J2000)                      [deg] \n')
    fh.writelines('#   3 DELTA_J2000            Declination of barycenter (J2000)                          [deg] \n')
    fh.writelines('#   4 A_WORLD                Profile RMS along major axis (world units)                 [arcsec] \n')
    fh.writelines('#   5 B_WORLD                Profile RMS along minor axis (world units)                 [arcsec] \n')
    fh.writelines('#   6 CLASS_STAR             S/G classifier output \n')
    fh.writelines('#   7 MAG_ISO_u              Isophotal magnitude in u-band                              [mag] \n')
    fh.writelines('#   8 MAGERR_ISO_u           Isophotal magnitude error in u-band                        [mag] \n')
    fh.writelines('#   9 MAG_ISO_g              Isophotal magnitude in g-band (median)                     [mag] \n')
    fh.writelines('#  10 MAGERR_ISO_g           1-sigma error for isophotal magnitude                      [mag] \n')
    fh.writelines('#  11 MAG_ISO_r              Isophotal magnitude in r-band (median)                     [mag] \n')
    fh.writelines('#  12 MAGERR_ISO_r           1-sigma error for isophotal magnitude                      [mag] \n')

    stk_all = ascii.read(sex_dir + clusters[k] + '_rsi.cat')
    # stk_gal = (stk_all['MAG_ISO'] < 90) & (stk_all['CLASS_STAR'] < class_star_lim) & \
    #             (stk_all['A_WORLD']*3600 > 2) & (stk_all['B_WORLD']*3600 > 2)
    stk_gal = (stk_all['MAG_ISO'] < 90) & (stk_all['CLASS_STAR'] < class_star_lim)
    num_of_gal = len(stk_all[stk_gal])
    stk_coords = SkyCoord(stk_all['ALPHA_J2000'][stk_gal], stk_all['DELTA_J2000'][stk_gal], unit='deg')

    f = listdir(sex_dir+clusters[k])

    u_derived = np.zeros((num_of_gal, num_of_param))
    g_derived = np.zeros((num_of_gal, num_of_param))
    r_derived = np.zeros((num_of_gal, num_of_param))

    u_f, g_f, r_f = [s for s in f if '_ui_' in s], [s for s in f if '_gi_' in s], [s for s in f if '_ri_' in s]

    u_indi_expos = np.empty((num_of_gal, num_of_param, len(u_f)))
    g_indi_expos = np.empty((num_of_gal, num_of_param, len(g_f)))
    r_indi_expos = np.empty((num_of_gal, num_of_param, len(r_f)))

    u_indi_expos[:] = np.nan
    g_indi_expos[:] = np.nan
    r_indi_expos[:] = np.nan

    for i in range(0, len(u_f)):
        time = u_f[i].split('_')[3]
        time_with_colon = time[:2] + ':' + time[2:]
        if '19aug' in u_f[i]:
            a, b = 3.349 - 5.0, -0.318
            am_match = am_cat.loc[(am_cat['date'] == 19) & (am_cat['time'] == time_with_colon)]
        elif '20aug' in u_f[i]:
            a, b = 3.127 - 5.0, -0.054
            am_match = am_cat.loc[(am_cat['date'] == 20) & (am_cat['time'] == time_with_colon)]
        else:
            a, b = 3.200 - 5.0, -0.141
            am_match = am_cat.loc[(am_cat['date'] == 21) & (am_cat['time'] == time_with_colon)]
        X = am_match.iloc[0, 4]
        t = am_match.iloc[0, 3]

        sex_indi = ascii.read(sex_dir+clusters[k]+'/'+u_f[i])
        indi_coords = SkyCoord(sex_indi['ALPHA_J2000'], sex_indi['DELTA_J2000'], unit='deg')

        idx, d2d, d3d = indi_coords.match_to_catalog_sky(stk_coords)
        sep_constraint = (d2d.arcsec < max_sep)
        indi_matches = sex_indi[sep_constraint]
        stk_matches = stk_all[stk_gal][idx[sep_constraint]]

        for j in range(0, num_of_param):
            if j == 0:
                u_indi_expos[idx[sep_constraint], j, i] = indi_matches[params[j]] + 2.5 * np.log10(t) + a + b * X
            elif j == 6:
                u_indi_expos[idx[sep_constraint], j, i] = t
            else:
                u_indi_expos[idx[sep_constraint], j, i] = indi_matches[params[j]]

    for i in range(0, len(g_f)):
        time = g_f[i].split('_')[3]
        time_with_colon = time[:2] + ':' + time[2:]
        if '19aug' in g_f[i]:
            a, b = 2.857 - 2.5, -0.183
            am_match = am_cat.loc[(am_cat['date'] == 19) & (am_cat['time'] == time_with_colon)]
        elif '20aug' in g_f[i]:
            a, b = 2.835 - 2.5, -0.148
            am_match = am_cat.loc[(am_cat['date'] == 20) & (am_cat['time'] == time_with_colon)]
        else:
            a, b = 2.756 - 2.5, -0.095
            am_match = am_cat.loc[(am_cat['date'] == 21) & (am_cat['time'] == time_with_colon)]
        X = am_match.iloc[0, 4]
        t = am_match.iloc[0, 3]

        sex_indi = ascii.read(sex_dir+clusters[k]+'/'+g_f[i])
        indi_coords = SkyCoord(sex_indi['ALPHA_J2000'], sex_indi['DELTA_J2000'], unit='deg')

        idx, d2d, d3d = indi_coords.match_to_catalog_sky(stk_coords)
        sep_constraint = (d2d.arcsec < max_sep)
        indi_matches = sex_indi[sep_constraint]
        stk_matches = stk_all[stk_gal][idx[sep_constraint]]

        for j in range(0, num_of_param):
            if j == 0:
                g_indi_expos[idx[sep_constraint], j, i] = indi_matches[params[j]] + 2.5 * np.log10(t) + a + b * X
            elif j == 6:
                g_indi_expos[idx[sep_constraint], j, i] = t
            else:
                g_indi_expos[idx[sep_constraint], j, i] = indi_matches[params[j]]

    for i in range(0, len(r_f)):
        time = r_f[i].split('_')[3]
        time_with_colon = time[:2] + ':' + time[2:]
        if '19aug' in r_f[i]:
            a, b = 3.024 - 2.5, -0.117
            am_match = am_cat.loc[(am_cat['date'] == 19) & (am_cat['time'] == time_with_colon)]
        elif '20aug' in r_f[i]:
            a, b = 2.949 - 2.5, -0.065
            am_match = am_cat.loc[(am_cat['date'] == 20) & (am_cat['time'] == time_with_colon)]
        else:
            a, b = 3.006 - 2.5, -0.095
            am_match = am_cat.loc[(am_cat['date'] == 21) & (am_cat['time'] == time_with_colon)]
        X = am_match.iloc[0, 4]
        t = am_match.iloc[0, 3]

        sex_indi = ascii.read(sex_dir+clusters[k]+'/'+r_f[i])
        indi_coords = SkyCoord(sex_indi['ALPHA_J2000'], sex_indi['DELTA_J2000'], unit='deg')

        idx, d2d, d3d = indi_coords.match_to_catalog_sky(stk_coords)
        sep_constraint = (d2d.arcsec < max_sep)
        indi_matches = sex_indi[sep_constraint]
        stk_matches = stk_all[stk_gal][idx[sep_constraint]]

        for j in range(0, num_of_param):
            if j == 0:
                r_indi_expos[idx[sep_constraint], j, i] = indi_matches[params[j]] + 2.5 * np.log10(t) + a + b * X
            elif j == 6:
                r_indi_expos[idx[sep_constraint], j, i] = t
            else:
                r_indi_expos[idx[sep_constraint], j, i] = indi_matches[params[j]]

    for j in range(0, num_of_param):
        if j == 0:    # 'MAG_ISO' Galactic extinction correction
            u_derived[:, j] = np.nanmin(u_indi_expos[:, j, :], axis=1) - for_gal_ext_u[k]
            g_derived[:, j] = np.nanmin(g_indi_expos[:, j, :], axis=1) - for_gal_ext_g[k]
            r_derived[:, j] = np.nanmin(r_indi_expos[:, j, :], axis=1) - for_gal_ext_r[k]
        if j == 1:    # 'MAGERR_ISO'
            u_derived[:, j] = np.nanstd(u_indi_expos[:, j - 1, :], axis=1)
            g_derived[:, j] = np.nanstd(g_indi_expos[:, j - 1, :], axis=1)
            r_derived[:, j] = np.nanstd(r_indi_expos[:, j - 1, :], axis=1)
        if j == 6:      # 'TOTAL_EXP_TIME'
            u_derived[:, j] = np.nansum(u_indi_expos[:, j, :], axis=1)
            g_derived[:, j] = np.nansum(g_indi_expos[:, j, :], axis=1)
            r_derived[:, j] = np.nansum(r_indi_expos[:, j, :], axis=1)

    for i in range(0, num_of_gal):
        if np.sum(~np.isnan(u_indi_expos[i, 0, :])) == 1:
            u_derived[i, 1] = u_derived[i, 3] = 0.5
        if np.sum(~np.isnan(g_indi_expos[i, 0, :])) == 1:
            g_derived[i, 1] = g_derived[i, 3] = 0.5
        if np.sum(~np.isnan(r_indi_expos[i, 0, :])) == 1:
            r_derived[i, 1] = r_derived[i, 3] = 0.5
        fh.writelines("{:4}".format(i+1) + ' ' +                                         # NUMBER   \
            "{:12.7f}".format(stk_coords[i].ra.value)+' '+                          # RA       \
            "{:12.7f}".format(stk_coords[i].dec.value) + ' ' +                      # DEC      \
            "{:7.3f}".format(stk_all['A_WORLD'][stk_gal][i]*3600) + ' ' +              # A_WORLD       \
            "{:7.3f}".format(stk_all['B_WORLD'][stk_gal][i]*3600) + ' ' +              # B_WORLD       \
            "{:4.2f}".format(stk_all['CLASS_STAR'][stk_gal][i]) + ' ' +                # CLASS_STAR    \
            "{:7.3f}".format(u_derived[i, 0]) + ' ' +      # u-band 'MAG_ISO'              \
            "{:7.3f}".format(u_derived[i, 1]) + ' ' +      # u-band 'MAGERR_ISO'           \
            # "{:7.3f}".format(u_derived[i, 2]) + ' ' +      # u-band 'MAG_APER'             \
            # "{:7.3f}".format(u_derived[i, 3]) + ' ' +      # u-band 'MAGERR_APER'          \
            "{:3}".format(np.sum(~np.isnan(u_indi_expos[i, 0, :]))) + ' ' +                  # number of u-band exposure     \
            "{:5}".format(u_derived[i, 6]) + ' ' +          # u-band 'TOTAL_EXP_TIME'           \
            "{:7.3f}".format(g_derived[i, 0]) + ' ' +      # g-band 'MAG_ISO'              \
            "{:7.3f}".format(g_derived[i, 1]) + ' ' +      # g-band 'MAGERR_ISO'           \
            # "{:7.3f}".format(g_derived[i, 2]) + ' ' +      # g-band 'MAG_APER'             \
            # "{:7.3f}".format(g_derived[i, 3]) + ' ' +      # g-band 'MAGERR_APER'          \
            "{}".format(np.sum(~np.isnan(g_indi_expos[i, 0, :]))) + ' ' +                  # number of g-band exposure     \
            "{:5}".format(g_derived[i, 6]) + ' ' +          # g-band 'TOTAL_EXP_TIME'           \
            "{:7.3f}".format(r_derived[i, 0]) + ' ' +      # r-band 'MAG_ISO'              \
            "{:7.3f}".format(r_derived[i, 1]) + ' ' +      # r-band 'MAGERR_ISO'           \
            # "{:7.3f}".format(r_derived[i, 2]) + ' ' +      # r-band 'MAG_APER'             \
            # "{:7.3f}".format(r_derived[i, 3]) + ' ' +      # r-band 'MAGERR_APER'          \
            "{:4}".format(np.sum(~np.isnan(r_indi_expos[i, 0, :]))) + ' ' +                    # number of r-band exposure)    \
            "{:5}".format(r_derived[i, 6]) + ' \n'             # r-band 'TOTAL_EXP_TIME'           \
            )

    fh.close()
