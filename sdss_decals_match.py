#####################################################
# Python script that match and compare between our cat vs. DECaLS and SDSS
# written by Duho Kim (10 Jun 22)
######################################################
import os
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
from astropy.coordinates import SkyCoord
from astropy import coordinates as coords
import astropy.units as u
import abell_cluster_module as ab
import importlib
importlib.reload(ab)
from astropy import wcs
from astropy.io import fits
import pickle
from astropy.table import Table, vstack

work_dir=("/Users/duhokim/work/abell/")

# clusters = ['A754', 'A2399', 'A2670', 'A3716']
clusters = ['A2399', 'A2670', 'A3716']
bands = ['u', 'g', 'r']

m1_our = [[25.8, 27.3, 26.8],
          [26.4, 27.1, 26.9],
          [26.4, 26.9, 27.1],
          [25.9, 26.8, 26.6]]

m1_SDSS = [[22.5, 23.4, 23.0],
           [22.4, 23.8, 22.9]]

m1_DECaLS = [[26.0, 25.0],
             [25.8, 25.4],
             [26.2, 25.6]]

max_sep = 1.0
class_star_lim = 0.9
mag_sys_sex = 'MAG_AUTO'
mag_sys_sdss = 'petroMag'
mag_lim = 99

alpha = 0.1
size = 1
fs = 20     # fontsize

fig1, axs1 = plt.subplots(3, 2, tight_layout = True, figsize = (7, 8))
fig2, axs2 = plt.subplots(2, 3, tight_layout = True, figsize = (8, 7))

for i in range(0, len(clusters)):
    our_cat = ascii.read(work_dir + f'catalogue/{clusters[i]}_merged_{ab.ver}.txt')
    our_coords = SkyCoord(our_cat['ALPHA_J2000'], our_cat['DELTA_J2000'], unit='deg')

    leg_total = Table()
    for file in os.listdir(f'{work_dir}cat/legacy/{clusters[i]}'):
        if file.endswith('.fits'):
            leg_tract = Table.read(f'{work_dir}cat/legacy/{clusters[i]}/' + file)
            leg_total = vstack([leg_total, leg_tract])

    leg_coords = SkyCoord(leg_total['ra'], leg_total['dec'], unit='deg')

    idx, d2d, d3d = our_coords.match_to_catalog_sky(leg_coords)
    sep_constraint = d2d.arcsec < max_sep
    our_matches = our_cat[sep_constraint]
    leg_matches = leg_total[idx[sep_constraint]]

    for j in range(1, 3):
        is_brighter_than_m1 = ((22.5 - 2.5 * np.log10(leg_matches[f'flux_{bands[j]}'])) < m1_DECaLS[i-1][j-1]) &    \
                    (our_matches[f'MAG_BEST_{bands[j]}'] < m1_our[i][j]) & (leg_matches['type'] != 'PSF')
        axs1[i, j-1].scatter(our_matches[f'MAG_BEST_{bands[j]}'][is_brighter_than_m1],
                           our_matches[f'MAG_BEST_{bands[j]}'][is_brighter_than_m1] -
                           (22.5 - 2.5 * np.log10(leg_matches[f'flux_{bands[j]}'][is_brighter_than_m1])),
                           alpha = alpha, s = size)
        axs1[i, j-1].plot([14, 25], [0, 0], '--', alpha=0.5)
        axs1[i, j-1].plot([m1_DECaLS[i - 1][j-1] - 1, m1_DECaLS[i - 1][j-1] + 1], [-1, 1], ':', alpha=0.5)
        axs1[i, j-1].set_xlim([14, 25])
        axs1[i, j-1].set_ylim([-1, 1])

        axs1[i, j-1].tick_params(direction='in', top=True, right=True)

        if (i == 1) and (j == 1):
            axs1[i, j-1].set_ylabel(r'$m_{our} - m_{DECaLS}$', fontsize=fs)

        if i == 0:
            axs1[i, j-1].text(0.1, 0.8, rf'${ab.bands[j]}\prime$', transform=axs1[i, j-1].transAxes, size=30)

        if j == 2:
            axs1[i, j - 1].yaxis.set_label_position('right')
            axs1[i, j-1].set_ylabel(f'{clusters[i]}', fontsize=fs, rotation=270, labelpad=20)

        if i == 2:
            axs1[i, j-1].set_xlabel(r'$m_{our}$', fontsize=fs)


    if (i == 0) or (i == 1):
        sdss_cat = ascii.read(work_dir + 'cat/sdss_galaxy_' + clusters[i] + '_cModel.csv')
        sdss_coords = coords.SkyCoord(sdss_cat['ra'], sdss_cat['dec'], unit=(u.deg, u.deg))

        idx, d2d, d3d = our_coords.match_to_catalog_sky(sdss_coords)
        sep_constraint = d2d.arcsec < max_sep
        our_matches = our_cat[sep_constraint]
        sdss_matches = sdss_cat[idx[sep_constraint]]

        for j in range(0, 3):
            is_brighter_than_m1 = (sdss_matches[f'cModelMag_{bands[j]}'] < m1_SDSS[i - 1][j]) & \
                                  (our_matches[f'MAG_BEST_{bands[j]}'] < m1_our[i][j])
            axs2[i, j].scatter(our_matches[f'MAG_BEST_{bands[j]}'][is_brighter_than_m1],
                               our_matches[f'MAG_BEST_{bands[j]}'][is_brighter_than_m1] -
                                 sdss_matches[f'cModelMag_{bands[j]}'][is_brighter_than_m1],
                               alpha=alpha, s=size)
            axs2[i, j].plot([14, 25], [0, 0], '--', alpha=0.5)
            axs2[i, j].plot([m1_SDSS[i-1][j]-1, m1_SDSS[i-1][j]+1], [-1, 1], ':', alpha=0.5)
            axs2[i, j].set_xlim([14, 25])
            axs2[i, j].set_ylim([-1, 1])

            axs2[i, j].tick_params(direction='in', top=True, right=True)

            if i == 0:
                axs2[i, j].text(0.1, 0.8, rf'${ab.bands[j]}\prime$', transform=axs2[i, j].transAxes, size=30)

            if j == 0:
                if i == 0:
                    axs2[i, j].set_ylabel(r'$       - m_{SDSS}$', fontsize=fs)
                    axs2[i, j].yaxis.set_label_coords(-0.3, 0.08)
                elif i == 1:
                    axs2[i, j].set_ylabel(r'$m_{our}       $', fontsize=fs)
                    axs2[i, j].yaxis.set_label_coords(-0.3, 0.9)
                # axs2[i, j].set_ylabel(' ', fontsize=fs)
                # axs2[i, j].text(12, 1, r'$m_{our} - m_{SDSS}$', fontsize=fs, rotation=90)
                # axs2[i, j].yaxis.set_label_coords(1.0, 1.5)

            if j == 2:
                axs2[i, j].yaxis.set_label_position('right')
                axs2[i, j].set_ylabel(f'{clusters[i]}', fontsize=fs, rotation=270, labelpad=20)

            if i == 1:
                axs2[i, j].set_xlabel(r'$m_{our}$', fontsize=fs)



fig1.savefig(work_dir + f'plot/comp_DECaLS_{ab.ver}_excl_psf.png')
fig2.savefig(work_dir + f'plot/comp_SDSS_{ab.ver}_model.png')
