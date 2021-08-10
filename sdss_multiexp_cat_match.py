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


work_dir=("/Users/duhokim/work/abell/")
plot_dir=("/Users/duhokim/work/abell/plot/")

clusters = ['A2399', 'A2670']
bands = ['u', 'g', 'r']
stacked = '_stack'
max_sep = 1.0
class_star_lim = 0.9
mag_sys_sex = 'MAG_AUTO'
mag_sys_sdss = 'petroMag'
mag_lim = 99

mw_ext = [[0.159, 0.124, 0.086],
               [0.188, 0.146, 0.101]]

fig1, axs1 = plt.subplots(2, 3, tight_layout = True, figsize = (12, 8))
cmap = plt.cm.rainbow
norm = matplotlib.colors.Normalize(vmin=0, vmax=9000)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

fig2, axs2 = plt.subplots(2, 3, tight_layout = True, figsize = (12, 8))

fig3, axs3 = plt.subplots(2, 2, tight_layout = True, figsize = (8, 8))

alpha3 = 0.1
size3 = 0.5
sdsscol3 = 'orange'
sexcol3 = 'blue'

for i in range(0, len(clusters)):
    merged_cat = ascii.read(work_dir + f'sex/cat/DECam_merged_SEx_cat_{clusters[i]}_Gal_ext_corrected_{ab.ver}.txt')
    merged_coords = SkyCoord(merged_cat['ALPHA_J2000'], merged_cat['DELTA_J2000'], unit='deg')

    sdss_cat = ascii.read(work_dir + 'cat/sdss_galaxy_' + clusters[i] + '.csv')
    sdss_coords = coords.SkyCoord(sdss_cat['ra'], sdss_cat['dec'], unit=(u.deg, u.deg))

    # in_mosaic = np.full(len(sdss_cat), False, dtype=bool)
    # ### Histogram for all SDSS sources inside our r-band mosaic
    # with fits.open(ab.work_dir + 'fits/stacked/' + clusters[i] + '_rsi.fits') as hdu_r:
    #     for j in range(0, len(sdss_cat)):
    #         for jj in range(1, len(hdu_r)):
    #             w = wcs.WCS(hdu_r[jj].header)
    #             cel_coord = [[sdss_cat['ra'][j], sdss_cat['dec'][j]], [0, 0]]
    #             pixcrd = w.wcs_world2pix(cel_coord, 1)
    #             if (pixcrd[0][0] > 0) & (pixcrd[0][0] < hdu_r[jj].shape[1]) & \
    #                     (pixcrd[0][1] > 0) & (pixcrd[0][1] < hdu_r[jj].shape[0]):
    #                 in_mosaic[j] = True
    # f = open('2399_in_mosaic', 'wb')
    # pickle.dump(in_mosaic, f)
    # f.close()

    f = open(f'{clusters[i][1:]}_in_mosaic', 'rb')
    in_mosaic = pickle.load(f)
    f.close()

    idx, d2d, d3d = merged_coords.match_to_catalog_sky(sdss_coords)
    sep_constraint = d2d.arcsec < max_sep
    merged_matches = merged_cat[sep_constraint]
    sdss_matches = sdss_cat[idx[sep_constraint]]

    for j in range(0, len(bands)):
        if j == 1:
            axs1[i, j].set_title(clusters[i])
            axs2[i, j].set_title(clusters[i])
        axs1[i, j].scatter(sdss_matches[mag_sys_sdss+'_'+bands[j]],
                           merged_matches[mag_sys_sex+'_'+bands[j]] + mw_ext[i][j],
                           alpha=0.2, s=1)
        axs1[i, j].plot([0, 30], [0,30], '--', alpha=0.1)
        axs1[i, j].set_xlim([30, 10])
        axs1[i, j].set_ylim([30, 10])
        axs1[i, j].set_xlabel('SDSS Galaxy catalog '+mag_sys_sdss+'_'+bands[j])
        axs1[i, j].set_ylabel('DECam '+mag_sys_sex+'_'+bands[j])

        axs2[i, j].hist(sdss_cat[mag_sys_sdss+'_'+bands[j]][in_mosaic], bins=range(15,31), label='SDSS Galaxy Catalog '+mag_sys_sdss,
                        histtype='step', color='orange')
        axs2[i, j].hist(merged_cat[mag_sys_sex+'_'+bands[j]], bins=range(15,31),
                        label='DECam SEx Catalog '+mag_sys_sex, histtype='step', color='blue')
        axs2[i, j].set_xlabel(bands[j])
        axs2[i, j].set_ylabel('N')
        axs2[i, j].set_xlim([15, 30])
        if i == 0 and j == 0:
            axs2[i, j].legend()

    good_g = (merged_cat[mag_sys_sex+'_r'] < mag_lim) & (
            merged_cat[mag_sys_sex+'_g'] < mag_lim)

    axs3[0, i].scatter(sdss_cat[mag_sys_sdss+'_r'],
                       sdss_cat[mag_sys_sdss+'_g'] - sdss_cat[mag_sys_sdss+'_r'],
                       alpha=alpha3, s=size3, label='SDSS Galaxy Catalog (all)', color=sdsscol3)
    axs3[0, i].scatter(merged_cat[mag_sys_sex+'_r'][good_g],
                       merged_cat[mag_sys_sex+'_g'][good_g] -
                       merged_cat[mag_sys_sex+'_r'][good_g],
                       alpha=alpha3, s=size3, label='DECam (CLASS_STAR < {})'.format(class_star_lim), color=sexcol3)

    good_u = (merged_cat[mag_sys_sex + '_r'] < mag_lim) & (
                merged_cat[mag_sys_sex + '_u'] < mag_lim)

    axs3[1, i].scatter(sdss_cat[mag_sys_sdss+'_r'],
                       sdss_cat[mag_sys_sdss+'_u'] - sdss_cat[mag_sys_sdss+'_r'],
                       alpha=alpha3, s=size3, label='SDSS Galaxy Catalog (all)', color=sdsscol3)
    axs3[1, i].scatter(merged_cat[mag_sys_sex+'_r'][good_u],
                       merged_cat[mag_sys_sex+'_u'][good_u] -
                       merged_cat[mag_sys_sex+'_r'][good_u],
                       alpha=alpha3, s=size3, label='DECam (CLASS_STAR < {})'.format(class_star_lim), color=sexcol3)


    axs3[0, i].set_title(clusters[i])
    axs3[0, i].set_xlim([13, 22])
    axs3[1, i].set_xlim([13, 22])
    axs3[0, i].set_ylim([-1, 2])
    axs3[1, i].set_ylim([-2, 5])
    #axs3[0, i].set_xlabel('r\'')
    axs3[1, i].set_xlabel('r\'')
    # plt.gca().invert_xaxis()
    if i == 0:
        axs3[0, i].legend()
        axs3[0, i].set_ylabel('g\'-r\'')
        axs3[1, i].set_ylabel('u\'-r\'')

# cbar = fig1.colorbar(sm)
# cbar.ax.set_ylabel('Total exposure time [s]')

fig1.savefig(work_dir + f'plot/SDSS_cat_vs_DECam_mag_auto_stacked_1_on_1_{ab.ver}.png')
fig2.savefig(work_dir + f'plot/SDSS_cat_vs_DECam_mag_auto_stacked_hist_{ab.ver}.png')
fig3.savefig(work_dir + f'plot/SDSS_cat_vs_DECam_mag_auto_stacked_cmd_{ab.ver}.png')
