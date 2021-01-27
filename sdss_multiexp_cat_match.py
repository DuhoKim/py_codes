from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
from astropy.coordinates import SkyCoord
from astropy import coordinates as coords
import astropy.units as u

work_dir=("/Users/duhokim/work/abell/")
plot_dir=("/Users/duhokim/work/abell/plot/")

clusters = ['A2399', 'A2670', 'A3716']
bands = ['u', 'g', 'r']
stacked = '_stack'
max_sep = 1.0
class_star_lim = 0.9
mag_sys_sex = 'MAG_AUTO'
mag_sys_sdss = 'petroMag'
mag_lim = 99

mw_ext = [[0.159, 0.124, 0.086],
               [0.188, 0.146, 0.101],
               [0.157, 0.122, 0.085]]

fig1, axs1 = plt.subplots(2, 3, tight_layout = True, figsize = (12, 8))
cmap = plt.cm.rainbow
norm = matplotlib.colors.Normalize(vmin=0, vmax=9000)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

fig2, axs2 = plt.subplots(3, 3, tight_layout = True, figsize = (10, 10))

fig3, axs3 = plt.subplots(2, 3, tight_layout = True, figsize = (12, 8))

alpha3 = 0.1
size3 = 0.5
sdsscol3 = 'orange'
sexcol3 = 'blue'

for i in range(0, len(clusters)):
    sex_result_decam = ascii.read(work_dir + 'sex/cat/DECam_19_21_aug_2014_stacked_SEx_cat_' + clusters[i] +
                                  '_match_rad_1as_Gal_ext_corrected_v1.2.txt')
    if i != 2:
        sex_coords_decam = SkyCoord(sex_result_decam['ALPHA_J2000'], sex_result_decam['DELTA_J2000'], unit='deg')

        sdss_galaxy_cat = ascii.read(work_dir + 'cat/sdss_galaxy_' + clusters[i] + '.csv')
        cat_coords = coords.SkyCoord(sdss_galaxy_cat['ra'], sdss_galaxy_cat['dec'], unit=(u.deg, u.deg))

        idx, d2d, d3d = sex_coords_decam.match_to_catalog_sky(cat_coords)
        sep_constraint = d2d.arcsec < max_sep
        sex_matches = sex_result_decam[sep_constraint]
        sdss_matches = sdss_galaxy_cat[idx[sep_constraint]]

        for j in range(0, len(bands)):
            if j == 1:
                axs1[i, j].set_title(clusters[i])
                axs2[i, j].set_title(clusters[i])
            axs1[i, j].scatter(sdss_matches[mag_sys_sdss+'_'+bands[j]],
                               sex_matches[mag_sys_sex+'_'+bands[j]+stacked] + mw_ext[i][j],
                               alpha=0.2, s=1)
            # axs1[i, j].scatter(sdss_matches[mag_sys_sdss + '_' + bands[j]],
            #                    sex_matches[mag_sys_sex + '_' + bands[j] + stacked], alpha=0.2, s=1)
            axs1[i, j].plot([0, 30], [0,30], '--', alpha=0.1)
            axs1[i, j].set_xlim([30, 10])
            axs1[i, j].set_ylim([30, 10])
            axs1[i, j].set_xlabel('SDSS Galaxy catalog '+mag_sys_sdss+'_'+bands[j])
            axs1[i, j].set_ylabel('DECam '+mag_sys_sex+'_'+bands[j])
            # axs1[i, j].gca().invert_xaxis()
            # axs1[i, j].gca().invert_yaxis()

            axs2[i, j].hist(sdss_galaxy_cat[mag_sys_sdss+'_'+bands[j]], bins=20, label='SDSS Galaxy Catalog '+mag_sys_sdss,
                            histtype='step', color='orange')
            axs2[i, j].hist(sex_result_decam[mag_sys_sex+'_'+bands[j]+stacked][sex_result_decam[mag_sys_sex+'_'+bands[j]+stacked] < 100],
                            bins=20, label='DECam SEx Catalog '+mag_sys_sex, histtype='step', color='blue')
            axs2[i, j].set_xlabel(bands[j])
            axs2[i, j].set_ylabel('N')
            axs2[i, j].set_xlim([15, 30])
            if i == 0 and j == 0:
                axs2[i, j].legend()

        good_g = (sex_result_decam[mag_sys_sex+'_r'+stacked] < mag_lim) & (
                sex_result_decam[mag_sys_sex+'_g'+stacked] < mag_lim)

        axs3[0, i].scatter(sdss_galaxy_cat[mag_sys_sdss+'_r'],
                           sdss_galaxy_cat[mag_sys_sdss+'_g'] - sdss_galaxy_cat[mag_sys_sdss+'_r'],
                           alpha=alpha3, s=size3, label='SDSS Galaxy Catalog (all)', color=sdsscol3)
        axs3[0, i].scatter(sex_result_decam[mag_sys_sex+'_r'+stacked][good_g],
                           sex_result_decam[mag_sys_sex+'_g'+stacked][good_g] -
                           sex_result_decam[mag_sys_sex+'_r'+stacked][good_g],
                           alpha=alpha3, s=size3, label='DECam (CLASS_STAR < {})'.format(class_star_lim), color=sexcol3)

        good_u = (sex_result_decam[mag_sys_sex + '_r' + stacked] < mag_lim) & (
                    sex_result_decam[mag_sys_sex + '_u' + stacked] < mag_lim)

        axs3[1, i].scatter(sdss_galaxy_cat[mag_sys_sdss+'_r'],
                           sdss_galaxy_cat[mag_sys_sdss+'_u'] - sdss_galaxy_cat[mag_sys_sdss+'_r'],
                           alpha=alpha3, s=size3, label='SDSS Galaxy Catalog (all)', color=sdsscol3)
        axs3[1, i].scatter(sex_result_decam[mag_sys_sex+'_r'+stacked][good_u],
                           sex_result_decam[mag_sys_sex+'_u'+stacked][good_u] -
                           sex_result_decam[mag_sys_sex+'_r'+stacked][good_u],
                           alpha=alpha3, s=size3, label='DECam (CLASS_STAR < {})'.format(class_star_lim), color=sexcol3)

    else:
        for j in range(0, len(bands)):
            if j == 1:
                axs2[i, j].set_title(clusters[i])
            axs2[i, j].hist(sex_result_decam[mag_sys_sex+'_'+bands[j]+stacked][sex_result_decam[mag_sys_sex+'_'+bands[j]+stacked] < 100],
                            bins=20, label='DECam SEx Catalog '+mag_sys_sex, histtype='step', color='blue')
            axs2[i, j].set_xlabel(bands[j])
            axs2[i, j].set_ylabel('N')
            axs2[i, j].set_xlim([15, 30])

        good_g = (sex_result_decam[mag_sys_sex+'_r'+stacked] < mag_lim) & (
                sex_result_decam[mag_sys_sex+'_g'+stacked] < mag_lim)
        good_u = (sex_result_decam[mag_sys_sex+'_r'+stacked] < mag_lim) & (
                sex_result_decam[mag_sys_sex+'_u'+stacked] < mag_lim)
        axs3[0, i].scatter(sex_result_decam[mag_sys_sex+'_r'+stacked][good_g],
                           sex_result_decam[mag_sys_sex+'_g'+stacked][good_g] -
                           sex_result_decam[mag_sys_sex+'_r'+stacked][good_g],
                           alpha=alpha3, s=size3, label='DECam (CLASS_STAR < {})'.format(class_star_lim), color=sexcol3)
        axs3[1, i].scatter(sex_result_decam[mag_sys_sex+'_r'+stacked][good_u],
                           sex_result_decam[mag_sys_sex+'_u'+stacked][good_u] -
                           sex_result_decam[mag_sys_sex+'_r'+stacked][good_u],
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

fig1.savefig(work_dir + 'plot/SDSS_cat_vs_DECam_mag_auto_stacked_1_on_1_v1.2.png')
fig2.savefig(work_dir + 'plot/SDSS_cat_vs_DECam_mag_auto_stacked_hist_v1.2.png')
fig3.savefig(work_dir + 'plot/SDSS_cat_vs_DECam_mag_auto_stacked_cmd_v1.2.png')
