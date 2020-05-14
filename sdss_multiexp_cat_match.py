from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
from astropy.coordinates import SkyCoord
from astropy import coordinates as coords
import astropy.units as u

data_dir=("/Users/dkim108/Documents/work/sex/gals/")
plot_dir=("/Users/dkim108/Documents/work/plot/")
cat_dir=("/Users/dkim108/Documents/work/cat/")

clusters = ['A2399', 'A2670', 'A3716']
bands = ['u', 'g', 'r']
max_sep = 1.0

fig1, axs1 = plt.subplots(2, 3, tight_layout = True, figsize = (12, 8))
cmap = plt.cm.rainbow
norm = matplotlib.colors.Normalize(vmin=0, vmax=9000)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

fig2, axs2 = plt.subplots(2, 3, tight_layout = True, figsize = (12, 8))

fig3, axs3 = plt.subplots(2, 3, tight_layout = True, figsize = (12, 8))

for i in range(0, len(clusters)):
    sex_result_decam = ascii.read(data_dir + 'DECam_19_21_aug_2014_multi_exp_SEx_cat_' + clusters[i] +
                                  '_CLASS_STAR_lt_05_match_rad_1as_Gal_ext_corrected_Min_as_mag.txt')
    if i != 2:
        sex_coords_decam = SkyCoord(sex_result_decam['ALPHA_J2000'], sex_result_decam['DELTA_J2000'], unit='deg')

        sdss_galaxy_cat = ascii.read(cat_dir + 'sdss_galaxy_' + clusters[i] + '.csv')
        cat_coords = coords.SkyCoord(sdss_galaxy_cat['ra'], sdss_galaxy_cat['dec'], unit=(u.deg, u.deg))

        idx, d2d, d3d = sex_coords_decam.match_to_catalog_sky(cat_coords)
        sep_constraint = d2d.arcsec < max_sep
        sex_matches = sex_result_decam[sep_constraint]
        sdss_matches = sdss_galaxy_cat[idx[sep_constraint]]

        for j in range(0, len(bands)):
            if j == 1:
                axs1[i, j].set_title(clusters[i])
                axs2[i, j].set_title(clusters[i])
            axs1[i, j].scatter(sdss_matches['petroMag_'+bands[j]], sex_matches['MAG_ISO_'+bands[j]], alpha=0.5, s=2,
                        color=cmap(norm(sex_matches['TOTAL_EXP_TIME_'+bands[j]])))
            axs1[i, j].plot([0, 30], [0,30], '--', alpha=0.1)
            axs1[i, j].set_xlim([30, 10])
            axs1[i, j].set_ylim([30, 10])
            axs1[i, j].set_xlabel('SDSS Galaxy catalog petroMag_'+bands[j])
            axs1[i, j].set_ylabel('DECam MAG_ISO_'+bands[j]+' standardized (Gal ext corrected)')
            # axs1[i, j].gca().invert_xaxis()
            # axs1[i, j].gca().invert_yaxis()

            axs2[i, j].hist(sdss_galaxy_cat['petroMag_' + bands[j]], bins=20, label='SDSS Galaxy Catalog', density=True,
                            histtype='step')
            axs2[i, j].hist(sex_result_decam['MAG_ISO_' + bands[j]][~np.isnan(sex_result_decam['MAG_ISO_' + bands[j]])],
                            bins=20, label='DECam SEx Catalog', density=True, histtype='step')
            axs2[i, j].set_xlabel(bands[j])
            axs2[i, j].set_ylabel('N')
            if i == 0 and j == 0:
                axs2[i, j].legend()

        axs3[0, i].scatter(sdss_matches['petroMag_r'], sdss_matches['petroMag_g'] - sdss_matches['petroMag_r'],
                           alpha=0.8, s=2, label='SDSS Galaxy Catalog (matched)', color='darkgreen')
        axs3[0, i].scatter(sex_matches['MAG_ISO_r'], sex_matches['MAG_ISO_g'] - sex_matches['MAG_ISO_r'],
                           alpha=0.8, s=2, label='DECam (matched)', color='maroon')
        axs3[0, i].scatter(sdss_galaxy_cat['petroMag_r'], sdss_galaxy_cat['petroMag_g'] - sdss_galaxy_cat['petroMag_r'],
                           alpha=0.1, s=2, label='SDSS Galaxy Catalog (all)', color='lime')
        axs3[0, i].scatter(sex_result_decam['MAG_ISO_r'], sex_result_decam['MAG_ISO_g'] - sex_result_decam['MAG_ISO_r'],
                           alpha=0.1, s=2, label='DECam (CLASS_STAR < 0.5, size > 2\")', color='lightcoral')

        axs3[1, i].scatter(sdss_matches['petroMag_r'], sdss_matches['petroMag_u'] - sdss_matches['petroMag_r'],
                           alpha=0.8, s=2, label='SDSS Galaxy Catalog (matched)', color='darkgreen')
        axs3[1, i].scatter(sex_matches['MAG_ISO_r'], sex_matches['MAG_ISO_u'] - sex_matches['MAG_ISO_r'],
                           alpha=0.8, s=2, label='DECam (matched)', color='maroon')
        axs3[1, i].scatter(sdss_galaxy_cat['petroMag_r'], sdss_galaxy_cat['petroMag_u'] - sdss_galaxy_cat['petroMag_r'],
                           alpha=0.1, s=2, label='SDSS Galaxy Catalog (all)', color='lime')
        axs3[1, i].scatter(sex_result_decam['MAG_ISO_r'], sex_result_decam['MAG_ISO_u'] - sex_result_decam['MAG_ISO_r'],
                           alpha=0.1, s=2, label='DECam (CLASS_STAR < 0.5, size > 2\")', color='lightcoral')

    else:
        axs3[0, i].scatter(sex_result_decam['MAG_ISO_r'], sex_result_decam['MAG_ISO_g'] - sex_result_decam['MAG_ISO_r'],
                           alpha=0.1, s=2, label='DECam (CLASS_STAR < 0.5, size > 2\")', color='lightcoral')
        axs3[1, i].scatter(sex_result_decam['MAG_ISO_r'], sex_result_decam['MAG_ISO_u'] - sex_result_decam['MAG_ISO_r'],
                           alpha=0.1, s=2, label='DECam (CLASS_STAR < 0.5, size > 2\")', color='lightcoral')

    axs3[0, i].set_title(clusters[i])
    axs3[0, i].set_xlim([13, 23])
    axs3[1, i].set_xlim([13, 23])
    axs3[0, i].set_ylim([0, 2])
    axs3[1, i].set_ylim([0, 4])
    #axs3[0, i].set_xlabel('r\'')
    axs3[1, i].set_xlabel('r\'')
    # plt.gca().invert_xaxis()
    if i == 0:
        axs3[0, i].legend()
        axs3[0, i].set_ylabel('g\'-r\'')
        axs3[1, i].set_ylabel('u\'-r\'')

cbar = fig1.colorbar(sm)
cbar.ax.set_ylabel('Total exposure time [s]')

fig1.savefig(plot_dir + 'SDSS_cat_vs_DECam_standardized_any_size_min_mag.pdf')
fig2.savefig(plot_dir + 'SDSS_cat_vs_DECam_standardized_hist_norm_any_size_min_mag.pdf')
fig3.savefig(plot_dir + 'SDSS_cat_vs_DECam_standardized_cmd_any_size_min_mag.pdf')
