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


clusters = ['A2399', 'A2670']
bands = ['u', 'g', 'r']
fn = [['300_0236_20', '60_0239_19', '300_0142_19'],
      ['300_0409_21', '60_0420_21', '300_0351_19'],
      ['300_0008_21', '300_0030_21', '300_0418_19']]
exp = [[300, 60, 300],
       [300, 60, 300],
       [300, 300, 300]]
seeing = [['1.76', '1.0', '0.83'],
          ['1.28', '1.28', '0.92'],
          ['1.07', '1.07', '0.82']]
X = [[1.26, 1.26, 1.53],
           [1.3, 1.26, 1.41],
           [1.38, 1.31, 1.09]]
a = [ [-1.257, 0.381, 0.567],
      [-1.796, 0.179, 0.567],
      [-1.796, 0.179, 0.567]]
b = [ [-0.592, -0.186, -0.136],
      [-0.132, -0.036, -0.136],
      [-0.132, -0.036, -0.136]]
# gal_ext = [ [0.159, 0.124, 0.086], # irsa.ipac.caltech.edu, S and F (2011)
#             [0.188, 0.146, 0.101],
#             [0.157, 0.122, 0.085]]
gal_ext = [ [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0]]
mag_lim = [20, 21, 21]    # limit magnitude for analysis
star_lim = [0.5, 0.5, 0.5]
max_sep = 1.0

fig1, axs1 = plt.subplots(2, 3, tight_layout = True, figsize = (12, 8))
cmap = plt.cm.rainbow
norm = matplotlib.colors.Normalize(vmin=0, vmax=9000)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

fig2, axs2 = plt.subplots(2, 3, tight_layout = True, figsize = (12, 8))

# fig3, axs3 = plt.subplots(2, 3, tight_layout = True, figsize = (12, 8))

for i in range(0, len(clusters)):
    if i != 2:
        for j in range(0, len(bands)):
            # sdss_galaxy_all = ascii.read(cat_dir + 'sdss_galaxy_' + clusters[i] + '_cModel.csv')
            # sdss_galaxy_flag = sdss_galaxy_all['cModelMag_'+bands[j]] < mag_lim[j]
            sdss_galaxy_all = ascii.read(cat_dir + 'sdss_galaxy_' + clusters[i] + '.csv')
            sdss_galaxy_flag = sdss_galaxy_all['petroMag_'+bands[j]] < mag_lim[j]
            sdss_galaxy_cat = sdss_galaxy_all[sdss_galaxy_flag]
            cat_coords = coords.SkyCoord(sdss_galaxy_cat['ra'], sdss_galaxy_cat['dec'], unit=(u.deg, u.deg))

            sex_all = ascii.read(data_dir + clusters[i] + '/' + clusters[i]+ '_' + bands[j] + 'a_' + fn[i][j]
                                          + 'aug.cat')
            sex_flag = (sex_all['FLAGS'] == 0) & (sex_all['MAG_PETRO']+a[i][j]+b[i][j]*X[i][j]+2.5*np.log10(exp[i][j])
                       - gal_ext[i][j] < mag_lim[j]) & (sex_all['CLASS_STAR'] < star_lim[j])
            sex_result_decam = sex_all[sex_flag]
            sex_coords_decam = SkyCoord(sex_result_decam['ALPHA_J2000'], sex_result_decam['DELTA_J2000'], unit='deg')

            idx, d2d, d3d = sex_coords_decam.match_to_catalog_sky(cat_coords)
            sep_constraint = d2d.arcsec < max_sep
            sex_matches = sex_result_decam[sep_constraint]
            sdss_matches = sdss_galaxy_cat[idx[sep_constraint]]

            with open(cat_dir+clusters[i]+'_'+bands[j]+'_petroMag_MAG_PETRO.cat', 'w') as f:
                for k in range(0, len(sex_matches)):
                    f.write(str(sex_matches['ALPHA_J2000'][k])+'\t'+str(sex_matches['DELTA_J2000'][k])+'\t'+
                            str(sex_matches['MAG_PETRO'][k]+a[i][j]+b[i][j]*X[i][j]+2.5*np.log10(exp[i][j]))+'\t'+str(sex_matches['MAGERR_PETRO'][k])+'\t'+
                            str(sdss_matches['petroMag_'+bands[j]][k])+'\t'+str(sdss_matches['petroMagErr_'+bands[j]][k])+'\n')

            if j == 1:
                axs1[i, j].set_title(clusters[i])
                axs2[i, j].set_title(clusters[i])
    #        axs1[i, j].scatter(sdss_matches['petroMag_'+bands[j]],
    #                            sex_matches['MAG_PETRO']+a[i][j]+b[i][j]*X[i][j]+2.5*np.log10(exp[i][j])-gal_ext[i][j],
    #                            alpha=0.2, s=2)
            axs1[i, j].plot([0, 30], [0,30], '--', alpha=0.5)
            axs1[i, j].errorbar(sdss_matches['petroMag_'+bands[j]],
                                sex_matches['MAG_PETRO']+a[i][j]+b[i][j]*X[i][j]+2.5*np.log10(exp[i][j])-gal_ext[i][j],
                                yerr = sex_matches['MAGERR_PETRO'],
                                xerr = sdss_matches['petroMagErr_'+bands[j]],
                                fmt='o',
                                ms=1,
                                alpha=0.1)
            axs1[i, j].set_xlim([22, 10])
            axs1[i, j].set_ylim([22, 10])
            axs1[i, j].set_xlabel('SDSS Galaxy petroMag_'+bands[j])
            axs1[i, j].set_ylabel('DECam MAG_PETRO_'+bands[j]+' standardized (Gal ext corrected)')
            axs1[i, j].text(21, 11, fn[i][j][-2:]+'Aug')
            axs1[i, j].text(21, 12, str(exp[i][j])+'s')
            axs1[i, j].text(21, 13, 'seeing: '+seeing[i][j])
            axs1[i, j].text(21, 14, 'X: ' + str(X[i][j]))
            axs1[i, j].text(21, 15, 'N: ' + str(len(sdss_matches)))
            # axs1[i, j].invert_xaxis()
            # axs1[i, j].invert_yaxis()

            axs2[i, j].hist(sdss_galaxy_cat['petroMag_' + bands[j]], bins=20, label='SDSS Galaxy Catalog', density=True,
                            histtype='step')
            axs2[i, j].hist(sex_result_decam['MAG_PETRO']+a[i][j]+b[i][j]*X[i][j]+2.5*np.log10(exp[i][j])-gal_ext[i][j],
                            bins=20, label='DECam SEx Catalog', density=True, histtype='step')
            axs2[i, j].set_xlabel(bands[j])
            axs2[i, j].set_ylabel('N')
            if i == 0 and j == 0:
                axs2[i, j].legend()

        # axs3[0, i].scatter(sdss_matches['cModelMag_r'], sdss_matches['cModelMag_g'] - sdss_matches['cModelMag_r'],
        #                    alpha=0.8, s=2, label='SDSS Galaxy Catalog (matched)', color='darkgreen')
        # axs3[0, i].scatter(sex_matches['MAG_ISO'], sex_matches['MAG_ISO_g'] - sex_matches['MAG_ISO_r'],
        #                    alpha=0.8, s=2, label='DECam (matched)', color='maroon')
        # axs3[0, i].scatter(sdss_galaxy_cat['petroMag_r'], sdss_galaxy_cat['petroMag_g'] - sdss_galaxy_cat['petroMag_r'],
        #                    alpha=0.1, s=2, label='SDSS Galaxy Catalog (all)', color='lime')
        # axs3[0, i].scatter(sex_result_decam['MAG_ISO_r'], sex_result_decam['MAG_ISO_g'] - sex_result_decam['MAG_ISO_r'],
        #                    alpha=0.1, s=2, label='DECam (CLASS_STAR < 0.5, size > 2\")', color='lightcoral')
        #
        # axs3[1, i].scatter(sdss_matches['petroMag_r'], sdss_matches['petroMag_u'] - sdss_matches['petroMag_r'],
        #                    alpha=0.8, s=2, label='SDSS Galaxy Catalog (matched)', color='darkgreen')
        # axs3[1, i].scatter(sex_matches['MAG_ISO_r'], sex_matches['MAG_ISO_u'] - sex_matches['MAG_ISO_r'],
        #                    alpha=0.8, s=2, label='DECam (matched)', color='maroon')
        # axs3[1, i].scatter(sdss_galaxy_cat['petroMag_r'], sdss_galaxy_cat['petroMag_u'] - sdss_galaxy_cat['petroMag_r'],
        #                    alpha=0.1, s=2, label='SDSS Galaxy Catalog (all)', color='lime')
        # axs3[1, i].scatter(sex_result_decam['MAG_ISO_r'], sex_result_decam['MAG_ISO_u'] - sex_result_decam['MAG_ISO_r'],
        #                    alpha=0.1, s=2, label='DECam (CLASS_STAR < 0.5, size > 2\")', color='lightcoral')

    # else:
    #     axs3[0, i].scatter(sex_result_decam['MAG_ISO_r'], sex_result_decam['MAG_ISO_g'] - sex_result_decam['MAG_ISO_r'],
    #                        alpha=0.1, s=2, label='DECam (CLASS_STAR < 0.5, size > 2\")', color='lightcoral')
    #     axs3[1, i].scatter(sex_result_decam['MAG_ISO_r'], sex_result_decam['MAG_ISO_u'] - sex_result_decam['MAG_ISO_r'],
    #                        alpha=0.1, s=2, label='DECam (CLASS_STAR < 0.5, size > 2\")', color='lightcoral')

    # axs3[0, i].set_title(clusters[i])
    # axs3[0, i].set_xlim([13, 23])
    # axs3[1, i].set_xlim([13, 23])
    # axs3[0, i].set_ylim([0, 2])
    # axs3[1, i].set_ylim([0, 4])
    # #axs3[0, i].set_xlabel('r\'')
    # axs3[1, i].set_xlabel('r\'')
    # # plt.gca().invert_xaxis()
    # if i == 0:
    #     axs3[0, i].legend()
    #     axs3[0, i].set_ylabel('g\'-r\'')
    #     axs3[1, i].set_ylabel('u\'-r\'')

# cbar = fig1.colorbar(sm)
# cbar.ax.set_ylabel('Total exposure time [s]')

fig1.savefig(plot_dir + 'SDSS_cat_vs_DECam_standardized_best_single_exp_aper_corr_petro_only_bright_no_Gal_ext.png')
fig2.savefig(plot_dir + 'SDSS_cat_vs_DECam_standardized_hist_norm_best_single_aper_corr_petro_only_bright_no_Gal_ext.png')
# fig3.savefig(plot_dir + 'SDSS_cat_vs_DECam_standardized_cmd_any_size_min_mag.pdf')
