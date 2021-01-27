from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
from astropy.coordinates import SkyCoord
from astropy import coordinates as coords
import astropy.units as u
from astropy.io import fits
from astropy import wcs

sex_dir=("/Users/duhokim/work/abell/sex/cat/")
plot_dir=("/Users/duhokim/work/abell/plot/")
cat_dir=("/Users/duhokim/work/abell/cat/")
work_dir=("/Users/duhokim/work/abell/")

ver = 'v1.1'
sdss_mag_sys = 'petroMag_'
sex_mag_sys = 'MAG_PETRO'
class_star_lim = 0.9
mag_lim = [99, 99, 99]    # limit magnitude for analysis
mag_ran = [13, 30]

clusters = ['A2399', 'A2670', 'A3716']
bands = ['u', 'g', 'r']
fn = [['300_0236', '60_0239', '300_0142'],
      ['300_0409', '60_0420', '300_0351'],
      ['300_2357', '300_0030', '300_0418']]
single_dqm_fn = [['A2399_ud_300_0236_20aug.fits', 'A2399_gd_60_0239_19aug.fits', 'A2399_rd_300_0142_19aug.fits'],
             ['A2670_ud_300_0409_21aug.fits', 'A2670_gd_60_0420_21aug.fits', 'A2670_rd_300_0351_19aug.fits'],
             ['A3716_ud_300_2357_21aug.fits', 'A3716_gd_300_0030_21aug.fits', 'A3716_rd_300_0418_19aug.fits']]

exp = [[300, 60, 300],
       [300, 60, 300],
       [300, 300, 300]]
seeing = [['1.76', '1.0', '0.83'],
          ['1.28', '1.28', '0.92'],
          ['1.07', '1.07', '0.82']]
X = [[1.26, 1.26, 1.53],
           [1.3, 1.26, 1.41],
           [1.42, 1.31, 1.09]]
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

max_sep = 1.0

ccd_x = 2046
ccd_y = 4094

fig1, axs1 = plt.subplots(2, 3, tight_layout = True, figsize = (12, 8))
cmap = plt.cm.rainbow
norm = matplotlib.colors.Normalize(vmin=0, vmax=9000)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

fig2, axs2 = plt.subplots(3, 3, tight_layout = True, figsize = (12, 12))

# fig3, axs3 = plt.subplots(2, 3, tight_layout = True, figsize = (12, 8))

for i in range(0, len(clusters)):
    if i != 2:
        sdss_all = ascii.read(cat_dir + 'sdss_galaxy_' + clusters[i] + '.csv')
        # sex_all = ascii.read(sex_dir + 'DECam_19_21_aug_2014_single_best_exposure_SEx_cat_' + clusters[i] + \
        #                      '_match_rad_1as_Gal_ext_corrected_' + ver + '.txt')
        for j in range(0, len(bands)):
            cat_coords = coords.SkyCoord(sdss_all['ra'], sdss_all['dec'], unit=(u.deg, u.deg))

            sex_all = ascii.read(sex_dir + 'best_single/'+single_dqm_fn[i][j][:7]+'s_'+fn[i][j]+'.cat')
            is_gal = sex_all['CLASS_STAR'] < class_star_lim  # only galaxies
            sex_gal_cat = sex_all[is_gal]

            # read Data Quality Mask (DQM) fits file to check saturation
            hdu_r = fits.open(work_dir + 'fits/best_single/' + single_dqm_fn[i][2])

            # make a boolean array to set FLAGS to SATURATED = 4 (could not directly)
            is_dqm_r_zero = np.ones(len(sex_gal_cat), dtype=bool)

            for k in range(0, len(sex_gal_cat)):
                # read the celestial coordinate
                cel_coord = [[sex_gal_cat['ALPHA_J2000'][k], sex_gal_cat['DELTA_J2000'][k]], [0, 0]]
                sky_coord = SkyCoord(sex_gal_cat['ALPHA_J2000'][i], sex_gal_cat['DELTA_J2000'][i], unit='deg')
                # for every CCD
                for k in range(1, 61):
                    # read WCS
                    w = wcs.WCS(hdu_r[k].header)
                    pixcrd = w.wcs_world2pix(cel_coord, 1)
                    # if w.footprint_contains(sky_coord):
                    if (pixcrd[0][0] > 0) & (pixcrd[0][0] < ccd_x) & (pixcrd[0][1] > 0) & (pixcrd[0][1] < ccd_y):
                        dqm_fits = hdu_r[k].data
                        # check the value of DQM
                        if dqm_fits[int(pixcrd[0][1])][int(pixcrd[0][0])] % 128:
                            # toggle on saturated bool array
                            is_dqm_r_zero[i] = False

            sex_coords = SkyCoord(sex_gal_cat['ALPHA_J2000'][is_dqm_r_zero], sex_gal_cat['DELTA_J2000'][is_dqm_r_zero], unit='deg')

            idx, d2d, d3d = sex_coords.match_to_catalog_sky(cat_coords)
            sep_constraint = d2d.arcsec < max_sep
            sex_matches = sex_gal_cat[is_dqm_r_zero][sep_constraint]
            sdss_matches = sdss_all[idx[sep_constraint]]

            # with open(cat_dir+clusters[i]+'_'+bands[j]+'_petroMag_MAG_PETRO.cat', 'w') as f:
            #     for k in range(0, len(sex_matches)):
            #         f.write(str(sex_matches['ALPHA_J2000'][k])+'\t'+str(sex_matches['DELTA_J2000'][k])+'\t'+
            #                 str(sex_matches['MAG_PETRO'][k]+a[i][j]+b[i][j]*X[i][j]+2.5*np.log10(exp[i][j]))+'\t'+str(sex_matches['MAGERR_PETRO'][k])+'\t'+
            #                 str(sdss_matches['petroMag_'+bands[j]][k])+'\t'+str(sdss_matches['petroMagErr_'+bands[j]][k])+'\n')

            if j == 1:
                axs1[i, j].set_title(clusters[i])
                axs2[i, j].set_title(clusters[i])
            axs1[i, j].scatter(sdss_matches[sdss_mag_sys+bands[j]],
                               sex_matches[sex_mag_sys]+a[i][j]+b[i][j]*X[i][j]+2.5*np.log10(exp[i][j])-gal_ext[i][j],
                               alpha=0.2, s=2)
            axs1[i, j].plot(mag_ran, mag_ran, '--', alpha=0.5)
            # axs1[i, j].errorbar(sdss_matches[sdss_mag_sys+bands[j]],
            #                     sex_matches['MAG_ISO_'+bands[j]],
            #                     yerr = sex_matches['MAGERR_ISO_'+bands[j]],
            #                     xerr = sdss_matches['cModelMagErr_'+bands[j]],
            #                     fmt='o',
            #                     ms=1,
            #                     alpha=0.1)
            axs1[i, j].set_xlim(mag_ran.reverse())
            axs1[i, j].set_ylim(mag_ran.reverse())
            axs1[i, j].set_xlabel('SDSS Galaxy '+sdss_mag_sys+bands[j])
            axs1[i, j].set_ylabel('DECam '+sex_mag_sys+'_'+bands[j]+' standardized')
            axs1[i, j].text(mag_ran[1]-1, mag_ran[0]+1, single_dqm_fn[i][j][-2:]+'Aug')
            axs1[i, j].text(mag_ran[1]-1, mag_ran[0]+2, str(exp[i][j])+'s')
            axs1[i, j].text(mag_ran[1]-1, mag_ran[0]+3, 'seeing: '+seeing[i][j])
            axs1[i, j].text(mag_ran[1]-1, mag_ran[0]+4, 'X: ' + str(X[i][j]))
            axs1[i, j].text(mag_ran[1]-1, mag_ran[0]+5, '#: ' + str(len(sdss_matches)))
            axs1[i, j].set_xlim(mag_ran)
            axs1[i, j].set_ylim(mag_ran)
            axs1[i, j].invert_xaxis()
            axs1[i, j].invert_yaxis()


            axs2[i, j].hist(sdss_all[sdss_mag_sys + bands[j]], bins = np.arange(mag_ran[0], mag_ran[1], 0.5),
                            label=r'SDSS Galaxy Catalog (~4deg$^{2}$)',
                            color='teal',
                            histtype='step')
            axs2[i, j].hist(sex_gal_cat[sex_mag_sys]+a[i][j]+b[i][j]*X[i][j]+2.5*np.log10(exp[i][j]),
                            bins = np.arange(mag_ran[0], mag_ran[1], 0.5),
                            label=r'DECam SEx Catalog (~3deg$^{2}$)',
                            color='orangered',
                            histtype='step')
            axs2[i, j].set_xlabel(bands[j])
            axs2[i, j].set_ylabel('N')
            axs2[i, j].set_xlim(mag_ran)
            if j == 0:
                axs1[i, j].axhline(y=21, linestyle='--', color='k', alpha=0.5)
                axs2[i, j].axvline(x=21, linestyle='--', color='k', alpha=0.5)
                # axs2[i, j].plot([21, 21], [0, plt.axvline], '--', alpha=0.5, clip_on=True)
            else:
                axs1[i, j].axhline(y=22, linestyle='--', color='k', alpha=0.5)
                axs2[i, j].axvline(x=22, linestyle='--', color='k', alpha=0.5)
                # axs2[i, j].plot([22, 22], [0, plt.axvline], '--', alpha=0.5, clip_on=True)
            if i == 0 and j == 0:
                axs2[i, j].legend(fontsize='small')
    else:
        for j in range(0, len(bands)):
            sex_all = ascii.read(sex_dir + 'best_single/' + single_dqm_fn[i][j][:7] + 's_' + fn[i][j] + '.cat')
            is_gal = sex_all['CLASS_STAR'] < class_star_lim  # only galaxies
            sex_gal_cat = sex_all[is_gal]

            # read Data Quality Mask (DQM) fits file to check saturation
            hdu_r = fits.open(work_dir + 'fits/best_single/' + single_dqm_fn[i][2])

            # make a boolean array to set FLAGS to SATURATED = 4 (could not directly)
            is_dqm_r_zero = np.ones(len(sex_gal_cat), dtype=bool)

            for k in range(0, len(sex_gal_cat)):
                # read the celestial coordinate
                cel_coord = [[sex_gal_cat['ALPHA_J2000'][k], sex_gal_cat['DELTA_J2000'][k]], [0, 0]]
                sky_coord = SkyCoord(sex_gal_cat['ALPHA_J2000'][i], sex_gal_cat['DELTA_J2000'][i], unit='deg')
                # for every CCD
                for l in range(1, 61):
                    # read WCS
                    w = wcs.WCS(hdu_r[l].header)
                    pixcrd = w.wcs_world2pix(cel_coord, 1)
                    # if w.footprint_contains(sky_coord):
                    if (pixcrd[0][0] > 0) & (pixcrd[0][0] < ccd_x) & (pixcrd[0][1] > 0) & (pixcrd[0][1] < ccd_y):
                        dqm_fits = hdu_r[l].data
                        # check the value of DQM
                        if dqm_fits[int(pixcrd[0][1])][int(pixcrd[0][0])] % 128:
                            # toggle on saturated bool array
                            is_dqm_r_zero[i] = False

            axs2[i, j].hist(sex_gal_cat[sex_mag_sys][is_dqm_r_zero] + a[i][j] + b[i][j] * X[i][j] + 2.5 * np.log10(exp[i][j]),
                            bins=np.arange(mag_ran[0], mag_ran[1], 0.5),
                            label=r'DECam SEx Catalog (~3deg$^{2}$)',
                            color='orangered',
                            histtype='step')
            axs2[i, j].set_xlabel(bands[j])
            axs2[i, j].set_ylabel('N')
            axs2[i, j].set_xlim(mag_ran)
            if j == 0:
                axs2[i, j].axvline(x=21, linestyle='--', color='k', alpha=0.5)
            else:
                axs2[i, j].axvline(x=22, linestyle='--', color='k', alpha=0.5)
            if j == 1:
                axs2[i, j].set_title(clusters[i])


fig1.savefig(plot_dir + 'SDSS_cat_vs_DECam_standardized_one_to_one_best_single_aper_corr_no_Gal_ext_'+ver+'.png')
fig2.savefig(plot_dir + 'SDSS_cat_vs_DECam_standardized_hist_best_single_no_Gal_ext_'+ver+'.png')

