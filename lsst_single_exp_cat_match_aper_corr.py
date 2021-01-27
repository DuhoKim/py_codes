from astropy.io import ascii
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import matplotlib.colors
from astropy.coordinates import SkyCoord
from astropy import coordinates as coords
import astropy.units as u
from astropy.io import fits
from astropy import wcs
from astropy.table import Table, vstack

sex_dir=("/Users/duhokim/work/abell/sex/cat/")
plot_dir=("/Users/duhokim/work/abell/plot/")
cat_dir=("/Users/duhokim/work/abell/cat/")
work_dir=("/Users/duhokim/work/abell/")
# lsst_dir=("/Users/duhokim/lsst_stack/demo_data_newinstall_decam/DATA/rerun/coaddForcedPhot/deepCoadd-results/")
# sav_dir=("/Users/duhokim/lsst_stack/demo_data_newinstall_decam/")
lsst_dir=("/Users/duhokim/lsst_stack/demo_data_orig_20_intersect/DATA/rerun/coaddForcedPhot/deepCoadd-results/")
sav_dir=("/Users/duhokim/lsst_stack/demo_data_orig_20_intersect/")

visit_nums = ['0350110', '0350835']
'''
patches = [
    ['0,2', '0,3', '0,4', '1,3', '1,4', '1,5',        '2,2',               '3,1', '3,3',        '4,1', '4,2', '4,3',
     '4,5', '5,1', '5,2', '5,3', '5,5', '6,2', '6,3'],
    ['0,2', '0,3', '0,4', '1,3', '1,4', '1,5', '2,1', '2,2', '2,6', '3,0', '3,1', '3,3', '4,0', '4,1', '4,2', '4,3',
     '4,5',        '5,2', '5,3', '5,5', '6,2', '6,3']
    ]

patches_all = ['0,2', '0,3', '0,4', '1,3', '1,4', '1,5', '2,2', '3,1', '3,3', '4,1', '4,2', '4,3', '4,5', '5,2', '5,3',
               '5,5', '6,2', '6,3']
'''
patches = [
    ['0,2', '0,3',        '1,2',       '2,2', '2,3',              '3,2',               '4,2',        '5,2', '6,2'],
    ['0,2', '0,3', '1,1', '1,2', '2,0',             '3,0', '3,1', '3,2', '4,0', '4,1', '4,2', '5,1', '5,2']
    ]

patches_all = ['0,2', '0,3', '1,2', '3,2', '4,2', '5,2']

ver = 'v1.1'
sdss_mag_sys = 'petroMag_'
sex_mag_sys = 'MAG_PETRO'
class_star_lim = 0.9
mag_lim = [99, 99, 99]    # limit magnitude for analysis
mag_ran = [25, 13]

clusters = ['A2670']
bands = ['g', 'r']
fn =  ['60_0420', '300_0351']
single_dqm_fn = ['A2670_gd_60_0420_21aug.fits', 'A2670_rd_300_0351_19aug.fits']

exp = [60, 300]
seeing = ['1.28', '0.92']
X = [1.26, 1.41]
a = [0.179, 0.567]
b = [-0.036, -0.136]
# gal_ext = [ [0.159, 0.124, 0.086], # irsa.ipac.caltech.edu, S and F (2011)
#             [0.188, 0.146, 0.101],
#             [0.157, 0.122, 0.085]]
gal_ext = [0, 0]

max_sep = 1.0

ccd_x = 2046
ccd_y = 4094

fig, axs = plt.subplots(2, 2, tight_layout = True, figsize = (10, 10))

sex_cat = ascii.read(sex_dir + 'DECam_19_21_aug_2014_single_best_exposure_SEx_cat_A2670_match_rad_1as_'
                                                                  'Gal_ext_corrected_' + ver + '.txt')
sex_coords = SkyCoord(sex_cat['ALPHA_J2000'], sex_cat['DELTA_J2000'], unit='deg')

for j in range(0, len(bands)):
    for i in range(0, len(patches[j])):
        rMags = np.load(sav_dir+bands[j]+'mag'+patches[j][i][0]+patches[j][i][2]+'.npy')
        lsst_all = Table.read(lsst_dir + bands[j] + '/0/' + patches[j][i] + '/forcedSrc-'+bands[j]+'-0-'+
                              patches[j][i]+'.fits')

        is_not_nan = [(not np.isnan(rMags[x][0])) & (not np.isnan(rMags[x][1])) & (rMags[x][1] < 0.1) for x in range(0, len(rMags))]
        rMags = rMags[is_not_nan]
        lsst_all = lsst_all[is_not_nan]

        lsst_coords = SkyCoord(lsst_all['coord_ra'], lsst_all['coord_dec'], unit=(u.rad, u.rad))

        idx, d2d, d3d = sex_coords.match_to_catalog_sky(lsst_coords)
        sep_constraint = d2d.arcsec < max_sep
        sex_matches = sex_cat[sep_constraint]
        lsst_matches = rMags[idx[sep_constraint]]

        axs[0, j].scatter(sex_matches['MAG_AUTO_'+bands[j]],
                          lsst_matches[:, 0],
                          alpha=0.2, s=1)

        axs[0, j].scatter(sex_matches['MAG_AUTO_'+bands[j]],
                          lsst_matches[:, 0],
                          alpha=0.2, s=1)

for i in range(0, len(patches_all)):
    gMags = np.load(sav_dir+'gmag'+patches_all[i][0]+patches_all[i][2]+'.npy')
    rMags = np.load(sav_dir+'rmag'+patches_all[i][0]+patches_all[i][2]+'.npy')

    is_not_nan = [(not np.isnan(gMags[x][0])) & (not np.isnan(gMags[x][1])) & (not np.isnan(rMags[x][0])) &
                  (not np.isnan(rMags[x][1])) & (gMags[x][1] < 0.1) & (rMags[x][1] < 0.1) for x in range(0, len(rMags))]

    if i == 0:
        lsst_g = gMags[is_not_nan, 0]
        lsst_r = rMags[is_not_nan, 0]
    else:
        lsst_g = np.append(lsst_g, gMags[is_not_nan, 0], axis=0)
        lsst_r = np.append(lsst_r, rMags[is_not_nan, 0], axis=0)

axs[1, 0].scatter(sex_cat['MAG_AUTO_'+bands[1]],
                  sex_cat['MAG_AUTO_'+bands[0]] - sex_cat['MAG_AUTO_'+bands[1]],
                  alpha=0.2, s=1, label='SEx (CLASS_STAR < 0.9 & DQM_cen = 0,128)')

axs[1, 1].scatter(lsst_r,
                  lsst_g - lsst_r,
                  alpha=0.2, s=1, label='LSST (merr < 0.1)')

axs[0, 0].plot(mag_ran, mag_ran, linestyle=':', alpha=0.5)
axs[0, 1].plot(mag_ran, mag_ran, linestyle=':', alpha=0.5)

axs[0, 0].text(27.5, 13, 'A2670', fontsize=24)
axs[0, 0].text(27.5, 13, 'A2670', fontsize=24)
axs[0, 0].set_title('g band (60s exposure)')
axs[0, 1].set_title('r band (300s exposure)')

axs[0, 0].set_xlabel('SEx (MAG_AUTO)')
axs[0, 1].set_xlabel('SEx (MAG_AUTO)')
axs[0, 0].set_ylabel('LSST (base_PsfFlux) (merr < 0.1)')
axs[0, 1].set_ylabel('LSST (base_PsfFlux) (merr < 0.1)')

axs[0, 0].set_xlim(mag_ran)
axs[0, 1].set_xlim(mag_ran)
axs[0, 0].set_ylim(mag_ran)
axs[0, 1].set_ylim(mag_ran)

axs[1, 0].set_xlabel('r')
axs[1, 0].set_ylabel('g-r')
axs[1, 0].set_xlim(mag_ran)
axs[1, 0].set_ylim([-3, 5])
axs[1, 0].legend()

axs[1, 1].set_xlabel('r')
axs[1, 1].set_ylabel('g-r')
axs[1, 1].set_xlim(mag_ran)
axs[1, 1].set_ylim([-3, 5])
axs[1, 1].legend()

fig.savefig(plot_dir + 'lsst_rMag_sav_vs_sex_single_exp_merr_lt_01_orig_20_intersect.png')


