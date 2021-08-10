from astropy.io import ascii
import numpy as np
import time
from astropy.coordinates import SkyCoord
import my_module
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from datetime import datetime
from astropy import units as u
from reproject.mosaicking import find_optimal_celestial_wcs, reproject_and_coadd
from reproject import reproject_interp
from image_registration import chi2_shift
from image_registration.fft_tools import shift
from astropy.visualization import make_lupton_rgb
import matplotlib.pyplot as plt

work_dir=("/Users/duhokim/work/abell/")

# clusters = ['A2399', 'A2670', 'A3716']
clusters = ['A2670']
fn = [['A2670_usi.fits', 'A2670_gsi.fits', 'A2670_rsi.fits']]
coords_cl_cen = SkyCoord(358.557083, -10.418889, unit='deg')
thumb = 5000   # half size of rgb mosaic

# stacked zero point (ZP; mag = mag0 (ZP_inst=25) + ZP)
stack_a = [
    [4.05, 6.06, 6.38],
    [3.96, 6.04, 6.22],
    [3.96, 6.15, 6.42]
]

for k in range(0, len(clusters)):
    with fits.open(work_dir + 'fits/stacked/' + fn[k][0]) as hdu_u, \
            fits.open(work_dir + 'fits/stacked/' + fn[k][1]) as hdu_g, \
            fits.open(work_dir + 'fits/stacked/' + fn[k][2]) as hdu_r:

        # ax1 = plt.subplot(1, 2, 1, projection=WCS(hdu_r[5].header))
        # ax1.imshow(hdu_r[5].data, origin='lower', vmin=1000., vmax=2000.)
        # ax1.coords['ra'].set_axislabel('Right Ascension')
        # ax1.coords['dec'].set_axislabel('Declination')
        # ax1.set_title('r')
        #
        # ax2 = plt.subplot(1, 2, 2, projection=WCS(hdu_u[1].header))
        # ax2.imshow(hdu_u[1].data, origin='lower', vmin=20, vmax=100)
        # # print(f'{np.min(hdu_u[5].data)} {np.max(hdu_u[5].data)} {np.nanmedian(hdu_u[5].data)}')
        # ax2.coords['ra'].set_axislabel('Right Ascension')
        # ax2.coords['dec'].set_axislabel('Declination')
        # ax2.coords['dec'].set_axislabel_position('r')
        # ax2.coords['dec'].set_ticklabel_position('r')
        # ax2.set_title('u')


        array, footprint = reproject_interp(hdu_u[5], hdu_r[5].header)
        #
        ax1 = plt.subplot(1,2,1, projection = WCS(hdu_r[5].header))
        ax1.imshow(array, origin='lower', vmin = 1000, vmax=2000)
        ax1.coords['ra'].set_axislabel('RA')
        ax1.coords['dec'].set_axislabel('Dec')
        ax1.set_title('Reprojected u')

        ax2 = plt.subplot(1,2,2, projection = WCS(hdu_r[5].header))
        ax2.imshow(footprint, origin='lower', vmin = 0, vmax=1.5)
        ax2.coords['ra'].set_axislabel('RA')
        ax2.coords['dec'].set_axislabel('Dec')
        ax2.coords['dec'].set_axislabel_position('r')
        ax2.coords['dec'].set_ticklabel_position('r')
        ax2.set_title('footprint')

        # fig.savefig(work_dir + 'fits/stacked/rgb_shape_out_new.jpg')

        # plt.imshow(array_rgb)

        # plt.figure(figsize=(12, 12))
        # wcs_out_u, shape_out_u = find_optimal_celestial_wcs(hdu_u[1:],
        #                                                     reference=coords_cl_cen)
        # wcs_out_g, shape_out_g = find_optimal_celestial_wcs(hdu_g[1:],
        #                                                     reference=coords_cl_cen)
        # wcs_out_r, shape_out_r = find_optimal_celestial_wcs(hdu_r[1:],
        #                                                     reference=coords_cl_cen)
        # array_u, footprint_u = reproject_and_coadd(hdu_u[1:],
        #                                            wcs_out_u,
        #                                            shape_out=shape_out_r,
        #                                            reproject_function=reproject_interp)
        # array_g, footprint_g = reproject_and_coadd(hdu_g[1:],
        #                                            wcs_out_g,
        #                                            shape_out=shape_out_r,
        #                                            reproject_function=reproject_interp)
        # array_r, footprint_r = reproject_and_coadd(hdu_r[1:],
        #                                            wcs_out_r,
        #                                            shape_out=shape_out_r,
        #                                            reproject_function=reproject_interp)

        # xoff_u, yoff_u, exoff_u, eyoff_u = chi2_shift(array_r[:, :array_u.shape[1]],
        #                                               array_u[:array_r.shape[0], :],
        #                                               return_error=True,
        #                                               upsample_factor='auto')
        # xoff_g, yoff_g, exoff_g, eyoff_g = chi2_shift(array_r[:array_g.shape[0], :array_g.shape[1]],
        #                                               array_g,
        #                                               return_error=True,
        #                                               upsample_factor='auto')
        #
        # mosaic_r = array_r[int(array_r.shape[0] / 2) - thumb: int(array_r.shape[0] / 2) + thumb,
        #            int(array_r.shape[1] / 2) - thumb: int(array_r.shape[1] / 2) + thumb]
        # mosaic_g = array_g[int(array_g.shape[0] / 2) - thumb - xoff_g: int(array_g.shape[0] / 2) + thumb - xoff_g,
        #            int(array_g.shape[1] / 2) - thumb - yoff_g: int(array_g.shape[1] / 2) + thumb - yoff_g]
        # mosaic_b = array_u[int(array_u.shape[0] / 2) - thumb - xoff_u: int(array_u.shape[0] / 2) + thumb - xoff_u,
        #            int(array_u.shape[1] / 2) - thumb - yoff_u: int(array_u.shape[1] / 2) + thumb - yoff_u]

        # mosaic_r = array_r[int(array_r.shape[0] / 2) - thumb: int(array_r.shape[0] / 2) + thumb,
        #            int(array_r.shape[1] / 2) - thumb: int(array_r.shape[1] / 2) + thumb]
        # mosaic_g = array_g[int(array_g.shape[0] / 2) - thumb: int(array_g.shape[0] / 2) + thumb,
        #            int(array_g.shape[1] / 2) - thumb: int(array_g.shape[1] / 2) + thumb]
        # mosaic_b = array_u[int(array_u.shape[0] / 2) - thumb: int(array_u.shape[0] / 2) + thumb,
        #            int(array_u.shape[1] / 2) - thumb: int(array_u.shape[1] / 2) + thumb]

        # r_ran = [np.nanmedian(mosaic_r), np.nanmedian(mosaic_r) + np.nanstd(mosaic_r)]
        # g_ran = [np.nanmedian(mosaic_g), np.nanmedian(mosaic_g) + np.nanstd(mosaic_g)]
        # b_ran = [np.nanmedian(mosaic_b), np.nanmedian(mosaic_b) + np.nanstd(mosaic_b)]
        #
        # mosaic_r[mosaic_r < r_ran[0]] = 0
        # mosaic_g[mosaic_g < g_ran[0]] = 0
        # mosaic_b[mosaic_b < b_ran[0]] = 0
        #
        # ind_r = np.where((mosaic_r > r_ran[0]) & (mosaic_r < r_ran[1]))
        # ind_g = np.where((mosaic_g > g_ran[0]) & (mosaic_g < g_ran[1]))
        # ind_b = np.where((mosaic_b > b_ran[0]) & (mosaic_b < b_ran[1]))
        #
        # mosaic_r[ind_r] = np.log(mosaic_r[ind_r] - r_ran[0]) / np.log(r_ran[1]-mosaic_r[ind_r])
        # mosaic_g[ind_g] = np.log(mosaic_g[ind_g] - g_ran[0]) / np.log(g_ran[1]-mosaic_g[ind_g])
        # mosaic_b[ind_b] = np.log(mosaic_b[ind_b] - b_ran[0]) / np.log(b_ran[1]-mosaic_b[ind_b])
        #
        # mosaic_r[mosaic_r > r_ran[1]] = 1
        # mosaic_g[mosaic_g > g_ran[1]] = 1
        # mosaic_b[mosaic_b > b_ran[1]] = 1

        # array_rgb = make_lupton_rgb(mosaic_r,
        #                             mosaic_g,
        #                             mosaic_b,
        #                             Q=10,
        #                             stretch=0.5,
        #                             minimum=[np.nanmedian(array_r), np.nanmedian(array_g), np.nanmedian(array_u)])
        # corrected_image_u = shift.shiftnd(offset_image, -yoff, -xoff)
        #

        # fig = plt.gca(projection=wcs_out_r)
        # plt.imshow(array_rgb)
        # fig.savefig(work_dir + 'fits/stacked/rgb_shape_out_new.jpg')


