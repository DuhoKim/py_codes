import numpy as np
import my_module as mm
from astropy.io import fits
import matplotlib.pyplot as plt
from reproject.mosaicking import find_optimal_celestial_wcs, reproject_and_coadd
from reproject import reproject_interp
from astropy.visualization.wcsaxes import SphericalCircle
from astropy import units as u
import abell_cluster_module as ab
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.visualization import make_lupton_rgb
from astropy.io import ascii
from astropy import wcs
import my_module as mm
import sys
from astropy.convolution import Gaussian2DKernel
from scipy.signal import convolve as scipy_convolve
from astropy.convolution import convolve


ns = range(1, 15)
n = 7



with fits.open(ab.work_dir+'pics/A754_cutout/cutout_u.fits') as hdu_u,  \
    fits.open(ab.work_dir + 'pics/A754_cutout/cutout_g.fits') as hdu_g,   \
    fits.open(ab.work_dir + 'pics/A754_cutout/cutout_r.fits') as hdu_r:

    cutout_u = hdu_u[0].data
    cutout_g = hdu_g[0].data
    cutout_r = hdu_r[0].data

    calib_u = cutout_u * 10 ** (-ab.stack_a[0][0]/2.5)
    calib_g = cutout_g * 10 ** (-ab.stack_a[0][1]/2.5)
    calib_r = cutout_r * 10 ** (-ab.stack_a[0][2]/2.5)

    # conv_u = scipy_convolve(cutout_u, kernel, mode='same', method='direct')
    # conv_g = scipy_convolve(cutout_g, kernel, mode='same', method='direct')
    # conv_r = scipy_convolve(cutout_r, kernel, mode='same', method='direct')

    plt.figure(figsize=(5, 5))

    # qs = [0.2 * exp for exp in range(1, 6)]
    # stretches = [0.007 * exp for exp in range(1, 6)]

    # qs = [8 * exp for exp in range(1, 2)]
    # stretches = [0.08 * exp for exp in range(1, 2)]
    # x_stddev = [1.0 * exp for exp in range(1, 2)]
    #
    # for x in x_stddev:
    #     kernel = Gaussian2DKernel(x_stddev=x)
    #
    #     conv_u = convolve(calib_u, kernel)
    #     conv_g = convolve(calib_g, kernel)
    #     conv_r = convolve(calib_r, kernel)
    #
    #     rgb_u = conv_u
    #     rgb_g = conv_g
    #     rgb_r = conv_r
    #
    #     for q in qs:
    #         for stretch in stretches:
    #             array_rgb = make_lupton_rgb(rgb_r,
    #                                         rgb_g,
    #                                         rgb_u,
    #                                         Q=q,
    #                                         stretch=stretch,
    #                                         minimum=[np.nanmedian(rgb_r),
    #                                                  np.nanmedian(rgb_g),
    #                                                  np.nanmedian(rgb_u)])
    #             plt.imshow(array_rgb, origin='lower')
    #             plt.savefig(ab.work_dir + f'pics/A754_cutout/lupton_test/test_rgb_Q{q}_stretch{stretch}_smooth{x}.png')

    for n in ns:

        stretch_r = mm.jarrett(cutout_r - np.nanmedian(cutout_r), np.nanstd(cutout_r), n)
        stretch_g = mm.jarrett(cutout_g - np.nanmedian(cutout_g), np.nanstd(cutout_g), n)
        stretch_u = mm.jarrett(cutout_u - np.nanmedian(cutout_u), np.nanstd(cutout_u), n)

        # stretch_r[np.isnan(stretch_r)] = 0
        # stretch_g[np.isnan(stretch_g)] = 0
        # stretch_u[np.isnan(stretch_u)] = 0

        plt.imshow(stretch_r,
                   origin = 'lower'
                   )
        plt.savefig(ab.work_dir + f'pics/A754_cutout/lupton_test/jarret_r_n{n}.png')

        plt.imshow(stretch_g,
                   origin = 'lower'
                   )
        plt.savefig(ab.work_dir + f'pics/A754_cutout/lupton_test/jarret_g_n{n}.png')

        plt.imshow(stretch_u,
                   origin = 'lower'
                   )
        plt.savefig(ab.work_dir + f'pics/A754_cutout/lupton_test/jarret_u_n{n}.png')

    plt.close()
