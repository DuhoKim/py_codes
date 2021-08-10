import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')     # https://stackoverflow.com/questions/28757348/how-to-clear-memory-completely-of-all-matplotlib-plots
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
import img_scale
import importlib
importlib.reload(img_scale)
import pylab as py
import sys
from astropy.convolution import Gaussian2DKernel
from scipy.signal import convolve as scipy_convolve
from astropy.convolution import convolve
from astropy.cosmology import Planck18 as Cosmo
from scipy import constants as const

qs = [8 + exp * 4 for exp in range(3, 4)]
# qs = [10]
stretches = [0.08 + 0.02 * exp for exp in range(-2, -1)]
# stretches = [0.08]
x_stddev = [0.5 * exp for exp in range(1, 2)]
n = 3

is_rgb_only = False
is_grey_only = False
is_r_grey_only = False
is_u_g_too = False
save_dir = ab.work_dir + 'pics/A3558_cutout_all/'
# save_dir = ab.work_dir + 'gui/A3558/'
cmap = 'gray'

def tile_check(headers, headers_d, coord, hthumb):
    # read the celestial coordinate
    cel_coord = [[coord.ra.value,
                  coord.dec.value], [0, 0]]
    for jj in range(1, len(headers)):
        # read WCS
        w = wcs.WCS(headers[jj].header)
        pixcrd = w.wcs_world2pix(cel_coord, 1)
        # if w.footprint_contains(sky_coord):
        if (pixcrd[0][0] > hthumb) & (pixcrd[0][0] < headers[jj].shape[1] - hthumb) & \
                (pixcrd[0][1] > hthumb) & (pixcrd[0][1] < headers[jj].shape[0] - hthumb):
            if headers_d[headers[jj].name].data[int(pixcrd[0][1]), int(pixcrd[0][0])] != 2:
                return headers[jj].name, int(pixcrd[0][1]), int(pixcrd[0][0])   # FITS file [y, x]

    return 0, 0, 0

with fits.open(ab.work_dir + 'fits/stacked/A3558_usi.fits') as hdu_u,  \
    fits.open(ab.work_dir + 'fits/stacked/A3558_usd.fits') as hdu_ud,  \
    fits.open(ab.work_dir + 'fits/stacked/A3558_gsi.fits') as hdu_g,   \
    fits.open(ab.work_dir + 'fits/stacked/A3558_gsd.fits') as hdu_gd,   \
    fits.open(ab.work_dir + 'fits/stacked/A3558_rsi.fits') as hdu_r,    \
    fits.open(ab.work_dir + 'fits/stacked/A3558_rsd.fits') as hdu_rd:

    cat = ascii.read("/Users/duhokim/work/abell/spec/Shapley/catalog.dat")

    for i in range(3100, len(cat)):
    # for i in range(3322, 3323):

        c = SkyCoord(f"{cat['col3'][i]}:{cat['col4'][i]}:{cat['col5'][i]} "
                     f"{cat['col6'][i]}:{cat['col7'][i]}:{cat['col8'][i]}", unit=(u.hourangle, u.deg))

        half_thumb = int(np.sqrt(cat['col13'][i])) * 4      # 4 x major axis as a box size
        z = cat['col18'][i]/const.c*1e3                     # redshift
        kpc_arcmin = Cosmo.kpc_proper_per_arcmin(z)         # xxx kpc / arcmin
        arcmin_5kpc = 5.0 / kpc_arcmin                      # xxx arcmin / 5 kpc
        pixel_5kpc = arcmin_5kpc * 60.0 / 0.2637              # xxx pixel / 5 kpc
        # frac_5kpc = arcmin_5kpc * 100.0 / (4 * half_thumb)       # calculate fraction of a length of 5 kpc in fig size

        tile_u, x_u, y_u = tile_check(hdu_u, hdu_ud, c, half_thumb)
        tile_g, x_g, y_g = tile_check(hdu_g, hdu_gd, c, half_thumb)
        tile_r, x_r, y_r = tile_check(hdu_r, hdu_rd, c, half_thumb)

        tile_u_big, x_u_big, y_u_big = tile_check(hdu_u, hdu_ud, c, half_thumb*2)
        tile_g_big, x_g_big, y_g_big = tile_check(hdu_g, hdu_gd, c, half_thumb*2)
        tile_r_big, x_r_big, y_r_big = tile_check(hdu_r, hdu_rd, c, half_thumb*2)

        print(f'{half_thumb} {tile_u} {x_u} {y_u} {tile_g} {x_g} {y_g} {tile_r} {x_r} {y_r}')

        # if is_rgb_only and (not tile_u_big or not tile_g_big or not tile_r_big):
        #     continue

        if tile_r:
            if not is_rgb_only:
                fig = plt.figure(figsize=(5, 5))
            r_img = hdu_r[tile_r].data
            # calib_r = r_img * 10 ** (-ab.stack_a[3][2] / 2.5)
            cutout_r = r_img[x_r - half_thumb: x_r + half_thumb, y_r - half_thumb: y_r + half_thumb] - np.nanmedian(r_img)
            calib_r = cutout_r * 10 ** (-ab.stack_a[3][2] / 2.5)
            if tile_r_big:
                cutout_r_big = r_img[x_r_big - half_thumb*2: x_r_big + half_thumb*2, y_r_big - half_thumb*2: y_r_big + half_thumb*2] - np.nanmedian(
                        r_img)
                calib_r_big = cutout_r_big * 10 ** (-ab.stack_a[3][2] / 2.5)
            if not is_rgb_only:
                stretch_r = mm.jarrett(cutout_r, np.nanstd(cutout_r), n)
                stretch_r[np.isnan(stretch_r)] = 0
                plt.imshow(stretch_r, origin='lower', cmap=cmap)
                plt.savefig(save_dir + f'{i + 1}_r_jar.png')
                stretch_r = img_scale.sqrt(cutout_r)
                stretch_r[np.isnan(stretch_r)] = 0
                plt.imshow(stretch_r, origin='lower', cmap=cmap)
                plt.savefig(save_dir + f'{i + 1}_r_sqrt.png')
                stretch_r = img_scale.asinh(cutout_r)
                stretch_r[np.isnan(stretch_r)] = 0
                plt.imshow(stretch_r, origin='lower', cmap=cmap)
                plt.savefig(save_dir + f'{i + 1}_r_asinh.png')
                plt.close(fig)
            if is_r_grey_only:
                continue
            if tile_g:
                if not is_rgb_only and is_u_g_too:
                    fig = plt.figure(figsize=(5, 5))
                g_img = hdu_g[tile_g].data
                g_array, footprint = reproject_interp(hdu_g[tile_g], hdu_r[tile_r].header)
                cutout_g = g_array[x_r - half_thumb: x_r + half_thumb, y_r - half_thumb: y_r + half_thumb] - np.nanmedian(g_img)
                calib_g = cutout_g * 10 ** (-ab.stack_a[3][1] / 2.5)
                if tile_g_big:
                    cutout_g_big = g_array[x_r_big - half_thumb*2: x_r_big + half_thumb*2,
                               y_r_big - half_thumb*2: y_r_big + half_thumb*2] - np.nanmedian(g_img)
                    calib_g_big = cutout_g_big * 10 ** (-ab.stack_a[3][1] / 2.5)
                if not is_rgb_only and is_u_g_too:
                    stretch_g = mm.jarrett(cutout_g - np.nanmedian(cutout_g), np.nanstd(cutout_g), n)
                    stretch_g[np.isnan(stretch_g)] = 0
                    plt.imshow(stretch_g, origin='lower', cmap=cmap)
                    plt.savefig(save_dir + f'{i + 1}_g.png')
                    plt.close(fig)
                if tile_u:
                    if not is_rgb_only and is_u_g_too:
                        fig = plt.figure(figsize=(5, 5))
                    u_img = hdu_u[tile_u].data
                    u_array, footprint = reproject_interp(hdu_u[tile_u], hdu_r[tile_r].header)
                    cutout_u = u_array[x_r - half_thumb: x_r + half_thumb, y_r - half_thumb: y_r + half_thumb]  - np.nanmedian(u_img)
                    calib_u = cutout_u * 10 ** (-ab.stack_a[3][0] / 2.5)
                    if tile_u_big:
                        cutout_u_big = u_array[x_r_big - half_thumb*2: x_r_big + half_thumb*2,
                                   y_r_big - half_thumb*2: y_r_big + half_thumb*2] - np.nanmedian(u_img)
                        calib_u_big = cutout_u_big * 10 ** (-ab.stack_a[3][0] / 2.5)
                    if not is_rgb_only and is_u_g_too:
                        stretch_u = mm.jarrett(cutout_u - np.nanmedian(cutout_u), np.nanstd(cutout_u), n)
                        stretch_u[np.isnan(stretch_u)] = 0
                        plt.imshow(stretch_u, origin='lower', cmap=cmap)
                        plt.savefig(save_dir + f'{i + 1}_u.png')
                        plt.close(fig)
                    if not is_grey_only:
                        ## test
                        fig = plt.figure(figsize=(5, 5))
                        # img_sqrt_0_1 = np.zeros((calib_r.shape[0], calib_r.shape[1], 3), dtype=float)
                        img_sqrt_0_3 = np.zeros((calib_r.shape[0], calib_r.shape[1], 3), dtype=float)
                        if tile_u_big and tile_g_big and tile_r_big:
                            img_sqrt_0_01 = np.zeros((calib_r_big.shape[0], calib_r_big.shape[1], 3), dtype=float)
                        img_jarret = np.zeros((calib_r.shape[0], calib_r.shape[1], 3), dtype=float)
                        minval = 0
                        maxval = 1
                        # img[:,:,0] = img_scale.sqrt(calib_r)
                        # img[:,:,1] = img_scale.sqrt(calib_g)
                        # img[:,:,2] = img_scale.sqrt(calib_u)
                        if tile_u_big and tile_g_big and tile_r_big:
                            img_sqrt_0_01[:,:,0] = img_scale.sqrt(calib_r_big,
                                                        scale_min=0,
                                                        scale_max=0.1)
                            img_sqrt_0_01[:,:,1] = img_scale.sqrt(calib_g_big,
                                                        scale_min=0,
                                                        scale_max=0.1)
                            img_sqrt_0_01[:,:,2] = img_scale.sqrt(calib_u_big,
                                                        scale_min=0,
                                                        scale_max=0.1)
                        img_sqrt_0_3[:,:,0] = img_scale.sqrt(calib_r,
                                                    scale_min=0,
                                                    scale_max=3)
                        img_sqrt_0_3[:,:,1] = img_scale.sqrt(calib_g,
                                                    scale_min=0,
                                                    scale_max=3)
                        img_sqrt_0_3[:,:,2] = img_scale.sqrt(calib_u,
                                                    scale_min=0,
                                                    scale_max=3)
                        # img_log[:,:,0] = img_scale.log(calib_r,
                        #                                scale_min = 0.01,
                        #                                scale_max = 0.1)
                        # img_log[:,:,1] = img_scale.log(calib_g,
                        #                                scale_min = 0.01,
                        #                                scale_max = 0.1)
                        # img_log[:,:,2] = img_scale.log(calib_u,
                        #                                scale_min = 0.01,
                        #                                scale_max = 0.1)
                        img_jarret[:,:,0] = mm.jarrett(calib_r,
                                                       np.nanstd(calib_r),
                                                       n)
                        img_jarret[:,:,1] = mm.jarrett(calib_g,
                                                       np.nanstd(calib_g),
                                                       n)
                        img_jarret[:,:,2] = mm.jarrett(calib_u,
                                                       np.nanstd(calib_u),
                                                       n)
                        ###
                        # py.clf()

                        plt.imshow(img_sqrt_0_3, aspect='equal', origin='lower')
                        plt.savefig(save_dir + f'{i + 1}_rgb_sqrt_0_3.png')
                        # plt.imshow(img_log, aspect='equal', origin='lower')
                        # plt.savefig(save_dir + f'{i + 1}_rgb_log_001_01.png')
                        plt.imshow(img_jarret, aspect='equal', origin='lower')
                        plt.savefig(save_dir + f'{i + 1}_rgb_jar_{n}.png')
                        if tile_u_big and tile_g_big and tile_r_big:
                            plt.imshow(img_sqrt_0_01, aspect='equal', origin='lower')
                            plt.plot([10, 10 + pixel_5kpc.value], [10, 10], color='white')
                            plt.text(10, 15, f'5kpc at z={z:5.3f}',color='white')
                            plt.savefig(save_dir + f'{i + 1}_rgb_sqrt_0_01_big.png')
                        plt.close(fig)
                        ###
                    else:
                        for q in qs:
                            for stretch in stretches:
                                fig = plt.figure(figsize=(5, 5))
                                array_rgb = make_lupton_rgb(calib_r,
                                                            calib_g,
                                                            calib_u,
                                                            Q=q,
                                                            stretch=stretch)
                                plt.imshow(array_rgb, origin='lower')
                                plt.savefig(save_dir + f'{i + 1}_rgb_Q{q}_str{stretch}.png')
                                plt.close(fig)
                    # for x in x_stddev:
                    #     kernel = Gaussian2DKernel(x_stddev=x)
                    #     conv_u = convolve(calib_u, kernel)
                    #     conv_g = convolve(calib_g, kernel)
                    #     conv_r = convolve(calib_r, kernel)
                    #     rgb_u = conv_u
                    #     rgb_g = conv_g
                    #     rgb_r = conv_r
                    #     for q in qs:
                    #         for stretch in stretches:
                    #             fig = plt.figure(figsize=(5, 5))
                    #             array_rgb = make_lupton_rgb(rgb_r,
                    #                                         rgb_g,
                    #                                         rgb_u,
                    #                                         Q=q,
                    #                                         stretch=stretch)
                    #                                         # minimum=[np.nanmedian(rgb_r),
                    #                                         #          np.nanmedian(rgb_g),
                    #                                         #          np.nanmedian(rgb_u)])
                    #
                    #             plt.imshow(array_rgb, origin='lower')
                    #             plt.savefig(save_dir + f'{i+1}_rgb_Q{q}_str{stretch}_x{x}.png')
                    #             plt.close(fig)
        elif tile_g and not is_rgb_only and not is_r_grey_only and is_u_g_too:
            fig = plt.figure(figsize=(5, 5))
            g_img = hdu_g[tile_g].data
            cutout_g = g_img[x_g - half_thumb: x_g + half_thumb, y_g - half_thumb: y_g + half_thumb]
            stretch_g = mm.jarrett(cutout_g - np.nanmedian(cutout_g), np.nanstd(cutout_g), n)
            stretch_g[np.isnan(stretch_g)] = 0
            plt.imshow(stretch_g, origin='lower', cmap=cmap)
            plt.savefig(save_dir + f'{i + 1}_g.png')
            plt.close(fig)
        elif tile_u and not is_rgb_only and not is_r_grey_only and is_u_g_too:
            fig = plt.figure(figsize=(5, 5))
            u_img = hdu_u[tile_u].data
            cutout_u = u_array[x_u - half_thumb: x_u + half_thumb, y_u - half_thumb: y_u + half_thumb]
            stretch_u = mm.jarrett(cutout_u - np.nanmedian(cutout_u), np.nanstd(cutout_u), n)
            stretch_u[np.isnan(stretch_u)] = 0
            plt.imshow(stretch_u, origin='lower', cmap=cmap)
            plt.savefig(save_dir + f'{i + 1}_u.png')
            plt.close(fig)

