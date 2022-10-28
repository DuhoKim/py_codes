import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
# import matplotlib
# matplotlib.use('TkAgg')     # https://stackoverflow.com/questions/28757348/how-to-clear-memory-completely-of-all-matplotlib-plots
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
importlib.reload(mm)
import pylab as py
import sys
from astropy.convolution import Gaussian2DKernel
from scipy.signal import convolve as scipy_convolve
from astropy.convolution import convolve
from astropy.cosmology import Planck18 as Cosmo
from scipy import constants as const
import os
import calc_kcor as kcor

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

# save_dir = ab.work_dir + 'gui/A3558/'
cmap = 'PuOr'


def tile_check_old(headers, headers_d, coord, hthumb):
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

def tile_check(headers, headers_d, coord, hthumb):
    # read the celestial coordinate
    cel_coord = [[coord.ra.value,
                  coord.dec.value], [0, 0]]
    for jj in range(1, len(headers)):
        # read WCS
        w = wcs.WCS(headers[jj].header)
        pixcrd = w.wcs_world2pix(cel_coord, 1)
        # if w.footprint_contains(sky_coord):

        yy = int(pixcrd[0][0])
        xx = int(pixcrd[0][1])
        xb = headers[jj].shape[0]
        yb = headers[jj].shape[1]
        ht = hthumb

        if (xx > 0) & (xx < xb) & (yy > 0) & (yy < yb):
            y1 = yy - ht if yy - ht > 0 else 0
            y2 = yy + ht if yy + ht < yb else yb
            x1 = xx - ht if xx - ht > 0 else 0
            x2 = xx + ht if xx + ht < xb else xb
            cutout = headers[jj].data[x1:x2, y1:y2]
            # print(f'{ht} {x1} {x2} {y1} {y2} {xb} {yb}')
            return cutout, headers[jj].name, xx, yy   # FITS file [y, x]

    return 0, 0, 0

# for k in range(0, len(ab.clusters)):
for k in range(2, 3):
    save_dir = ab.work_dir + f'pics/{ab.clusters[k]}_Franco/'
    if not os.path.isdir(save_dir):
        os.makedirs(save_dir)

    with fits.open(ab.work_dir + f'fits/stacked/{ab.clusters[k]}_usi.fits') as hdu_u,  \
        fits.open(ab.work_dir + f'fits/stacked/{ab.clusters[k]}_usd.fits') as hdu_ud,  \
        fits.open(ab.work_dir + f'fits/stacked/{ab.clusters[k]}_gsi.fits') as hdu_g,   \
        fits.open(ab.work_dir + f'fits/stacked/{ab.clusters[k]}_gsd.fits') as hdu_gd,   \
        fits.open(ab.work_dir + f'fits/stacked/{ab.clusters[k]}_rsi.fits') as hdu_r,    \
        fits.open(ab.work_dir + f'fits/stacked/{ab.clusters[k]}_rsd.fits') as hdu_rd, \
        open(ab.work_dir + f'spec/{ab.clusters[k]}_spec_ind.npy', 'rb') as spec_ind, \
        open(ab.work_dir+f'spec/{ab.clusters[k]}_rs_each_ind.npy', 'rb') as rs_ind,  \
        open(save_dir+f'{ab.clusters[k]}_rs_or_spec_tmp.vis', 'w') as vis:

        # cat = ascii.read("/Users/duhokim/work/abell/spec/Shapley/catalog.dat")
        # cat = ascii.read("/Users/duhokim/work/abell/cat/sheen_2012_table_6.txt")
        sex_cat = ascii.read(ab.sex_dir + f'DECam_merged_SEx_cat_{ab.clusters[k]}_Gal_ext_corrected_{ab.ver}.txt')
        ind_spec = np.load(spec_ind)
        ind_rs = np.load(rs_ind)
        ind_xor = ind_spec ^ ind_rs
        ind_spec_only = ind_xor & ind_spec
        ind_spec_or_rs = ind_spec | ind_rs
        cat_rs = sex_cat[ind_spec_or_rs]

        # to calculate M_r < -20 limit
        sex_cat[ab.mag_sys + '_g'] -= kcor.calc_kcor('g', ab.redshifts[k], 'g - r',
                                                     sex_cat[ab.mag_sys + '_g'] - sex_cat[ab.mag_sys + '_r'])
        sex_cat[ab.mag_sys + '_r'] -= kcor.calc_kcor('r', ab.redshifts[k], 'g - r',
                                                     sex_cat[ab.mag_sys + '_g'] - sex_cat[ab.mag_sys + '_r'])
        r_abs_kcor = sex_cat[ab.mag_sys + '_r'] - ab.distmod[k]
        ind_mag_cut = r_abs_kcor < -20
        cat_rs = sex_cat[ind_mag_cut]

        ####### FOR TEMPOPARY USE ONLY ###########
        # select color region only for
        color_rs = ["" for x in range(len(cat_rs))]
        for i in range(0, len(cat_rs)):
            if cat_rs[ab.mag_sys + '_g'][i] - cat_rs[ab.mag_sys + '_r'][i] > -0.03197 * cat_rs[
                ab.mag_sys + '_r'][i] + 1.275 + 0.1188:
                color_rs[i] = 'ER'     # 'ER' extemely red
            elif (cat_rs[ab.mag_sys + '_g'][i] - cat_rs[ab.mag_sys + '_r'][i] < -0.03197 * cat_rs[
                ab.mag_sys + '_r'][i] + 1.275 + 0.1188) & \
                    (cat_rs[ab.mag_sys + '_g'][i] - cat_rs[ab.mag_sys + '_r'][i] > -0.03197 * cat_rs[
                        ab.mag_sys + '_r'][i] + 1.275 - 0.1188):
                color_rs[i] = 'RS'     # 'RS'  red sequence
            else:
                color_rs[i] = 'BC'     # 'BC'  blue cloud
        ####### FOR TEMPOPARY USE ONLY ###########

        nan_id = []
        # for i in range(0, 1):
        for i in range(0, len(cat_rs)):
            vis.writelines(f'{cat_rs["NUMBER"][i]}\n')
            c = SkyCoord(f"{cat_rs['ALPHA_J2000'][i]} {cat_rs['DELTA_J2000'][i]}", unit='deg')
            if os.path.exists(save_dir + f'{cat_rs["NUMBER"][i]}_rgb_sqrt_0_001_big.png'):
                mm.copy_gala_result(ab.clusters[k], hdu_r, c, save_dir, cat_rs["NUMBER"][i])
                continue
            half_thumb = int(10 * cat_rs['A_WORLD'][i] / 0.2637)  # 4 x major axis as a box size
            z = ab.redshifts[k]  # redshift
            kpc_arcmin = Cosmo.kpc_proper_per_arcmin(z)         # xxx kpc / arcmin
            arcmin_5kpc = 5.0 / kpc_arcmin                      # xxx arcmin / 5 kpc
            pixel_5kpc = arcmin_5kpc * 60.0 / 0.2637              # xxx pixel / 5 kpc
            # frac_5kpc = arcmin_5kpc * 100.0 / (4 * half_thumb)       # calculate fraction of a length of 5 kpc in fig size

            sky_u = mm.get_gala_sky(ab.clusters[k], 'u', hdu_u, c)
            sky_g = mm.get_gala_sky(ab.clusters[k], 'g', hdu_g, c)
            sky_r = mm.get_gala_sky(ab.clusters[k], 'r', hdu_r, c)

            if np.isnan(sky_u) or np.isnan(sky_g) or np.isnan(sky_r):
                print(f'{ab.clusters[k]} sky is nan: {cat_rs["NUMBER"][i]}')
                nan_id.append(f'{cat_rs["NUMBER"][i]}')
                continue

            cutout_u, tile_u, x_u, y_u = tile_check(hdu_u, hdu_ud, c, half_thumb)
            cutout_g, tile_g, x_g, y_g = tile_check(hdu_g, hdu_gd, c, half_thumb)
            cutout_r, tile_r, x_r, y_r = tile_check(hdu_r, hdu_rd, c, half_thumb)

            cutout_u_big, tile_u_big, x_u_big, y_u_big = tile_check(hdu_u, hdu_ud, c, half_thumb*2)
            cutout_g_big, tile_g_big, x_g_big, y_g_big = tile_check(hdu_g, hdu_gd, c, half_thumb*2)
            cutout_r_big, tile_r_big, x_r_big, y_r_big = tile_check(hdu_r, hdu_rd, c, half_thumb*2)
            # print(f'{half_thumb} {tile_u} {x_u} {y_u} {tile_g} {x_g} {y_g} {tile_r} {x_r} {y_r} {sky_u} {sky_g} {sky_r}')
            # if is_rgb_only and (not tile_u_big or not tile_g_big or not tile_r_big):
            #     continue

            if tile_r:
                if not is_rgb_only:
                    fig = plt.figure(figsize=(5, 5))
                cutout_r = cutout_r - sky_r
                calib_r = cutout_r * 10 ** (-ab.stack_a[k][2] / 2.5)
                if tile_r_big:
                    cutout_r_big = cutout_r_big - sky_r
                    calib_r_big = cutout_r_big * 10 ** (-ab.stack_a[k][2] / 2.5)
                if not is_rgb_only:
                    stretch_r = mm.jarrett(cutout_r, np.nanstd(cutout_r), n)
                    stretch_r[np.isnan(stretch_r)] = 0
                    plt.imshow(stretch_r, origin='lower', cmap=cmap)
                    plt.savefig(save_dir + f'{cat_rs["NUMBER"][i]}_r_jar.png')
                    stretch_r = img_scale.sqrt(cutout_r, scale_min=0, scale_max=50)
                    stretch_r[np.isnan(stretch_r)] = 0
                    plt.imshow(stretch_r, origin='lower', cmap=cmap)
                    plt.savefig(save_dir + f'{cat_rs["NUMBER"][i]}_r_sqrt_50.png')
                    stretch_r = img_scale.asinh(cutout_r, scale_min=0, scale_max=100)
                    stretch_r[np.isnan(stretch_r)] = 0
                    plt.imshow(stretch_r, origin='lower', cmap=cmap)
                    plt.savefig(save_dir + f'{cat_rs["NUMBER"][i]}_r_asinh_100.png')
                    plt.close(fig)
                if is_r_grey_only:
                    continue
                if tile_g:
                    if not is_rgb_only and is_u_g_too:
                        fig = plt.figure(figsize=(5, 5))
                    g_array, footprint = reproject_interp(hdu_g[tile_g], hdu_r[tile_r].header)
                    cutout_g = g_array[ x_r - half_thumb if x_r - half_thumb > 0 else 0:
                                        x_r + half_thumb if x_r + half_thumb < g_array.shape[0] else g_array.shape[0], \
                                        y_r - half_thumb if y_r - half_thumb > 0 else 0:
                                        y_r + half_thumb if y_r + half_thumb < g_array.shape[1] else g_array.shape[1]]- sky_g
                    calib_g = cutout_g * 10 ** (-ab.stack_a[k][1] / 2.5)
                    if tile_g_big:
                        cutout_g_big = g_array[ x_r - 2*half_thumb if x_r - 2*half_thumb > 0 else 0:
                                                x_r + 2*half_thumb if x_r + 2*half_thumb < g_array.shape[0] else g_array.shape[0], \
                                                y_r - 2*half_thumb if y_r - 2*half_thumb > 0 else 0:
                                                y_r + 2*half_thumb if y_r + 2*half_thumb < g_array.shape[1] else g_array.shape[1]]- sky_g
                        calib_g_big = cutout_g_big * 10 ** (-ab.stack_a[k][1] / 2.5)
                    if not is_rgb_only and is_u_g_too:
                        stretch_g = mm.jarrett(cutout_g, np.nanstd(cutout_g), n)
                        stretch_g[np.isnan(stretch_g)] = 0
                        plt.imshow(stretch_g, origin='lower', cmap=cmap)
                        plt.savefig(save_dir + f'{cat_rs["NUMBER"][i]}_g.png')
                        plt.close(fig)
                    if tile_u:
                        if not is_rgb_only and is_u_g_too:
                            fig = plt.figure(figsize=(5, 5))
                        u_img = hdu_u[tile_u].data
                        u_array, footprint = reproject_interp(hdu_u[tile_u], hdu_r[tile_r].header)
                        cutout_u = u_array[ x_r - half_thumb if x_r - half_thumb > 0 else 0:
                                            x_r + half_thumb if x_r + half_thumb < u_array.shape[0] else u_array.shape[0], \
                                            y_r - half_thumb if y_r - half_thumb > 0 else 0:
                                            y_r + half_thumb if y_r + half_thumb < u_array.shape[1] else u_array.shape[1]]- sky_u
                        calib_u = cutout_u * 10 ** (-ab.stack_a[k][0] / 2.5)
                        if tile_u_big:
                            cutout_u_big = u_array[ x_r - 2*half_thumb if x_r - 2*half_thumb > 0 else 0:
                                                    x_r + 2*half_thumb if x_r + 2*half_thumb < u_array.shape[0] else u_array.shape[0], \
                                                    y_r - 2*half_thumb if y_r - 2*half_thumb > 0 else 0:
                                                    y_r + 2*half_thumb if y_r + 2*half_thumb < u_array.shape[1] else u_array.shape[1]]- sky_u
                            calib_u_big = cutout_u_big * 10 ** (-ab.stack_a[k][0] / 2.5)
                        if not is_rgb_only and is_u_g_too:
                            stretch_u = mm.jarrett(cutout_u, np.nanstd(cutout_u), n)
                            stretch_u[np.isnan(stretch_u)] = 0
                            plt.imshow(stretch_u, origin='lower', cmap=cmap)
                            plt.savefig(save_dir + f'{cat_rs["NUMBER"][i]}_u.png')
                            plt.close(fig)
                        if not is_grey_only:
                            ## test
                            fig = plt.figure(figsize=(5, 5))
                            # img_sqrt_0_1 = np.zeros((calib_r.shape[0], calib_r.shape[1], 3), dtype=float)
                            img_sqrt_0_3 = np.zeros((calib_r.shape[0], calib_r.shape[1], 3), dtype=float)
                            if tile_u_big and tile_g_big and tile_r_big:
                                img_sqrt_0_01 = np.zeros((calib_r_big.shape[0], calib_r_big.shape[1], 3), dtype=float)
                                img_sqrt_0_001 = np.zeros((calib_r_big.shape[0], calib_r_big.shape[1], 3), dtype=float)
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
                                img_sqrt_0_001[:,:,0] = img_scale.sqrt(calib_r_big,
                                                            scale_min=0,
                                                            scale_max=0.01)
                                img_sqrt_0_001[:,:,1] = img_scale.sqrt(calib_g_big,
                                                            scale_min=0,
                                                            scale_max=0.01)
                                img_sqrt_0_001[:,:,2] = img_scale.sqrt(calib_u_big,
                                                            scale_min=0,
                                                            scale_max=0.01)
                            img_sqrt_0_3[:,:,0] = img_scale.sqrt(calib_r,
                                                        scale_min=0,
                                                        scale_max=1)
                            img_sqrt_0_3[:,:,1] = img_scale.sqrt(calib_g,
                                                        scale_min=0,
                                                        scale_max=1)
                            img_sqrt_0_3[:,:,2] = img_scale.sqrt(calib_u,
                                                        scale_min=0,
                                                        scale_max=1)
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
                            plt.savefig(save_dir + f'{cat_rs["NUMBER"][i]}_rgb_sqrt_0_1.png')
                            # plt.imshow(img_log, aspect='equal', origin='lower')
                            # plt.savefig(save_dir + f'{i + 1}_rgb_log_001_01.png')
                            plt.imshow(img_jarret, aspect='equal', origin='lower')
                            ### Add Text 'ER', 'RS', 'BC'
                            plt.text(10, 15, color_rs[i], color='white')
                            plt.savefig(save_dir + f'{cat_rs["NUMBER"][i]}_rgb_jar_{n}.png')
                            if tile_u_big and tile_g_big and tile_r_big:
                                plt.imshow(img_sqrt_0_01, aspect='equal', origin='lower')
                                plt.plot([10, 10 + pixel_5kpc.value], [10, 10], color='white')
                                plt.text(10, 15, f'5kpc at z={z:5.3f}',color='white')
                                plt.savefig(save_dir + f'{cat_rs["NUMBER"][i]}_rgb_sqrt_0_01_big.png')
                                plt.imshow(img_sqrt_0_001, aspect='equal', origin='lower', cmap=cmap)
                                # plt.plot([10, 10 + pixel_5kpc.value], [10, 10], color='white')
                                # plt.text(10, 15, f'5kpc at z={z:5.3f}',color='white')
                                plt.savefig(save_dir + f'{cat_rs["NUMBER"][i]}_rgb_sqrt_0_001_big.png')
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
                                    plt.savefig(save_dir + f'{cat_rs["NUMBER"][i]}_rgb_Q{q}_str{stretch}.png')
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
                plt.savefig(save_dir + f'{cat_rs["NUMBER"][i]}_g.png')
                plt.close(fig)
            elif tile_u and not is_rgb_only and not is_r_grey_only and is_u_g_too:
                fig = plt.figure(figsize=(5, 5))
                u_img = hdu_u[tile_u].data
                cutout_u = u_array[x_u - half_thumb: x_u + half_thumb, y_u - half_thumb: y_u + half_thumb]
                stretch_u = mm.jarrett(cutout_u - np.nanmedian(cutout_u), np.nanstd(cutout_u), n)
                stretch_u[np.isnan(stretch_u)] = 0
                plt.imshow(stretch_u, origin='lower', cmap=cmap)
                plt.savefig(save_dir + f'{cat_rs["NUMBER"][i]}_u.png')
                plt.close(fig)

    with open(save_dir+f'{ab.clusters[k]}_rs_or_spec_tmp.vis', 'r') as f:
        lines = f.readlines()

    with open(save_dir+f'{ab.clusters[k]}_rs_or_spec.vis', 'w') as f:
        f.write('0\n')      # first line is reserved for the current index
        for line in lines:
            if line.strip('\n') not in nan_id:
                f.write(line)
