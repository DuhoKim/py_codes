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

# qs = [10 + exp * 2 for exp in range(1, 4)]
qs = [10]
stretches = [0.08 * exp for exp in range(1, 2)]
x_stddev = [1.0 * exp for exp in range(1, 2)]
n = 7

# for k in range(0, len(ab.clusters)):
for k in range(0, 1):
    with fits.open(ab.work_dir + 'fits/stacked/' + ab.clusters[k] + '_usi.fits') as hdu_u,  \
        fits.open(ab.work_dir + 'fits/stacked/' + ab.clusters[k] + '_gsi.fits') as hdu_g,   \
        fits.open(ab.work_dir + 'fits/stacked/' + ab.clusters[k] + '_rsi.fits') as hdu_r:

        cat = ascii.read("/Users/duhokim/work/abell/cat/WINGS_Cava+2009_A754.txt")
        sex_cat_r = ascii.read(ab.short_cat_dir + ab.short_cat_fn[k][2] + '_deblend.cat')

        coords_r = SkyCoord(sex_cat_r['ALPHA_J2000'], sex_cat_r['DELTA_J2000'], unit='deg')

        for i in range(0, len(cat)):
        # for i in range(1, 2):
            plt.figure(figsize=(5, 5))
            c = SkyCoord(f"{cat['col8'][i]}:{cat['col9'][i]}:{cat['col10'][i]} {cat['col11'][i]}:{cat['col12'][i]}:{cat['col13'][i]}",
                         unit=(u.hourangle, u.deg))

            idx_r2c, d2d_c2r, d3d = c.match_to_catalog_sky(coords_r)
            if d2d_c2r.arcsec < 1:
                thumb = int(sex_cat_r['A_WORLD'][idx_r2c] * 3600 / 0.2635) * 15     # cutout x15 of major axis
            else:
                thumb = 150
            cel_coord = [[c.ra.value, c.dec.value], [0, 0]]

            for jj in range(1, len(hdu_r)):
                w = wcs.WCS(hdu_r[jj].header)
                pixcrd = w.wcs_world2pix(cel_coord, 1)
                if (pixcrd[0][0] > thumb) & (pixcrd[0][0] < hdu_r[jj].shape[1]-thumb) & \
                        (pixcrd[0][1] > thumb) & (pixcrd[0][1] < hdu_r[jj].shape[0]-thumb):
                    r_img = hdu_r[jj].data
                    cutout_r = r_img[int(pixcrd[0][1]) - thumb: int(pixcrd[0][1]) + thumb,
                               int(pixcrd[0][0]) - thumb: int(pixcrd[0][0]) + thumb]
                    calib_r = cutout_r * 10 ** (-ab.stack_a[0][2] / 2.5)
                    stretch_r = mm.jarrett(cutout_r - np.nanmedian(cutout_r), np.nanstd(cutout_r), n)
                    stretch_r[np.isnan(stretch_r)] = 0
                    plt.imshow(stretch_r, origin='lower')
                    plt.savefig(ab.work_dir + f'pics/A754_cutout/best/{i+1}_r.png')

                    for gg in range(1, len(hdu_g)):
                        g_w = wcs.WCS(hdu_g[gg].header)
                        g_pixcrd = g_w.wcs_world2pix(cel_coord, 1)
                        if (g_pixcrd[0][0] > thumb) & (g_pixcrd[0][0] < hdu_g[gg].shape[1] - thumb) & \
                                (g_pixcrd[0][1] > thumb) & (g_pixcrd[0][1] < hdu_g[gg].shape[0] - thumb):
                            g_img = hdu_g[gg].data
                            g_array, footprint = reproject_interp(hdu_g[gg], hdu_r[jj].header)
                            cutout_g = g_array[int(pixcrd[0][1]) - thumb: int(pixcrd[0][1]) + thumb,
                                       int(pixcrd[0][0]) - thumb: int(pixcrd[0][0]) + thumb]
                            calib_g = cutout_g * 10 ** (-ab.stack_a[0][1] / 2.5)
                            stretch_g = mm.jarrett(cutout_g - np.nanmedian(cutout_g), np.nanstd(cutout_g), n)
                            stretch_g[np.isnan(stretch_g)] = 0
                            plt.imshow(stretch_g, origin='lower')
                            plt.savefig(ab.work_dir + f'pics/A754_cutout/best/{i+1}_g.png')

                            for uu in range(1, len(hdu_u)):
                                u_w = wcs.WCS(hdu_u[uu].header)
                                u_pixcrd = u_w.wcs_world2pix(cel_coord, 1)
                                if (u_pixcrd[0][0] > thumb) & (u_pixcrd[0][0] < hdu_u[uu].shape[1] - thumb) & \
                                        (u_pixcrd[0][1] > thumb) & (u_pixcrd[0][1] < hdu_u[uu].shape[0] - thumb):
                                    u_img = hdu_u[uu].data
                                    u_array, footprint = reproject_interp(hdu_u[uu], hdu_r[jj].header)
                                    cutout_u = u_array[int(pixcrd[0][1]) - thumb: int(pixcrd[0][1]) + thumb,
                                             int(pixcrd[0][0]) - thumb: int(pixcrd[0][0]) + thumb]
                                    calib_u = cutout_u * 10 ** (-ab.stack_a[0][0] / 2.5)
                                    stretch_u = mm.jarrett(cutout_u - np.nanmedian(cutout_u), np.nanstd(cutout_u), n)
                                    stretch_u[np.isnan(stretch_u)] = 0
                                    plt.imshow(stretch_u, origin='lower')
                                    plt.savefig(ab.work_dir + f'pics/A754_cutout/best/{i+1}_u_.png')

                                    for x in x_stddev:
                                        kernel = Gaussian2DKernel(x_stddev=x)
                                        conv_u = convolve(calib_u, kernel)
                                        conv_g = convolve(calib_g, kernel)
                                        conv_r = convolve(calib_r, kernel)
                                        rgb_u = conv_u
                                        rgb_g = conv_g
                                        rgb_r = conv_r
                                        for q in qs:
                                            for stretch in stretches:
                                                array_rgb = make_lupton_rgb(rgb_r,
                                                                            rgb_g,
                                                                            rgb_u,
                                                                            Q=q,
                                                                            stretch=stretch,
                                                                            minimum=[np.nanmedian(rgb_r),
                                                                                     np.nanmedian(rgb_g),
                                                                                     np.nanmedian(rgb_u)])

                                                plt.imshow(array_rgb, origin='lower')
                                                plt.savefig(
                                                    ab.work_dir + f'pics/A754_cutout/best/{i+1}_rgb.png')
            plt.close()

