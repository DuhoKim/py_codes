import numpy as np
import my_module as mm
from astropy.io import fits
import matplotlib.pyplot as plt
from reproject.mosaicking import find_optimal_celestial_wcs, reproject_and_coadd
from reproject import reproject_interp
from astropy.visualization.wcsaxes import SphericalCircle, Quadrangle
from astropy import units as u
import abell_cluster_module as ab
from astropy.cosmology import Planck18 as Cosmo
import img_scale
# from matplotlib.patches import Rectangle

# for k in range(1, len(ab.clusters)):
for k in range(2, 3):
    mm.print_time()

    with fits.open(ab.work_dir + 'fits/best_single/' + ab.short_sci_fn[k][2]) as hdu_r:
        plt.figure(figsize=(20, 20))
        wcs_out, shape_out = find_optimal_celestial_wcs(hdu_r[1:len(hdu_r)])  # has only CompImageHDU files
        array, footprint = reproject_and_coadd(hdu_r[1:len(hdu_r)],
                                               wcs_out,
                                               shape_out=shape_out,
                                               reproject_function=reproject_interp)

        ax = plt.gca(projection=wcs_out)

        # plt.imshow(mm.jarrett(array-sky_r, np.nanstd(array), 5), cmap='gray_r')
        stretch_r = img_scale.sqrt(array - np.nanmedian(array), scale_min=0, scale_max=100)
        # plt.imshow(stretch_r, origin='lower', cmap='PuOr')
        plt.imshow(stretch_r, origin='lower', cmap='binary')
        plt.xlabel('R.A.', fontsize=30)
        plt.ylabel('Decl.', fontsize=30)
        ax.tick_params(axis='both', which='major', labelsize=25)
        # ax.tick_params(axis='both', which='major', labelsize=20)
        # c = SphericalCircle((ab.coords_cl_cen[k].ra.value, ab.coords_cl_cen[k].dec.value) * u.deg,
        #                     0.8 * u.deg,
        #                     edgecolor='purple',
        #                     facecolor='none',
        #                     linestyle='--',
        #                     transform=ax.get_transform('fk5'),
        #                     linewidth=3)
        # ax.add_patch(c)

        z = ab.redshifts[k]  # redshift
        kpc_arcmin = Cosmo.kpc_proper_per_arcmin(z)  # xxx kpc / arcmin
        arcmin_Mpc = 1e3 / kpc_arcmin  # xxx arcmin / 5 Mpc
        pixel_Mpc = arcmin_Mpc * 60.0 / 0.2637  # xxx pixel / 5 Mpc

        # h2.text(0.1, 0.9, f'Sersic index = {ser_idx}', color='white', transform=h2.transAxes)

        plt.plot([1000, 1000 + pixel_Mpc.value], [1000, 1000], linewidth=3, color='purple')
        plt.text(1000, 1500, f'1Mpc at z={z:5.3f}', fontsize=25)

        plt.savefig(ab.work_dir + 'plot/' + ab.clusters[k] + '_short_mosaic_sqrt_binary.pdf')

    with fits.open(ab.work_dir + 'fits/stacked/' + ab.clusters[k] + '_rsi.fits') as hdu_r:
        sky_r = mm.get_gala_sky(ab.clusters[k], 'r', hdu_r, ab.coords_cl_cen[k])

        plt.figure(figsize=(20, 20))
        wcs_out, shape_out = find_optimal_celestial_wcs(hdu_r[1:len(hdu_r)])  # has only CompImageHDU files
        array, footprint = reproject_and_coadd(hdu_r[1:len(hdu_r)],
                                               wcs_out,
                                               shape_out=shape_out,
                                               reproject_function=reproject_interp)

        ax = plt.gca(projection=wcs_out)

        # plt.imshow(mm.jarrett(array-sky_r, np.nanstd(array), 5), cmap='gray_r')
        stretch_r = img_scale.sqrt(array - sky_r, scale_min=0, scale_max=100)
        # plt.imshow(stretch_r, origin='lower', cmap='PuOr')
        plt.imshow(stretch_r, origin='lower', cmap='binary')
        plt.xlabel('R.A.', fontsize=30)
        plt.ylabel('Decl.', fontsize=30)
        ax.tick_params(axis='both', which='major', labelsize=25)
        # ax.tick_params(axis='both', which='major', labelsize=20)

        r = Quadrangle((358.15, -10.65)*u.deg, 0.65*u.deg, 0.65*u.deg,
                      edgecolor='black', facecolor='none',
                      transform=ax.get_transform('fk5'))
        ax.add_patch(r)


        z = ab.redshifts[k]  # redshift
        kpc_arcmin = Cosmo.kpc_proper_per_arcmin(z)  # xxx kpc / arcmin
        arcmin_Mpc = 1e3 / kpc_arcmin  # xxx arcmin / 5 Mpc
        pixel_Mpc = arcmin_Mpc * 60.0 / 0.2637  # xxx pixel / 5 Mpc

        c = SphericalCircle((ab.coords_cl_cen[k].ra.value, ab.coords_cl_cen[k].dec.value) * u.deg,
                            ab.r200[k] * arcmin_Mpc.value * u.arcmin,
                            edgecolor='black',
                            facecolor='none',
                            linestyle='--',
                            transform=ax.get_transform('fk5'),
                            linewidth=3)
        ax.add_patch(c)

        # h2.text(0.1, 0.9, f'Sersic index = {ser_idx}', color='white', transform=h2.transAxes)

        plt.plot([1000, 1000 + pixel_Mpc.value], [1000, 1000], linewidth=3, color='purple')
        plt.text(1000, 1500, f'1Mpc at z={z:5.3f}', fontsize=25)

        plt.savefig(ab.work_dir + 'plot/' + ab.clusters[k] + '_stack_mosaic_sqrt_binary.pdf')
