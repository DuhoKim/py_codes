import numpy as np
import my_module as mm
from astropy.io import fits
import matplotlib.pyplot as plt
from reproject.mosaicking import find_optimal_celestial_wcs, reproject_and_coadd
from reproject import reproject_interp
from astropy.visualization.wcsaxes import SphericalCircle
from astropy import units as u
import abell_cluster_module as ab

for k in range(0, len(ab.clusters)):
    # with fits.open(work_dir + 'fits/stacked/' + clusters[k] + '_rsi.fits') as hdu_r:
    with fits.open(ab.work_dir + 'fits/best_single/' + ab.short_sci_fn[k][2]) as hdu_r:
        plt.figure(figsize=(12, 12))
        wcs_out, shape_out = find_optimal_celestial_wcs(hdu_r[1:len(hdu_r)])  # has only CompImageHDU files
        array, footprint = reproject_and_coadd(hdu_r[1:len(hdu_r)],
                                               wcs_out,
                                               shape_out=shape_out,
                                               reproject_function=reproject_interp)

        ax = plt.gca(projection=wcs_out)
        plt.imshow(mm.jarrett(array-np.nanmedian(array), np.nanstd(array), 5), cmap='gray_r')
        plt.xlabel('R.A.')
        plt.ylabel('Decl.')
        c = SphericalCircle((ab.coords_cl_cen[k].ra.value, ab.coords_cl_cen[k].dec.value) * u.deg,
                            0.8 * u.deg,
                            edgecolor='gray',
                            facecolor='none',
                            linestyle='--',
                            transform=ax.get_transform('fk5'))
        ax.add_patch(c)

        plt.savefig(ab.work_dir + 'plot/' + ab.clusters[k] + '_short_mosaic_sub_n5.pdf')
        plt.show()