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
import importlib
importlib.reload(ab)

# number_of_pixel for each tile is ~ 10,000 so number of mesh will be ~10 x 10 for each tile
mesh_sizes = [100, 500, 1000, 2000, 5000]

for k in range(0, len(ab.clusters)):
# for k in range(0, 1):
    print(ab.clusters[k] + '\'s surface brightness limit is:')
    with fits.open(ab.work_dir + 'fits/stacked/' + ab.clusters[k] + '_rsi.fits') as hdu_r, \
        fits.open(ab.check_dir + ab.clusters[k] + '_rsi_seg.fits') as hdu_r_seg:
        for mesh in mesh_sizes:
            print(f'For mesh size : {mesh}')
            sigmas = []
            for jj in range(1, len(hdu_r_seg)):
                w = wcs.WCS(hdu_r[jj].header)
                tile = hdu_r[jj].data
                seg_tile = hdu_r_seg[jj].data

                ind = np.where(seg_tile != 0)

                tile[ind] = np.nan

                x_size = w.array_shape[0]
                y_size = w.array_shape[1]

                nx = int(x_size / mesh)
                ny = int(y_size / mesh)

                for xi in range(0, nx):
                    for yi in range(0, ny):
                        if (xi + 1) * mesh < x_size and (yi + 1) * mesh < y_size:
                            # opposite convention between FITS & Python array indexing
                            sigmas.append(np.nanstd(tile[yi * mesh : (yi+1) * mesh, xi * mesh : (xi+1) * mesh]))
            print(f'sigma : {np.nanmedian(sigmas)}')
            print(f'SB = {-2.5 * np.log10(np.nanmedian(sigmas) * np.sqrt(1/(0.2637)**2)) + 25.0 + ab.stack_a[k][2]} [mag/arcsec^2]')

