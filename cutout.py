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

# for k in range(0, len(ab.clusters)):
for k in range(0, 1):
    with fits.open(ab.work_dir + 'fits/stacked/' + ab.clusters[k] + '_rsi.fits') as hdu_r:

        cat = ascii.read("/Users/duhokim/work/abell/cat/WINGS_Cava+2009_A754.txt")
        sex_cat_r = ascii.read(ab.short_cat_dir + ab.short_cat_fn[k][2] + '_deblend.cat')

        coords_r = SkyCoord(sex_cat_r['ALPHA_J2000'], sex_cat_r['DELTA_J2000'], unit='deg')

        for i in range(0, len(cat)):
            c = SkyCoord(f"{cat['col8'][i]}:{cat['col9'][i]}:{cat['col10'][i]} {cat['col11'][i]}:{cat['col12'][i]}:{cat['col13'][i]}",
                         unit=(u.hourangle, u.deg))

            idx_r2c, d2d_c2r, d3d = c.match_to_catalog_sky(coords_r)
            if d2d_c2r.arcsec < 1:
                thumb = int(sex_cat_r['A_WORLD'][idx_r2c] * 3600 / 0.2635) * 10
            else:
                thumb = 100
            cel_coord = [[c.ra.value, c.dec.value], [0, 0]]

            for jj in range(1, len(hdu_r)):
                w = wcs.WCS(hdu_r[jj].header)
                pixcrd = w.wcs_world2pix(cel_coord, 1)

                if (pixcrd[0][0] > thumb) & (pixcrd[0][0] < hdu_r[jj].shape[1]-thumb) & \
                        (pixcrd[0][1] > thumb) & (pixcrd[0][1] < hdu_r[jj].shape[0]-thumb):
                    r_img = hdu_r[jj].data
                    array = r_img[int(pixcrd[0][1]) - thumb : int(pixcrd[0][1]) + thumb,
                            int(pixcrd[0][0]) - thumb : int(pixcrd[0][0]) + thumb]

                    plt.figure(figsize=(5, 5))

                    plt.imshow(mm.jarrett(array-np.nanmedian(array), np.nanstd(array), 0.05), cmap='gray_r')
                    plt.savefig(ab.work_dir + 'pics/A754_cutout/'+str(i)+'_n005.png')
                    plt.close()

