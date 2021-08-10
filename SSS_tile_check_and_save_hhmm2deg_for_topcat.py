from astropy.io import fits
from astropy import units as u
import abell_cluster_module as ab
from astropy.coordinates import SkyCoord
from astropy.io import ascii
from astropy import wcs
import my_module as mm
from scipy import constants as const

fn = ab.work_dir+'spec/Shapley/inside_3574.cat'

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

print('start!')
mm.print_time()

with fits.open(ab.work_dir + 'fits/stacked/A3574_usi.fits') as hdu_u,  \
    fits.open(ab.work_dir + 'fits/stacked/A3574_usd.fits') as hdu_ud,  \
    fits.open(ab.work_dir + 'fits/stacked/A3574_gsi.fits') as hdu_g,   \
    fits.open(ab.work_dir + 'fits/stacked/A3574_gsd.fits') as hdu_gd,   \
    fits.open(ab.work_dir + 'fits/stacked/A3574_rsi.fits') as hdu_r,    \
    fits.open(ab.work_dir + 'fits/stacked/A3574_rsd.fits') as hdu_rd,   \
    open(fn, 'w') as f:

    cat = ascii.read("/Users/duhokim/work/abell/spec/Shapley/catalog.dat")

    for i in range(0, len(cat)):
    #for i in range(3200, 3300):

        c = SkyCoord(f"{cat['col3'][i]}:{cat['col4'][i]}:{cat['col5'][i]} "
                     f"{cat['col6'][i]}:{cat['col7'][i]}:{cat['col8'][i]}", unit=(u.hourangle, u.deg))

        # half_thumb = int(np.sqrt(cat['col13'][i])) * 4  # 4 x major axis as a box size
        half_thumb = 1  # This is for membership check not for rgb cutout

        tile_u, x_u, y_u = tile_check(hdu_u, hdu_ud, c, half_thumb)
        tile_g, x_g, y_g = tile_check(hdu_g, hdu_gd, c, half_thumb)
        tile_r, x_r, y_r = tile_check(hdu_r, hdu_rd, c, half_thumb)

        # print(f'{half_thumb} {tile_u} {x_u} {y_u} {tile_g} {x_g} {y_g} {tile_r} {x_r} {y_r}')

        if tile_u or tile_g or tile_r:
            f.write(f"{cat['col1'][i]} {c.ra.value} {c.dec.value} {cat['col18'][i]/const.c*1e3} \n")

print('end!')
mm.print_time()