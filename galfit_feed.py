from astropy.io import ascii
from astropy.wcs import WCS
from astropy import wcs
import abell_cluster_module as ab
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np

hdu_r_stack = fits.open(ab.work_dir + 'fits/stacked/A3558_rsi.fits')
cat = ascii.read("/Users/duhokim/work/abell/spec/Shapley/catalog.dat")

def tile_check(coord, hthumb):
    # read the celestial coordinate
    cel_coord = [[coord.ra.value,
                  coord.dec.value], [0, 0]]
    for jj in range(1, len(hdu_r_stack)):
        # read WCS
        w = wcs.WCS(hdu_r_stack[jj].header)
        pixcrd = w.wcs_world2pix(cel_coord, 1)
        # if w.footprint_contains(sky_coord):
        if (pixcrd[0][0] > hthumb) & (pixcrd[0][0] < hdu_r_stack[jj].shape[0] - hthumb) & \
                (pixcrd[0][1] > hthumb) & (pixcrd[0][1] < hdu_r_stack[jj].shape[1] - hthumb):
            return hdu_r_stack[jj].name, pixcrd[0][0], pixcrd[0][1]

    return 0, 0, 0

for i in range(0, len(cat)):
# for i in range(0, 1000):
    if cat['col17'][i] != 1:       # only galaxies
        continue

    c = SkyCoord(f"{cat['col3'][i]}:{cat['col4'][i]}:{cat['col5'][i]} "
                       f"{cat['col6'][i]}:{cat['col7'][i]}:{cat['col8'][i]}", unit=(u.hourangle, u.deg))

    half_thumb = int(np.sqrt(cat['col13'][i])) * 4      # 4 x major axis as a box size

    tile, x, y = tile_check(c, half_thumb)

    if tile:
        #####  xxx.feedme    ########
        fm = open(f"/Users/duhokim/work/abell/galfit/A3558_no_weight/{cat['col1'][i]}.feedme",'w')
        fm.writelines(f"A) {ab.work_dir + f'fits/stacked/A3558_rsi.fits[{tile}]'} \n")		# Input data Image (FITS file)
        fm.writelines(f"B) {cat['col1'][i]}.out.fits \n")	# Output data image block
        # fm.writelines(f"C) {ab.work_dir + f'fits/stacked/A3558_rsw.fits[{tile}]'} \n")	# Sigma image name
        fm.writelines('C) none \n')  # Sigma image name
        fm.writelines('D) none \n')
        fm.writelines('E) 1 \n')
        fm.writelines('F) none \n')	# Bad pixel mask
        fm.writelines('G) none \n')
        fm.writelines(f"H) {int(x) - half_thumb} {int(x) + half_thumb} {int(y) - half_thumb} {int(y) + half_thumb} \n")	# Image region to fit (xmin, xmax, ymin, ymax)
        fm.writelines(f"I) {half_thumb} {half_thumb} \n")	# Size of the convolution box (x, y)
        fm.writelines('J) 31.456 \n')					# Magnitude photometric ZP
        fm.writelines('K) 0.2637 0.2637 \n')				# Plate scale (dx dy)	[arcsec per pixel]
        fm.writelines('O) regular \n')					# Display type (regular, curses, both)
        fm.writelines('P) 0 \n\n')					# Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps

        #######  Object number : 1  #######
        fm.writelines('0) expdisk \n')				# object type
        fm.writelines(f"1) {x} {y} 1 1 \n")	# position x, y
        fm.writelines(f"3) {cat['col11'][i]} 1 \n")	# Integrated magnitude
        fm.writelines(f"4) {int(cat['col13'][i]/10)} 1 \n")			# R_e (half-light radius)	[pix]
        fm.writelines('5) 1 1 \n')				# Sersic index n
        fm.writelines(f"9) {cat['col15'][i]/cat['col14'][i]} 1 \n") 			# axis ratio (b/a)
        fm.writelines(f"10) {cat['col16'][i]} 1 \n")				# position angle (PA) [deg: Up=0, Left=90]
        fm.writelines('Z) 0 \n\n')				# output option (0 = resig., 1 = Don't subtract)

        #######  Object number : 2  #######
        fm.writelines('0) devauc \n')                          # object type
        fm.writelines(f"1) {x} {y} 1 1 \n")	# position x, y
        fm.writelines(f"3) {cat['col11'][i]} 1 \n")	# Integrated magnitude
        fm.writelines(f"4) {int(cat['col13'][i]/10)} 1 \n")			# R_e (half-light radius)	[pix]
        fm.writelines('5) 1 1 \n')                              # Sersic index n
        fm.writelines(f"9) {cat['col15'][i]/cat['col14'][i]} 1 \n") 			# axis ratio (b/a)
        fm.writelines(f"10) {cat['col16'][i]} 1 \n")				# position angle (PA) [deg: Up=0, Left=90]
        fm.writelines('Z) 0 \n\n')                                # output option (0 = resig., 1 = Don't subtract)

        #######  Object number : 3  #######
#                fm.writelines('0) psf \n')                          # object type
#                fm.writelines('1) '+str(int(i))+' '+str(int(j))+' 1 1 \n')        # position x, y
#                fm.writelines('3) '+str(ngc_match.iloc[0,16]+3)+' 1 \n')    # Integrated magnitude
#                fm.writelines('Z) 0 \n\n')                                # output option (0 = resig., 1 = Don't subtract)

        #######  Object number : 4  #######
        fm.writelines('0) sky \n')                          	# object type
        fm.writelines('1) 1050 1 \n')        			# sky BG at center of fitting regions [ADUs]
        fm.writelines('2) 0.0 0 \n')    			# dsky/dx (sky gradient in x)
        fm.writelines('3) 0.0 0 \n')    			# dsky/dy (sky gradient in y)
        fm.writelines('Z) 0')                                # output option (0 = resig., 1 = Don't subtract)

        fm.close()
