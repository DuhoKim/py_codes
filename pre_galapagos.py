#####################################################
# Python script that preparing GALAPAGOS
# written by Duho Kim (10 Aug 21)
######################################################
import os
from shutil import copyfile
from astropy.io import fits

sex_dir = '/Users/duhokim/work/abell/fits/extracted/'

cur_dir = os.getcwd()
os.chdir(sex_dir)

for i in range(2, 10):
	os.system(f"sex A754_rsi_{i}.fits")
	os.system(f"psfex prepsfex.cat")
	hdu = fits.open('prepsfex.fits')
	data = hdu[1].data[0][0][0]
	header = hdu[0].header
	fits.writeto(f'A754_rsi_{i}_psf.fits', data, header)

os.chdir(cur_dir)