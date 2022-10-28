#####################################################
# Python script that adds mock galaxies on DECam mosaic
# written by Duho Kim (22 May 20)
######################################################
import numpy as np
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table
import os
import matplotlib.pyplot as plt
from scipy.stats import mode
import math
import shutil
import sys
from astropy.modeling.models import Sersic2D, Gaussian2D

work_dir = '/Users/duhokim/work/abell/'

with fits.open(work_dir+'fits/best_single_extracted/A754_usi_1.fits') as hdu_orig,	\
	fits.open(work_dir+'fits/best_single_extracted/A754_usi_1_check.fits') as hdu_ch:

	img = hdu_orig[0].data
	x, y = np.meshgrid(np.arange(hdu_orig[0].shape[1]), np.arange(hdu_orig[0].shape[0]))
	mod = Gaussian2D(amplitude=1000, x_mean=50, y_mean=50, x_stddev=2, y_stddev=2)
	img_mod = mod(x, y)
	new_img = img + img_mod

	fits.writeto(f'{work_dir}fits/best_single_extracted_mock/A754_usi_1_1src_at_50_50.fits',
				 new_img,
				 hdu_orig[0].header)

# mod = Sersic2D(amplitude = 1, r_eff = 25, n=4, x_0=50, y_0=50,
# 			   ellip=.5, theta=-1)

# img = mod(x, y)
	log_img = np.log10(img)
	log_mod = np.log10(img_mod)
	log_mock = np.log10(img + img_mod)
#
# psf = fits.open(work_dir+'fits/best_single_extracted/A3716_u_59.psf')
# psf_img = psf[1].data[0][0][0]

	plt.figure()
	plt.imshow(log_mock, origin='lower', interpolation='nearest',
			   vmin=1, vmax=2)
	plt.xlabel('x')
	plt.ylabel('y')
	cbar = plt.colorbar()
	cbar.set_label('Log Brightness', rotation=270, labelpad=25)
	cbar.set_ticks([-1, 0, 1, 2], update_ticks=True)
	plt.show()


##### READ FITS #####
# hdu_orig = fits.open(work_dir+orig)
# hdu_seg = fits.open(work_dir+seg)

# for i in range(1, 61):
# 	fits_orig = hdu_orig[i].data
# 	fits_seg = hdu_seg[i].data
#
# 	# Adopt median of background(where SEGMENTATION=0) as background value (mode function takes too much time?)
# 	bg = np.median(fits_orig[np.where(fits_seg==0)])
# 	sig = np.std(fits_orig[np.where(fits_seg==0)])			# Adopt std of background image as sigma value
#
# 	print('median background value is %d' % bg)
# 	print('standard deviation of the background value is %d' % sig)
