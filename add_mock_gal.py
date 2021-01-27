#####################################################
# Python script that adds mock galaxies on DECam mosaic
# written by Duho Kim (22 May 20)
######################################################
import numpy as np
from pyraf import iraf
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table
import os
import matplotlib.pyplot as plt
from scipy.stats import mode
import math
import shutil
import sys

work_dir = '/Users/duhokim/work/abell/'
orig = 'fits/best_single/A2399_ri_300_0142.fits'
seg = 'sex/seg/A2399_rps_300_0142.check.fits'

##### READ FITS #####
hdu_orig = fits.open(work_dir+orig)
hdu_seg = fits.open(work_dir+seg)

for i in range(1, 61):
	fits_orig = hdu_orig[i].data
	fits_seg = hdu_seg[i].data

	# Adopt median of background(where SEGMENTATION=0) as background value (mode function takes too much time?)
	bg = np.median(fits_orig[np.where(fits_seg==0)])
	sig = np.std(fits_orig[np.where(fits_seg==0)])			# Adopt std of background image as sigma value

	print('median background value is %d' % bg)
	print('standard deviation of the background value is %d' % sig)
