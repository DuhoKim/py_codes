#####################################################
# Python script that generate SEx scripts
# written by Duho Kim (21 Apr 20)
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
import in_place

sex_dir = '/Users/dkim108/Documents/work/sex/'

dates = ['19', '20', '21']

for i in range(0, len(dates)):
	data_dir = '../../ssd_dhkim/DECam/'+dates[i]+'aug/'
	f_sex = open(dates[i]+'aug_sex.sh', 'w')
	f_sex.write('#!/bin/bash\n\n')
    with in_place.InPlace(sex_dir+dates[i]+'.list') as f:
        for line in f:
			if 'i_' in line:
				weight = line.replace('i_', 'w_')
				f_sex.write('sextractor '+data_dir+line+' -WEIGHT_IMAGE '+data_dir+
							weight+' -WEIGHT_TYPE MAP_WEIGHT -CATALOG_NAME '+dates[i])
            line = line.replace('*', '')
            f.write(line)

am_cat = pd.read('/Users/dkim108/Documents/work/cat/airmass.csv')

work_dir = '/Users/dkim108/Documents/work/DECam/ss/'
dates = ['19aug', '20aug', '21aug']

os.chdir(work_dir)

#### BASIC PARAMETERS ##########
pixel_size = 0.2631      # Pixel size of DECam data
decam_fwhm = 1.0         # Initial guess of FWHM of DECam data ["]
zmag = 25.0		# Zero point of magnitude scale

#num_psf_stars = 10 (SHA)
# num_psf_stars = 15      # Number of PSF stars
# psf_radius = 30         # Radius of PSF [pixel]
# nan_check = 2		# NaN check around a PSF stars with a radius of nan_check * psf_radius
tol_star_match = 1.0    # Tolerance of coordinate matching of PSF stars between sex catalog and IRAF tables

##### SELECT PSF STARS #########################################
# Select objects from the catalog resulted by Source Extractor
# above which is having 'FLAGS' 0 and 'CLASS_STAR' larger than 
# half and ELLIPTICITY (1-b/a) is less than 0.1 and FWHM_IMAGE
# is in between +- 1.5 pixels centered in mode value of FWHM_IMAGE
# values and sort with magnitudes 'MAG_AUTO' to build PSF from
################################################################
#class_range = 0.7 (SHA)
class_range = 0.8       # Maximum value of CLASS_STAR selecting PSF stars form sex catalog
ellip_range = 0.1       # Maximum value of ELLIPTICITY selecting PSF stars from sex catalog 

##### PARAMETER SETTING FOR PyRAF ##############################################
# referenced Manuel 'PSF photometry using DAOPHOT' written by H. Hwang & M. Lee
################################################################################
iraf.daophot()
iraf.daopars.function = 'auto'          # fit all 6 analytic functions and use a function having minimum norm scatter value
iraf.daopars.varorder = 0               # -1 results simple not-realistic model, 1 and 2 are noisy
iraf.daopars.nclean = 5                 # showed improvement with an increament until 10 but not much difference compared to 100
                                        # The number of additional iterations the PSF task performs to compute the PSF look-up tables.
                                        # If nclean is > 0, stars which contribute deviant residuals to the PSF look-up tables in the 
                                        # first iteration, will be down-weighted in succeeding iterations. (default=0)
iraf.daopars.matchrad = 3.0             # Object matching radius in scale units (default=3.0)
iraf.daopars.psfrad = 93        # Opt to use 60x60 size PSF image size (DK)

iraf.daopars.recenter = 'yes'           # Recenter stars during fit?
iraf.daopars.fitsky = 'yes'             # Recompute group sky value during fit? 
iraf.daopars.groupsky = 'yes'           # Use group rather than individual sky values? (default=yes)
iraf.daopars.sannulus = 0.0             # The inner radius of the sky annulus used by ALLSTAR to recompute the sky values. (default=0.0)
iraf.daopars.wsannulus = 11.0           # Width of sky fitting annulus in scale units (default=11.0)
iraf.daopars.critsnratio = 1.0          # The ratio of the model intensity of the brighter star computed at a distance of one fitting radius from the center
                                        # of the fainter star,  (default=1.0)
#iraf.findpars.threshold = 4.0           # Threshold in sigma for feature detection (default=4.0)
iraf.findpars.threshold = 50.0           # Threshold in sigma for feature detection (default=4.0)
iraf.findpars.nsigma = 1.5              # Width of convolution kernel in sigma  (default=1.5)
iraf.findpars.sharplo = 0.2             # Lower bound on sharpness for feature detection (default=0.2)
iraf.findpars.sharphi = 1.0             # Upper bound on sharpness for feature detection (default=1.0)

for k in range(0, len(dates)):
	os.chdir(work_dir+dates[k])
	ss_imgs = os.listdir(work_dir+dates[k])

	for x in range(0, len(ss_imgs)):
		if '.fits' in ss_imgs[x] and 'w' not in ss_imgs[x]:

			fn = ss_imgs[x]		# Image file name

			##### READ FITS #####
			hdu = fits.open(fn)
			for i in range(1, 61):
				fits_data = hdu[i].data

				bg = np.median(fits_data)			# Adopt median of background image as background value (mode function takes too much time)
				sig = np.std(fits_data)			# Adopt std of background image as sigma value

				fwhm = 4.2
				##### PARAMETER SETTING FOR PyRAF ##############################################
				# referenced Manuel 'PSF photometry using DAOPHOT' written by H. Hwang & M. Lee
				################################################################################
				iraf.datapars.fwhmpsf = fwhm
				iraf.datapars.sigma = sig

				iraf.centerpars.cbox = max(5, 2*fwhm)
				iraf.fitskypars.skyvalu = bg
				iraf.fitskypars.annulus = 4*fwhm
				iraf.fitskypars.dannulus = 3*fwhm
				#iraf.photpars.aperture = max(3, fwhm)
				iraf.photpars.aperture = max(3, 91)		# from 24" diameter from Smith+06
				iraf.photpars.zmag = zmag
				iraf.daopars.fitrad = fwhm              # Fitting radius in scale units. Only pixels within the fitting radius of the center of a star will contribute to
													# the fits computed by the PEAK, NSTAR and ALLSTAR tasks.

				##### Get Coordiates and Magnitudes Stars #################
				# Run DAOFIND to get coordinates and PHOT to get magnitudes
				###########################################################
				iraf.daofind(fn+'['+str(i)+']', '../'+dates[k]+'_phot/'+fn[:-5]+'.coo.'+str(i),verify='no',verbose='no')			# generate coo.1 file
				iraf.phot(fn+'['+str(i)+']', '../'+dates[k]+'_phot/'+fn[:-5]+'.coo.'+str(i),'../'+dates[k]+'_phot/'+fn[:-5]+'.mag.'+str(i),verify='no',verbose='no')	# generate mag.1 file

				##### MAKE PSF STAR LIST FILE 'pst.1' #######
				# read coordinates and magnitudes of the stars selected based
				# on the catalog generated by Source Extractor and make
				# 'pst.1' file to use for running PSF task
				#######################################################
				# coo1 = ascii.read(fn[:-5]+'.coo.'+str(i))
				# mag1 = ascii.read(fn[:-5]+'.mag.'+str(i))
				# coo1.sort('ID')
				# mag1.sort('ID')
				# pst1 = Table(names=('ID','XCENTER','YCENTER','MAG','MSKY'),dtype=('i4','f8','f8','f8','f8'))
				# for j in range(min(len(mag1),len(coo1))):
				# 	ind = np.where((abs(mag1['XCENTER'][j]-coo1['XCENTER']) < tol_star_match) & \
				# 		(abs(mag1['Ycenter'][j]-coo1['YCENTER']) < tol_star_match))[0]
				# 	pst1.add_row( ( mag1['ID'][ind[0]], mag1['XINIT'][ind[0]], mag1['YINIT'][ind[0]], mag1['MAG'][ind[0]], mag1['MSKY'][ind[0]]) )
				#
				# pst1.write(fn[:-5]+'.pst.'+str(i), format='ascii.commented_header', delimiter='\t', comment= '#N ID    XCENTER   YCENTER \
				# 	MAG         MSKY\n#U ##    pixels    pixels    magnitudes  counts\n#F %-9d  %-10.3f   %-10.3f   %-12.3f     %-15.7g\n#')










