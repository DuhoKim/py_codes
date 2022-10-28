#####################################################
# Python script that extract best single
# written by Duho Kim (10 Aug 21)
######################################################
import os
from shutil import copyfile
from astropy.io import fits
from splitfits import extract_extension as ext
import abell_cluster_module as ab
import weight2rms as wei
import importlib
importlib.reload(wei)
#from weight2rms import wei2rms

fits_dir = '/Users/duhokim/work/abell/fits/DECaLS/'

cur_dir = os.getcwd()
os.chdir(fits_dir)

clusters = ['A2399', 'A2670', 'A3716']
bands = ['g', 'r']

for cluster in clusters:			# for each cluster
	for band in bands:
		for k in range(1, 26):		# for each tile
			fn = f'{fits_dir}{cluster}_{k}.fits'

			os.system(f"sex {cluster}_{band}si_{k}.fits")
			os.system(f"psfex prepsfex.cat -PSF_SUFFIX .psf")
			os.system(f"mv prepsfex.psf {cluster}_{band}_{k}.psf")

os.chdir(cur_dir)