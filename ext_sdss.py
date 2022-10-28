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

sex_dir = '/Users/duhokim/work/abell/fits/SDSS_extracted/'

cur_dir = os.getcwd()
os.chdir(sex_dir)

clusters = ['A2399', 'A2670']

for cluster in clusters:
	for band in ab.bands:
		fn = f'/Users/duhokim/work/abell/fits/SDSS_stacked/{cluster}_{band}si'
		for i in range(1, 26):	# for each tile
			os.system(f"sex {cluster}_{band}si_{i}.fits")
			os.system(f"psfex prepsfex.cat -PSF_SUFFIX .psf")
			os.system(f"mv prepsfex.psf {cluster}_{band}_{i}.psf")

os.chdir(cur_dir)