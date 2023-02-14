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

# work_dir = ("/Users/user1/Library/CloudStorage/OneDrive-충남대학교/work/abell/")

sex_dir = '/Users/user1/Library/CloudStorage/OneDrive-충남대학교/work/abell/fits/extracted/'

cur_dir = os.getcwd()
os.chdir(sex_dir)

bands = ['r']

for cluster in ab.clusters:
	for band in bands:
		fn = f'/Users/user1/Library/CloudStorage/OneDrive-충남대학교/work/abell/fits/stacked/{cluster}_{band}si'
		hdu = fits.open(fn+'.fits')
		for i in range(1, len(hdu)):	# for each tile
			os.system(f"sex {cluster}_{band}si_{i}.fits")
			os.system(f"psfex prepsfex.cat -PSF_SUFFIX .fits")
			os.system(f"mv prepsfex.fits {cluster}_{band}_{i}.fits")

os.chdir(cur_dir)