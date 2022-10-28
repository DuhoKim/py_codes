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

sex_dir = '/Users/duhokim/work/abell/fits/best_single_extracted/'

cur_dir = os.getcwd()
os.chdir(sex_dir)

cluster = ab.clusters[0]
bands = ['u', 'g', 'r']

for j in range(2, 7):			# for each cluste
	cluster = ab.clusters[j]
	for k in range(0, 1):		# for each band
		band = bands[k]
		fn = f'/Users/duhokim/work/abell/fits/best_single/{ab.short_sci_fn[j][k]}'
		hdu = fits.open(fn)
		for i in range(1, len(hdu)):	# for each tile
			ext(fn, f'{sex_dir}{cluster}_{band}si_{i}.fits', i)
			os.system(f"sex {cluster}_{band}si_{i}.fits")
			os.system(f"psfex prepsfex.cat")
			os.system(f"mv prepsfex.psf {cluster}_{band}_{i}.psf")

os.chdir(cur_dir)