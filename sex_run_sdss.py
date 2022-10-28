#####################################################
# Python script that preparing GALAPAGOS
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
# from multiprocessing import Pool

fits_dir = '/Users/duhokim/work/abell/fits/SDSS_extracted/'
sex_dir = '/Users/duhokim/work/abell/sex/run_SDSS/'

cur_dir = os.getcwd()
os.chdir(sex_dir)

NUM_CPUS = None	# defaults to all available

def sex_worker(cl, ban, ii):
	os.system(f"sex {fits_dir}{cl}_{ban}si_{ii}.fits -PSF_NAME {fits_dir}{cl}_{ban}_{ii}.psf "
			  f"-CATALOG_NAME {cl}_{ban}si_{ii}_psf.cat -MAG_ZEROPOINT 22.5 "
			  f"-CHECKIMAGE_TYPE SEGMENTATION -CHECKIMAGE_NAME {fits_dir}{cl}_{ban}si_{ii}_check.fits")

def test_run(pool):
	# for j in range(0, len(ab.clusters)):
	for j in range(1, 3):
		cluster = ab.clusters[j]
		for k in range(0, len(ab.bands)):
		# for k in range(1, 2):
			band = ab.bands[k]
			for iii in range(1, 26):
				pool.apply_async(sex_worker, args=(cluster, band, iii))

if __name__ == "__main__":
	import multiprocessing as mp
	pool = mp.Pool(NUM_CPUS)
	test_run(pool)
	pool.close()
	pool.join()

# for cluster in ab.clusters:			# for each cluster
# 	for i in range(1, 10):	# for each tile
# 		os.system(f"sex {fits_dir}{cluster}_rsi_{i}.fits -c prepsfex.sex -CATALOG_NAME "
# 				  f"{cluster}_rsi_{i}_FITS_LDAC.cat -MAG_ZEROPOINT 25.0 -VERBOSE_TYPE QUIET")
# 		os.system(f"psfex {cluster}_rsi_{i}_FITS_LDAC.cat")


os.chdir(cur_dir)