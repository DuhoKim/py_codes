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

fits_dir = '/Users/duhokim/work/abell/fits/extracted/'
sex_dir = '/Users/duhokim/work/abell/sex/run/'

cur_dir = os.getcwd()
os.chdir(sex_dir)

NUM_CPUS = None	# defaults to all available

def sex_worker(cl, ii):
	os.system(f"sex {fits_dir}{cl}_rsi_{ii}.fits -PSF_NAME {cl}_rsi_{ii}_FITS_LDAC.psf "
			  f"-CATALOG_NAME {cl}_rsi_{ii}_psf.cat -MAG_ZEROPOINT 25.0 -VERBOSE_TYPE QUIET")

def test_run(pool):
	for cluster in ab.clusters:
		for iii in range(1, 10):
			pool.apply_async(sex_worker, args=(cluster, iii))

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