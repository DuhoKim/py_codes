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

fits_dir = '/Users/duhokim/work/abell/fits/best_single_extracted/'
mock_dir = '/Users/duhokim/work/abell/fits/best_single_extracted_mock/'
sex_dir = '/Users/duhokim/work/abell/sex/run_short_mock/'

cur_dir = os.getcwd()
os.chdir(sex_dir)

bands = ['u', 'g', 'r']

NUM_CPUS = None	# defaults to all available

def sex_worker(cl, ban, ii):
	fn = f'{cl}_{ban}si_{ii}_1src_at_50_50'
	os.system(f"sex {mock_dir}{fn}.fits -PSF_NAME {fits_dir}{cl}_{ban}_{ii}.psf "
			  f"-CATALOG_NAME {fn}.cat -MAG_ZEROPOINT 25.0 -VERBOSE_TYPE QUIET ")

def test_run(pool):
	for j in range(0, 1):
		cluster = ab.clusters[j]
		for k in range(0, 1):
			band = bands[k]
			fn = f'/Users/duhokim/work/abell/fits/best_single/{ab.short_sci_fn[j][k]}'
			hdu = fits.open(fn)
			# for iii in range(1, len(hdu)):		# for each tile
			for iii in range(1, 2):  # for each tile
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