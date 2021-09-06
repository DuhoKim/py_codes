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

sex_dir = '/Users/duhokim/work/abell/fits/extracted/'

cur_dir = os.getcwd()
os.chdir(sex_dir)

cluster = ab.clusters[0]
bands = ['u', 'g']

for j in range(6, 7):			# for each cluster
	cluster = ab.clusters[j]
	for k in range(0, 2):		# for each band
		band = bands[k]
		for i in range(1, 10):	# for each tile
			ext(f'/Users/duhokim/work/abell/fits/stacked/{cluster}_{band}si.fits',
				f'{sex_dir}{cluster}_{band}si_{i}.fits',
				i)

			ext(f'/Users/duhokim/work/abell/fits/stacked/{cluster}_{band}sw.fits',
				f'{sex_dir}{cluster}_{band}sw_{i}.fits',
				i)

			hdu = fits.open(f'{sex_dir}{cluster}_{band}sw_{i}.fits')
			data = hdu[0].data
			rms = wei.wei2rms(data)
			hdu_rms = fits.PrimaryHDU(rms)
			hdu_rms.writeto(f'{sex_dir}{cluster}_{band}sr_{i}.fits')

			os.system(f"sex {cluster}_{band}si_{i}.fits")
			os.system(f"psfex prepsfex.cat")
			hdu = fits.open('prepsfex.fits')
			data = hdu[1].data[0][0][0]
			header = hdu[0].header
			fits.writeto(f'{cluster}_{band}si_{i}_psf.fits', data, header)

os.chdir(cur_dir)