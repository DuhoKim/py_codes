import numpy as np
from astropy.io import fits
import glob
import abell_cluster_module as ab
import os.path
from os import path

def nan_dqm(dqmArray, sciArray):
	sciData = np.array(sciArray, copy=True)
	indices = np.where(dqmArray != 0)
	sciData[indices] = np.nan
	return sciData

fits_dir = '/Users/duhokim/work/abell/fits/extracted/'

for k in range(3, len(ab.clusters)):
	for i in range(0, len(ab.bands)):
		hdu_dqm = fits.open(ab.work_dir + 'fits/stacked/' + ab.stack_dqm_fn[k][i])
		for j in range(1, len(hdu_dqm)):
			fn = f'{fits_dir}{ab.clusters[k]}_{ab.bands[i]}si_{j}.fits'
			if path.exists(fn):
				with fits.open(fn) as hdu:
					hdr = hdu[0].header
					data = hdu[0].data
					new = nan_dqm(hdu_dqm[j].data, data)		# zero out pixels lower than a fifth of background
					fits.writeto(fn[:-5]+'_nan.fits', new, hdr, overwrite=True)