#
# Written by Duho Kim
# Korea Astronomy $ Space Science Institute
#
# You can freely use the code.
#
import numpy
from astropy.io import fits

def wei2rms(inputArray):
	"""Return 1/square root of the input numpy array.

	@type inputArray: numpy array
	@param inputArray: image data array
	@rtype: numpy array
	@return: image data array

	http://www.ifa.hawaii.edu/~rgal/science/sextractor_notes.html

	"""
	imageData = numpy.array(inputArray, copy=True)
	indices = numpy.where(imageData <= 0)
	imageData[indices] = 0.0
	indices = numpy.where(imageData > 0)
	imageData[indices] = 1/numpy.sqrt(imageData[indices])

	return imageData

fn = '/Users/duhokim/work/abell/fits/extracted/A3558_rsw_1.fits'
out_fn = '/Users/duhokim/work/abell/fits/extracted/A3558_rsr_1.fits'

if __name__ == "__main__":
	### READ FITS
	for i in range(1, 10):
		hdu = fits.open(f'/Users/duhokim/work/abell/fits/extracted/A754_rsw_{i}.fits')
		data = hdu[0].data
		rms = wei2rms(data)
		hdu_rms = fits.PrimaryHDU(rms)
		hdu_rms.writeto(f'/Users/duhokim/work/abell/fits/extracted/A754_rsr_{i}.fits')