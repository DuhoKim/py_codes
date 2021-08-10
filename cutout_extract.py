#
# Written by Duho Kim
# Korea Astronomy $ Space Science Institute
#
# You can freely use the code.
#
import numpy
from astropy.io import fits

fn1 = '/Users/duhokim/work/abell/fits/extracted/A3558_rsi_1.fits'
fn2 = '/Users/duhokim/work/abell/fits/extracted/A3558_rsr_1.fits'
fn3 = '/Users/duhokim/work/abell/fits/extracted/A3558_rsw_1.fits'

out_fn1 = '/Users/duhokim/work/abell/fits/extracted/A3558_rsi_1_ADU.fits'
out_fn2 = '/Users/duhokim/work/abell/fits/extracted/A3558_rsr_1_ADU.fits'
out_fn3 = '/Users/duhokim/work/abell/fits/extracted/A3558_rsw_1_ADU.fits'

size = 500

### READ FITS
for i in range(3, 10):
    hdu = fits.open(fn1[:-6]+f'{i}'+'.fits')
    data = hdu[0].data
    tiny = data * 10800.0
    hdu_tiny = fits.PrimaryHDU(tiny)
    hdu_tiny.header = hdu[0].header
    hdu_tiny.writeto(out_fn1[:-10]+f'{i}'+'_ADU.fits', overwrite=True)

# hdu = fits.open(fn2)
# data = hdu[0].data
# tiny = data[790-size:790+size, 10815-size:10815+size]
# hdu_tiny = fits.PrimaryHDU(tiny)
# hdu_tiny.header = hdu[0].header
# hdu_tiny.writeto(out_fn2, overwrite=True)
#
# hdu = fits.open(fn3)
# data = hdu[0].data
# tiny = data[790-size:790+size, 10815-size:10815+size]
# hdu_tiny = fits.PrimaryHDU(tiny)
# hdu_tiny.header = hdu[0].header
# hdu_tiny.writeto(out_fn3, overwrite=True)
