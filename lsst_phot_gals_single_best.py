from astropy.io import ascii
import numpy as np
import time
from astropy.coordinates import SkyCoord
import my_module
import pandas as pd
from astropy.io import fits
from astropy import wcs
from datetime import datetime
from os import listdir
from os.path import isfile, join
import math

fits_dir=("/Users/duhokim/work/abell/")
work_dir=("/Users/duhokim/lsst_stack/demo_data_orig_20/DATA/rerun/processCcdOutputs/")
cat_dir=("/Users/duhokim/work/abell/cat/")

visit_nums = ['0350110', '0350835']

ver = 'v1.1'
class_star_lim = 0.9
fwhm_lim = 0.0
mag_lim = [99, 99, 99]

max_sep = 1.0   # matching radius limit among each individual exposures
mag_sys = 'MAG_ISO'
magerr_sys = 'MAGERR_ISO'
mag_sys1 = 'MAG_AUTO'
magerr_sys1 = 'MAGERR_AUTO'
ccd_x = 2046
ccd_y = 4094

pix_scale = 0.263

clusters = ['A2670']
single_fn = ['A2670_gs_60_0420.cat', 'A2670_rs_300_0351.cat']
single_dqm_fn = ['A2670_gd_60_0420_21aug.fits', 'A2670_rd_300_0351_19aug.fits']

mag_exp_t_corr_300 = 2.5 * np.log10(300)
mag_exp_t_corr_60 = 2.5 * np.log10(60)

### Standardization parameters      mag = inst_mag [ZP:25 with flux/s] + a + b * AIRMASS
single_a = [0.179 + mag_exp_t_corr_60, 0.567 + mag_exp_t_corr_300]

single_b = [-0.036 * 1.26, -0.136 * 1.41]

# irsa.ipac.caltech.edu, S and F (2011)
mw_ext = [0.146, 0.101]

# am_cat = pd.read_csv(work_dir+'sex/cat/airmass.csv')

start_time = time.time()

fn = cat_dir+'lsst_19_21_aug_2014_single_best_exposure_SEx_cat_' + clusters[0] + \
     '_match_rad_1as_Gal_ext_corrected_'+ver

with open(fn + '.txt', 'w') as fh, open(fn + '.reg', 'w') as regFile:
    # read Data Quality Mask (DQM) fits file to check saturation
    hdu_r = fits.open(fits_dir + 'fits/best_single/' + single_dqm_fn[1])

    mypath = work_dir + visit_nums[0] + '/src/'
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

    id_num = 1
    for i in range(0, len(onlyfiles)):
        if int(onlyfiles[i][-7:-5]) == 62 or int(onlyfiles[i][-7:-5]) == 1:
            continue
        with fits.open(join(mypath, onlyfiles[i])) as hdu:
            cat = hdu[1].data
            w = wcs.WCS(hdu_r[int(onlyfiles[i][-7:-5])-1].header)
            for j in range(0, len(cat)):
                if np.isnan(cat['base_NaiveCentroid_x'][j]) or cat['base_PsfFlux_instFlux'][j] < 1e4:       # only bright sources
                    continue
                pixcrd = [[cat['base_NaiveCentroid_x'][j], cat['base_NaiveCentroid_y'][j]], [0, 0]]
                world = w.wcs_pix2world(pixcrd, 0)
                fh.writelines("{:4}".format(id_num) + ' ' +  # NUMBER
                              "{:12.7f}".format(world[0][0]) + ' ' +  # RA
                              "{:12.7f}".format(world[0][1]) + ' ' +  # DEC
                              "{:7.3f}".format(cat['base_PsfFlux_instFlux'][j]) + '\n'
                              )

                xx = cat['base_SdssShape_xx'][j]
                yy = cat['base_SdssShape_yy'][j]
                xy = cat['base_SdssShape_xy'][j]
                theta = 0.5 * math.atan2(2 * xy, xx - yy) * 180 / np.pi
                if np.isnan(theta):       # only bright sources
                    continue
                regFile.writelines("j2000; ellipse({:12.7f}, {:12.7f}, {:7.3f}\", {:7.3f}\", {:7.3f}) # text=\'{}\', "
                                   "color={}, dash=1 \n".format(
                    world[0][0],
                    world[0][1],
                    np.sqrt(xx) * pix_scale,
                    np.sqrt(yy) * pix_scale,
                    theta,
                    id_num,
                    'blue'))
                id_num = id_num + 1


print("--- %s minutes ---" % (((time.time() - start_time))/60.0))

my_module.print_time()