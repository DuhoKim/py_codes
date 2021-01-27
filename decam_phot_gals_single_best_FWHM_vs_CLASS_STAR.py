from astropy.io import ascii
import numpy as np
import time
from astropy.coordinates import SkyCoord
import my_module
import pandas as pd
import matplotlib.pyplot as plt

plot_dir=("/Users/duhokim/work/abell/plot/")
sex_dir=("/Users/duhokim/work/abell/sex/cat/")

pix_scale = 0.263
class_star_lim = 0.5
max_sep = 1.0   # matching radius limit among each individual exposures
flag_lim = 0    # no neighbor no blend no saturate
mag_sys = 'MAG_ISO'
magerr_sys = 'MAGERR_ISO'

clusters = ['A2399', 'A2670', 'A3716']
# single_fn = [['A2399_us_300_0236_20aug.cat', 'A2399_gs_60_0239_19aug.cat', 'A2399_rs_300_0142_19aug.cat'],
#              ['A2670_us_300_0409_21aug.cat', 'A2670_gs_60_0420_21aug.cat', 'A2670_rs_300_0351_19aug.cat'],
#              ['A3716_us_300_2357_21aug.cat', 'A3716_gs_300_0030_21aug.cat', 'A3716_rs_300_0418_19aug.cat']]
single_fn = [['A2399_us_300_0236.cat', 'A2399_gs_60_0239.cat', 'A2399_rs_300_0142.cat'],
             ['A2670_us_300_0409.cat', 'A2670_gs_60_0420.cat', 'A2670_rs_300_0351.cat'],
             ['A3716_us_300_2357.cat', 'A3716_gs_300_0030.cat', 'A3716_rs_300_0418.cat']]

mag_exp_t_corr_300 = 2.5 * np.log10(300)
mag_exp_t_corr_60 = 2.5 * np.log10(60)

### Standardization parameters      mag = inst_mag [ZP:25 with flux/s] + a + b * AIRMASS
single_a = [[-1.257 + mag_exp_t_corr_300, 0.381 + mag_exp_t_corr_60, 0.567 + mag_exp_t_corr_300],
            [-1.796 + mag_exp_t_corr_300, 0.179 + mag_exp_t_corr_60, 0.567 + mag_exp_t_corr_300],
            [-1.796 + mag_exp_t_corr_300, 0.179 + mag_exp_t_corr_300, 0.567 + mag_exp_t_corr_300]]

single_b = [[-0.592 * 1.26, -0.186 * 1.26, -0.136 * 1.53],
            [-0.132 * 1.3, -0.036 * 1.26, -0.136 * 1.41],
            [-0.132 * 1.42, -0.036 * 1.31, -0.136 * 1.09]]

# irsa.ipac.caltech.edu, S and F (2011)
mw_ext = [[0.159, 0.124, 0.086],
               [0.188, 0.146, 0.101],
               [0.157, 0.122, 0.085]]

am_cat = pd.read_csv('/Users/duhokim/work/abell/sex/cat/airmass.csv')

fig = plt.figure()

for k in range(0, len(clusters)):

    # read SEx cat
    sex_cat_u = ascii.read(sex_dir + 'best_single/' + single_fn[k][0])
    sex_cat_g = ascii.read(sex_dir + 'best_single/' + single_fn[k][1])
    sex_cat_r = ascii.read(sex_dir + 'best_single/' + single_fn[k][2])

    plt.scatter(sex_cat_u['CLASS_STAR'],
                      sex_cat_u['FWHM_IMAGE'] * pix_scale,
                      label=clusters[k]+' u',
                      alpha=0.1,
                      s=1)
    plt.scatter(sex_cat_g['CLASS_STAR'],
                      sex_cat_g['FWHM_IMAGE'] * pix_scale,
                      label=clusters[k]+' g',
                      alpha=0.1,
                      s=1)
    plt.scatter(sex_cat_r['CLASS_STAR'],
                      sex_cat_r['FWHM_IMAGE'] * pix_scale,
                      label=clusters[k]+' r',
                      alpha=0.1,
                      s=1)

plt.ylim([0, 4])
plt.ylabel('FWHM in arcsec')
plt.legend(loc='upper right')
fig.savefig(plot_dir + 'DECam_FWHM_vs_CLASS_STAR.png')
