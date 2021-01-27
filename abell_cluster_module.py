import numpy as np
from astropy.coordinates import SkyCoord

ver = '2013-2014_bright'
class_star_lim = 0.5
max_sep = 1.0   # matching radius limit among each individual exposures
mag_sys = 'MAG_AUTO'
magerr_sys = 'MAGERR_AUTO'
mag_lim = 22
rad_lim = 0.8 # only include detection inside [degree]

work_dir = ("/Users/duhokim/work/abell/")
plot_dir = work_dir+'plot/'
sex_dir = work_dir+'sex/cat/'
cat_dir = work_dir+'cat/'
short_cat_dir=("/Volumes/APPLE SSD/Server/work/sex/best_single/")

clusters = ['A754', 'A2399', 'A2670', 'A3558', 'A3574', 'A3659', 'A3716']
bands = ['u', 'g', 'r']

ccd_x = 2046
ccd_y = 4094

mag_exp_t_corr_300 = 2.5 * np.log10(300)
mag_exp_t_corr_60 = 2.5 * np.log10(60)

short_sci_fn = [['A754_ui_300_11apr.fits',          'A754_gi_300_11apr.fits',       'A754_ri_300_10apr.fits'],
                ['A2399_ui_300_0236_20aug.fits',    'A2399_gi_60_0239_19aug.fits',  'A2399_ri_300_0142_19aug.fits'],
                ['A2670_ui_300_0409_21aug.fits',    'A2670_gi_60_0420_21aug.fits',  'A2670_ri_300_0351_19aug.fits'],
                ['A3558_ui_300_11apr.fits',         'A3558_gi_300_11apr.fits',      'A3558_ri_300_10apr.fits'],
                ['A3574_ui_300_11apr.fits',         'A3574_gi_300_11apr.fits',      'A3574_ri_300_10apr.fits'],
                ['A3659_ui_300_11apr.fits',         'A3659_gi_300_11apr.fits',      'A3659_ri_300_11apr.fits'],
                ['A3716_ui_300_2357_21aug.fits',    'A3716_gi_300_0030_21aug.fits', 'A3716_ri_300_0418_19aug.fits']]

short_dqm_fn = [['A754_ud_300_11apr.fits', 'A754_gd_300_11apr.fits', 'A754_rd_300_10apr.fits'],
                ['A2399_ud_300_0236_20aug.fits', 'A2399_gd_60_0239_19aug.fits', 'A2399_rd_300_0142_19aug.fits'],
                ['A2670_ud_300_0409_21aug.fits', 'A2670_gd_60_0420_21aug.fits', 'A2670_rd_300_0351_19aug.fits'],
                ['A3558_ud_300_11apr.fits', 'A3558_gd_300_11apr.fits', 'A3558_rd_300_10apr.fits'],
                ['A3574_ud_300_11apr.fits', 'A3574_gd_300_11apr.fits', 'A3574_rd_300_10apr.fits'],
                ['A3659_ud_300_11apr.fits', 'A3659_gd_300_11apr.fits', 'A3659_rd_300_11apr.fits'],
                ['A3716_ud_300_2357_21aug.fits', 'A3716_gd_300_0030_21aug.fits', 'A3716_rd_300_0418_19aug.fits']]

short_cat_fn = [['A754_us_300',         'A754_gs_300',      'A754_rs_300'],
                ['A2399_us_300_0236',   'A2399_gs_60_0239', 'A2399_rs_300_0142'],
                ['A2670_us_300_0409',   'A2670_gs_60_0420', 'A2670_rs_300_0351'],
                ['A3558_us_300',        'A3558_gs_300',     'A3558_rs_300'],
                ['A3574_us_300',        'A3574_gs_300',     'A3574_rs_300'],
                ['A3659_us_300',        'A3659_gs_300',     'A3659_rs_300'],
                ['A3716_us_300_2357',   'A3716_gs_300_0030', 'A3716_rs_300_0418']]

coords_cl_cen = [SkyCoord(137.2850624, -9.7128782, unit='deg'),
                 SkyCoord(329.3575, -7.794722, unit='deg'),
                 SkyCoord(358.5291618, -10.3350941, unit='deg'),
                 SkyCoord(202.6336399, -31.6472662, unit='deg'),
                 SkyCoord(207.0703863, -30.3243432, unit='deg'),
                 SkyCoord(300.6202774, -30.0870595, unit='deg'),
                 SkyCoord(312.819614, -52.695408, unit='deg')]

short_a = [[-1.534 + mag_exp_t_corr_300, 0.426 + mag_exp_t_corr_300, 0.498 + mag_exp_t_corr_300],
            [-1.6863 + mag_exp_t_corr_300, 0.3501 + mag_exp_t_corr_60, 0.4448 + mag_exp_t_corr_300],
            [-1.7019 + mag_exp_t_corr_300, 0.2680 + mag_exp_t_corr_60, 0.4448 + mag_exp_t_corr_300],
            [-1.534 + mag_exp_t_corr_300, 0.426 + mag_exp_t_corr_300, 0.498 + mag_exp_t_corr_300],
            [-1.534 + mag_exp_t_corr_300, 0.426 + mag_exp_t_corr_300, 0.498 + mag_exp_t_corr_300],
            [-1.534 + mag_exp_t_corr_300, 0.426 + mag_exp_t_corr_300, 0.513 + mag_exp_t_corr_300],
            [-1.7019 + mag_exp_t_corr_300, 0.2680 + mag_exp_t_corr_300, 0.4448 + mag_exp_t_corr_300]]

short_b = [[-0.355 * 1.19,    -0.166 * 1.07, -0.061 * 1.07],
            [-0.1724 * 1.26,    -0.1631 * 1.26, -0.0544 * 1.53],
            [-0.2034 * 1.3,     -0.0846 * 1.26, -0.0544 * 1.41],
            [-0.355 * 1.04,     -0.166 * 1.0, -0.061 * 1.0],
            [-0.355 * 1.32,     -0.166 * 1.13, -0.061 * 1.23],
            [-0.355 * 1.1,     -0.166 * 1.2, -0.073 * 1.09],
            [-0.2034 * 1.42,    -0.0846 * 1.31, -0.0544 * 1.09]]

# # statsmodel + default SEx param
stack_a = [
    [3.962, 6.148, 6.458],
    [4.186, 6.1368, 6.4124],
    [3.9915, 6.1074, 6.2289],
    [4.020, 6.176, 6.456],
    [3.938, 6.204, 6.453],
    [3.754, 5.776, 6.167],
    [4.0240, 6.2455, 6.4456]
]

