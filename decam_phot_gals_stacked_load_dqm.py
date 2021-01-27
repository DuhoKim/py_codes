from astropy.io import ascii
from astropy.table import Table, vstack
import numpy as np
import time
from astropy.coordinates import SkyCoord
import my_module
import pandas as pd
from astropy.io import fits
from astropy import wcs
from datetime import datetime
from astropy import units as u
from multiprocessing import Process, Pool, Lock

work_dir=("/Users/duhokim/work/abell/")

ver = 'v1.2'
class_star_lim = 0.9
max_sep = 1.0   # matching radius limit among each individual exposures
mag_sys = 'MAG_ISO'
magerr_sys = 'MAGERR_ISO'
mag_sys1 = 'MAG_AUTO'
magerr_sys1 = 'MAGERR_AUTO'

clusters = ['A2399', 'A2670', 'A3716']

# single_fn = [['A2399_us_300_0236_20aug.cat', 'A2399_gs_60_0239_19aug.cat', 'A2399_rs_300_0142_19aug.cat'],
#              ['A2670_us_300_0409_21aug.cat', 'A2670_gs_60_0420_21aug.cat', 'A2670_rs_300_0351_19aug.cat'],
#              ['A3716_us_300_2357_21aug.cat', 'A3716_gs_300_0030_21aug.cat', 'A3716_rs_300_0418_19aug.cat']]
single_fn = [['A2399_us_300_0236.cat', 'A2399_gs_60_0239.cat', 'A2399_rs_300_0142.cat'],
             ['A2670_us_300_0409.cat', 'A2670_gs_60_0420.cat', 'A2670_rs_300_0351.cat'],
             ['A3716_us_300_2357.cat', 'A3716_gs_300_0030.cat', 'A3716_rs_300_0418.cat']]

single_dqm_fn = [['A2399_ud_300_0236_20aug.fits', 'A2399_gd_60_0239_19aug.fits', 'A2399_rd_300_0142_19aug.fits'],
             ['A2670_ud_300_0409_21aug.fits', 'A2670_gd_60_0420_21aug.fits', 'A2670_rd_300_0351_19aug.fits'],
             ['A3716_ud_300_2357_21aug.fits', 'A3716_gd_300_0030_21aug.fits', 'A3716_rd_300_0418_19aug.fits']]

stack_dqm_fn = [['A2399_usd.fits', 'A2399_gsd.fits', 'A2399_rsd.fits'],
             ['A2670_usd.fits', 'A2670_gsd.fits', 'A2670_rsd.fits'],
             ['A3716_usd.fits', 'A3716_gsd.fits', 'A3716_rsd.fits']]

mag_exp_t_corr_300 = 2.5 * np.log10(300)
mag_exp_t_corr_60 = 2.5 * np.log10(60)

### Standardization parameters      mag = inst_mag [ZP:25 with flux/s] + a + b * AIRMASS
single_a = [[-1.257 + mag_exp_t_corr_300, 0.381 + mag_exp_t_corr_60, 0.567 + mag_exp_t_corr_300],
            [-1.796 + mag_exp_t_corr_300, 0.179 + mag_exp_t_corr_60, 0.567 + mag_exp_t_corr_300],
            [-1.796 + mag_exp_t_corr_300, 0.179 + mag_exp_t_corr_300, 0.567 + mag_exp_t_corr_300]]

single_b = [[-0.592 * 1.26, -0.186 * 1.26, -0.136 * 1.53],
            [-0.132 * 1.3, -0.036 * 1.26, -0.136 * 1.41],
            [-0.132 * 1.42, -0.036 * 1.31, -0.136 * 1.09]]

# stacked zero point (ZP; mag = mag0 (ZP_inst=25) + ZP)
stack_a = [
    [4.05, 6.06, 6.38],
    [3.96, 6.04, 6.22],
    [3.96, 6.15, 6.42]
]

# irsa.ipac.caltech.edu, S and F (2011)
mw_ext = [[0.159, 0.124, 0.086],
               [0.188, 0.146, 0.101],
               [0.157, 0.122, 0.085]]

am_cat = pd.read_csv(work_dir+'sex/cat/airmass.csv')

# lock = Lock()

def dqm_check_single(n):
    prev_time = start_time
    is_dqm_zero = np.ones(ni, dtype=bool)
    hdu_r = fits.open(work_dir + 'fits/best_single/' + single_dqm_fn[k][2])
    for ii in range(0, ni):
        # read the celestial coordinate
        cel_coord = [[sex_cat_r[galaxy_r]['ALPHA_J2000'][ii+n*ni],
                      sex_cat_r[galaxy_r]['DELTA_J2000'][ii+n*ni]], [0, 0]]
        # sky_coord = SkyCoord(arr_for_coord['ALPHA_J2000'][i],
        #                      arr_for_coord['DELTA_J2000'][i], unit='deg')
        # for every CCD

        # lock.acquire()
        # try:
        for jj in range(1, len(hdu_r)):
            # read WCS
            w = wcs.WCS(hdu_r[jj].header)
            pixcrd = w.wcs_world2pix(cel_coord, 1)
            # if w.footprint_contains(sky_coord):
            if (pixcrd[0][0] > 0) & (pixcrd[0][0] < hdu_r[jj].shape[1]) & \
                    (pixcrd[0][1] > 0) & (pixcrd[0][1] < hdu_r[jj].shape[0]):
                dqm_fits = hdu_r[jj].data
                # check the value of DQM
                if dqm_fits[int(pixcrd[0][1])][int(pixcrd[0][0])] % 128:
                    # toggle on saturated bool array
                    is_dqm_zero[ii] = False
        # finally:
        #     lock.release()
        if ii % 1000 == 0:
            print("--- %06.2f minutes ---" % ((time.time() - prev_time) / 60.0))
            print("--- %d / %d - %d (single r)---" % (ii, n*ni, (n+1)*ni))
            my_module.print_time()
            prev_time = time.time()
    hdu_r.close()
    return is_dqm_zero

def dqm_check_stack(n):
    prev_time = time.time()
    is_dqm_zero = np.ones(ni, dtype=bool)
    hdu_r_stack = fits.open(work_dir + 'fits/stacked/' + stack_dqm_fn[k][2])
    for ii in range(0, ni):
        # read the celestial coordinate
        cel_coord = [[sex_cat_r_stack[galaxy_r_stack]['ALPHA_J2000'][ii+n*ni],
                      sex_cat_r_stack[galaxy_r_stack]['DELTA_J2000'][ii+n*ni]], [0, 0]]
        # sky_coord = SkyCoord(arr_for_coord['ALPHA_J2000'][i],
        #                      arr_for_coord['DELTA_J2000'][i], unit='deg')
        # for every CCD
        # lock.acquire()
        # try:
        for jj in range(1, len(hdu_r_stack)):
            # read WCS
            w = wcs.WCS(hdu_r_stack[jj].header)
            pixcrd = w.wcs_world2pix(cel_coord, 1)
            # if w.footprint_contains(sky_coord):
            if (pixcrd[0][0] > 0) & (pixcrd[0][0] < hdu_r_stack[jj].shape[1]) & \
                    (pixcrd[0][1] > 0) & (pixcrd[0][1] < hdu_r_stack[jj].shape[0]):
                dqm_fits = hdu_r_stack[jj].data
                # check the value of DQM
                if dqm_fits[int(pixcrd[0][1])][int(pixcrd[0][0])] % 128:
                    # toggle on saturated bool array
                    is_dqm_zero[ii] = False
        if ii % 1000 == 0:
            print("--- %06.2f minutes ---" % ((time.time() - prev_time) / 60.0))
            print("--- %d / %d - %d (stack r)---" % (ii, n*ni, (n+1)*ni))
            my_module.print_time()
            prev_time = time.time()
        # finally:
        #     lock.release()
    hdu_r_stack.close()
    return is_dqm_zero

# for k in range(0, len(clusters)):
for k in range(0, 1):
    start_time = time.time()
    fn = work_dir+'sex/cat/DECam_19_21_aug_2014_stacked_SEx_cat_{}_match_rad_1as_Gal_ext_corrected_{}'.format(
        clusters[k], ver)

    if k == 0:  # A2399
        coords_cl_cen = SkyCoord(329.3575, -7.794722, unit='deg')
    elif k == 1: # A2670
        coords_cl_cen = SkyCoord(358.557083, -10.418889, unit='deg')
    elif k == 2: # A3716
        coords_cl_cen = SkyCoord(312.819614, -52.695408, unit='deg')

    with open(fn + '.txt', 'w') as fh, open(fn + '.reg', 'w') as regFile:
        fh.writelines('############################################################################################\n')
        fh.writelines('#   This is a catalog of photometry of galaxies detected in either single-exposure or stac-#\n')
        fh.writelines('# ked exposure of Abell cluster '+clusters[k]+' mosaics in u, g, and r band taken by Dark Energy #\n')
        fh.writelines('# Camera (DECam) mounted at the prime focus of the Blanco 4-m telescope at CTIO over       #\n')
        fh.writelines('# 19th - 21st Aug 2014. The single exposures that I used for photometry are:               #\n')
        if k == 0:
            fh.writelines('#   A2399 u-band : 300s exposure taken at 02:36 in 20-Aug observing date                   #\n')
            fh.writelines('#   A2399 g-band :  60s exposure taken at 02:39 in 19-Aug observing date                   #\n')
            fh.writelines('#   A2399 r-band : 300s exposure taken at 01:42 in 19-Aug observing date                   #\n')
        if k == 1:
            fh.writelines('#   A2670 u-band : 300s exposure taken at 04:09 in 21-Aug observing date                   #\n')
            fh.writelines('#   A2670 g-band :  60s exposure taken at 04:20 in 21-Aug observing date                   #\n')
            fh.writelines('#   A2670 r-band : 300s exposure taken at 03:51 in 19-Aug observing date                   #\n')
        if k == 2:
            fh.writelines('#   A3716 u-band : 300s exposure taken at 23:57 in 21-Aug observing date                   #\n')
            fh.writelines('#   A3716 g-band : 300s exposure taken at 00:30 in 21-Aug observing date                   #\n')
            fh.writelines('#   A3716 r-band : 300s exposure taken at 04:18 in 19-Aug observing date                   #\n')
        fh.writelines('# , which were selected based on seeing and airmass.                                       #\n')
        if k == 0:
            fh.writelines('#   Stacked mosaics with 6300s, 4200s, and 9000s exposure time for u, g, and r,            #\n')
        if k == 1:
            fh.writelines(
                '#   Stacked mosaics with 6000s, 4200s, and 10500s exposure time for u, g, and r,            #\n')
        if k == 2:
            fh.writelines(
                '#   Stacked mosaics with 3300s, 2700s, and 4800s exposure time for u, g, and r,            #\n')
        fh.writelines('# respectively, are calibrated based on photometry of the single exposures and cataloged.  #\n')
        fh.writelines('#                                                                                          #\n')
        fh.writelines('#                                  ***  Standardization  ***                               #\n')
        fh.writelines('#   All 3 bands in 3 days of observing data are standardized using standard star observati-#\n')
        fh.writelines('# ons taken right before and after the science observations (Please refer to the log file  #\n')
        fh.writelines('# [noclouds_log.xlsx]). All FITS files are available on NOAO archive. I used Source Extrac-#\n')
        fh.writelines('# tor to measure magnitudes of standard stars and galaxies. I used MAG_APER at 14.8" diame-#\n')
        fh.writelines('# ter for standard star photometry and compared with Southern Standard Stars (https://www-s#\n')
        fh.writelines('# tar.fnal.gov/Southern_ugriz/New/index.html) to convert instrument magnitudes (ZP=25) to -#\n')
        fh.writelines('# calibrated AB magnitudes correcting for atmospheric extinction by fitting the airmass te-#\n')
        fh.writelines('# rm.                                                                                      #\n')
        fh.writelines('#                                    ***  Photometry  ***                                  #\n')
        fh.writelines('#   The MAG_AUTO and MAG_ISO of galaxies from the Source Extractor are added. I calibrated #\n')
        fh.writelines('# both magnitudes using the standardization equation that I derived above. The Galactic ex-#\n')
        fh.writelines('# tinction was corrected adopting the values from irsa.ipac.caltech.edu, S and F (2011).   #\n')
        fh.writelines('#                                                                                          #\n')
        fh.writelines('#                                  ***  Sample Selection  ***                              #\n')
        fh.writelines('#  I added galaxies into this catalog if they were detected in r-band stacked mosaic and   #\n')
        fh.writelines('# Source Extractor parameter \'CLASS_STAR\' < {:4.2f} (close to 0: likely to be galaxies,       #\n'
                      .format(class_star_lim))
        fh.writelines('# close to 1: likely to be stars), and central pixel values of \'DQM\' (Data Quality Mask;   #\n')
        fh.writelines('# DECam pipeline product which was useful for getting rid of bright stars) == 0 or 128. The#\n')
        fh.writelines('# u-band and g-band magnitudes are added if there were matched sources within 1\" radii. I  #\n')
        fh.writelines('# put arbitrary number 999.0 if there were no matched source.                              #\n')
        fh.writelines('#                                                                                          #\n')
        fh.writelines('#    If you have any question, please send me an email.                                    #\n')
        fh.writelines('#    Duho Kim [email:dhkim@kasi.re.kr]                                                     #\n')
        fh.writelines('#                                                                                          #\n')
        fh.writelines('#    This Version {} was written on {}                                           #\n'
                      .format(ver, datetime.today().strftime('%Y-%m-%d')))
        fh.writelines('#                                                                                          #\n')
        fh.writelines('############################################################################################\n')
        fh.writelines('#   1 NUMBER                 Running object number                                         #\n')
        fh.writelines('#   2 ALPHA_J2000            Right ascension of barycenter (J2000)                [deg]    #\n')
        fh.writelines('#   3 DELTA_J2000            Declination of barycenter (J2000)                    [deg]    #\n')
        fh.writelines('#   4 A_WORLD                Profile RMS along major axis (world units)           [arcsec] #\n')
        fh.writelines('#   5 B_WORLD                Profile RMS along minor axis (world units)           [arcsec] #\n')
        fh.writelines('#   6 CLASS_STAR             S/G classifier output                                         #\n')
        fh.writelines('#   7 MAG_AUTO_u             Kron-like elliptical aperture magnitude in u band    [mag]    #\n')
        fh.writelines('#   8 MAGERR_AUTO_u          RMS error for AUTO magnitude in u band               [mag]    #\n')
        fh.writelines('#   9 MAG_ISO_u              Isophotal magnitude in u band                        [mag]    #\n')
        fh.writelines('#  10 MAGERR_ISO_u           RMS error for isophotal magnitude in u band          [mag]    #\n')
        fh.writelines('#  11 MAG_AUTO_g             Kron-like elliptical aperture magnitude in g band    [mag]    #\n')
        fh.writelines('#  12 MAGERR_AUTO_g          RMS error for AUTO magnitude in g band               [mag]    #\n')
        fh.writelines('#  13 MAG_ISO_g              Isophotal magnitude in g band                        [mag]    #\n')
        fh.writelines('#  14 MAGERR_ISO_g           RMS error for isophotal magnitude in g band          [mag]    #\n')
        fh.writelines('#  15 MAG_AUTO_r             Kron-like elliptical aperture magnitude in r band    [mag]    #\n')
        fh.writelines('#  16 MAGERR_AUTO_r          RMS error for AUTO magnitude in r band               [mag]    #\n')
        fh.writelines('#  17 MAG_ISO_r              Isophotal magnitude in r band                        [mag]    #\n')
        fh.writelines('#  18 MAGERR_ISO_r           RMS error for isophotal magnitude in r band          [mag]    #\n')
        fh.writelines('#  19 MAG_AUTO_u_stack       MAG_AUTO_u for stacked mosaic                        [mag]    #\n')
        fh.writelines('#  20 MAGERR_AUTO_u_stack    MAGERR_AUTO_u for stacked mosaic                     [mag]    #\n')
        fh.writelines('#  21 MAG_ISO_u_stack        MAG_ISO_u for stacked mosaic                         [mag]    #\n')
        fh.writelines('#  22 MAGERR_ISO_u_stack     MAGERR_ISO_u for stacked mosaic                      [mag]    #\n')
        fh.writelines('#  23 MAG_AUTO_g_stack       MAG_AUTO_g for stacked mosaic                        [mag]    #\n')
        fh.writelines('#  24 MAGERR_AUTO_g_stack    MAGERR_AUTO_g for stacked mosaic                     [mag]    #\n')
        fh.writelines('#  25 MAG_ISO_g_stack        MAG_ISO_g for stacked mosaic                         [mag]    #\n')
        fh.writelines('#  26 MAGERR_ISO_g_stack     MAGERR_ISO_g for stacked mosaic                      [mag]    #\n')
        fh.writelines('#  27 MAG_AUTO_r_stack       MAG_AUTO_r for stacked mosaic                        [mag]    #\n')
        fh.writelines('#  28 MAGERR_AUTO_r_stack    MAGERR_AUTO_r for stacked mosaic                     [mag]    #\n')
        fh.writelines('#  29 MAG_ISO_r_stack        MAG_ISO_r for stacked mosaic                         [mag]    #\n')
        fh.writelines('#  30 MAGERR_ISO_r_stack     MAGERR_ISO_r for stacked mosaic                      [mag]    #\n')
        fh.writelines('############################################################################################\n')

        # read SEx cat
        sex_cat_u = ascii.read(work_dir + 'sex/cat/best_single/' + single_fn[k][0])
        sex_cat_g = ascii.read(work_dir + 'sex/cat/best_single/' + single_fn[k][1])
        sex_cat_r = ascii.read(work_dir + 'sex/cat/best_single/' + single_fn[k][2])

        sex_cat_u_stack = ascii.read(work_dir + 'sex/cat/stack/' + clusters[k] + '_usi.cat')
        sex_cat_g_stack = ascii.read(work_dir + 'sex/cat/stack/' + clusters[k] + '_gsi.cat')
        sex_cat_r_stack = ascii.read(work_dir + 'sex/cat/stack/' + clusters[k] + '_rsi.cat')

        # standardize magnitudes and Milky Way extinction correction
        sex_cat_u[mag_sys] = sex_cat_u[mag_sys] + single_a[k][0] + single_b[k][0] - mw_ext[k][0]
        sex_cat_g[mag_sys] = sex_cat_g[mag_sys] + single_a[k][1] + single_b[k][1] - mw_ext[k][1]
        sex_cat_r[mag_sys] = sex_cat_r[mag_sys] + single_a[k][2] + single_b[k][2] - mw_ext[k][2]

        sex_cat_u[mag_sys1] = sex_cat_u[mag_sys1] + single_a[k][0] + single_b[k][0] - mw_ext[k][0]
        sex_cat_g[mag_sys1] = sex_cat_g[mag_sys1] + single_a[k][1] + single_b[k][1] - mw_ext[k][1]
        sex_cat_r[mag_sys1] = sex_cat_r[mag_sys1] + single_a[k][2] + single_b[k][2] - mw_ext[k][2]

        sex_cat_u_stack[mag_sys] = sex_cat_u_stack[mag_sys] + stack_a[k][0] - mw_ext[k][0]
        sex_cat_g_stack[mag_sys] = sex_cat_g_stack[mag_sys] + stack_a[k][1] - mw_ext[k][1]
        sex_cat_r_stack[mag_sys] = sex_cat_r_stack[mag_sys] + stack_a[k][2] - mw_ext[k][2]

        sex_cat_u_stack[mag_sys1] = sex_cat_u_stack[mag_sys1] + stack_a[k][0] - mw_ext[k][0]
        sex_cat_g_stack[mag_sys1] = sex_cat_g_stack[mag_sys1] + stack_a[k][1] - mw_ext[k][1]
        sex_cat_r_stack[mag_sys1] = sex_cat_r_stack[mag_sys1] + stack_a[k][2] - mw_ext[k][2]

        # for every catalog item
        for i in range(0, len(sex_cat_u)):
            # make 'NUMBER' unique
            sex_cat_u['NUMBER'][i] = i
        for i in range(0, len(sex_cat_g)):
            sex_cat_g['NUMBER'][i] = i
        for i in range(0, len(sex_cat_r)):
            sex_cat_r['NUMBER'][i] = i
        for i in range(0, len(sex_cat_u_stack)):
            sex_cat_u_stack['NUMBER'][i] = i
        for i in range(0, len(sex_cat_g_stack)):
            sex_cat_g_stack['NUMBER'][i] = i
        for i in range(0, len(sex_cat_r_stack)):
            sex_cat_r_stack['NUMBER'][i] = i

        ### exclude stars ###
        galaxy_r = (sex_cat_r['CLASS_STAR'] < class_star_lim)
        galaxy_r_stack = (sex_cat_r_stack['CLASS_STAR'] < class_star_lim)
        # crit_r = sex_cat_r[criteria_r]

        ### exclude saturated sources in single-exposure r band ###
        num_process = 12
        num_gal_r = len(sex_cat_r[galaxy_r])
        ni = int(num_gal_r / num_process)
        with Pool(num_process) as p:
            is_dqm_r_zero = p.map(dqm_check_single, range(num_process))
        is_dqm_r_zero = np.array(is_dqm_r_zero).reshape(-1)

        num_gal_r_stack = len(sex_cat_r_stack[galaxy_r_stack])
        ni = int(num_gal_r_stack / num_process)
        with Pool(num_process) as p:
            is_dqm_r_zero_stack = p.map(dqm_check_stack, range(num_process))
        is_dqm_r_zero_stack = np.array(is_dqm_r_zero_stack).reshape(-1)

        ### make ref_r to match other bands ###
        coords_single = SkyCoord(sex_cat_r[galaxy_r][0:len(is_dqm_r_zero)][is_dqm_r_zero]['ALPHA_J2000'],
                                 sex_cat_r[galaxy_r][0:len(is_dqm_r_zero)][is_dqm_r_zero]['DELTA_J2000'], unit='deg')
        coords_stack = SkyCoord(sex_cat_r_stack[galaxy_r_stack][0:len(is_dqm_r_zero_stack)][is_dqm_r_zero_stack]['ALPHA_J2000'],
                                sex_cat_r_stack[galaxy_r_stack][0:len(is_dqm_r_zero_stack)][is_dqm_r_zero_stack]['DELTA_J2000'], unit='deg')
        idx_r_single, d2d, d3d = coords_single.match_to_catalog_sky(coords_stack)
        matched_r_single = (d2d.arcsec < max_sep)

        ref_r = vstack([sex_cat_r_stack[galaxy_r_stack][0:len(is_dqm_r_zero_stack)][is_dqm_r_zero_stack],
                          sex_cat_r[galaxy_r][0:len(is_dqm_r_zero)][is_dqm_r_zero][~matched_r_single]])
        num_stack = len(sex_cat_r_stack[galaxy_r_stack][0:len(is_dqm_r_zero_stack)][is_dqm_r_zero_stack])
        num_all = len(ref_r)

        # detected in both = 0, single = 1, stack = 2
        # single_or_stack = np.zeros(num_all)
        # for i in range(0, num_all):
        #     if i < num_stack:
        #         if i not in idx_r_single[matched_r_single]:
        #             single_or_stack[i] = 2
        #     else:
        #         single_or_stack[i] = 1
        #
        # ref_r.add_column(single_or_stack, name='detect', index=0)
        ref_r.sort(mag_sys1)

        for i in range(0, len(ref_r)):
            # make 'NUMBER' unique
            ref_r['NUMBER'][i] = i

        coords_r = SkyCoord(ref_r['ALPHA_J2000'], ref_r['DELTA_J2000'], unit='deg')

        ### match r band again ###
        coords_r_single = SkyCoord(sex_cat_r[galaxy_r][0:len(is_dqm_r_zero)][is_dqm_r_zero]['ALPHA_J2000'],
                                   sex_cat_r[galaxy_r][0:len(is_dqm_r_zero)][is_dqm_r_zero]['DELTA_J2000'], unit='deg')
        coords_r_stack = SkyCoord(sex_cat_r_stack[galaxy_r_stack][0:len(is_dqm_r_zero_stack)][is_dqm_r_zero_stack]['ALPHA_J2000'],
                                  sex_cat_r_stack[galaxy_r_stack][0:len(is_dqm_r_zero_stack)][is_dqm_r_zero_stack]['DELTA_J2000'], unit='deg')
        idx_r_single, d2d_r_single, d3d = coords_r_single.match_to_catalog_sky(coords_r)
        idx_r_stack, d2d_r_stack, d3d = coords_r_stack.match_to_catalog_sky(coords_r)
        matched_r_single = (d2d_r_single.arcsec < max_sep)
        matched_r_stack = (d2d_r_stack.arcsec < max_sep)

        ### match u band ###
        coords_u_single = SkyCoord(sex_cat_u['ALPHA_J2000'],
                                   sex_cat_u['DELTA_J2000'], unit='deg')
        coords_u_stack = SkyCoord(sex_cat_u_stack['ALPHA_J2000'],
                                  sex_cat_u_stack['DELTA_J2000'], unit='deg')
        idx_u_single, d2d_u_single, d3d = coords_u_single.match_to_catalog_sky(coords_r)
        idx_u_stack, d2d_u_stack, d3d = coords_u_stack.match_to_catalog_sky(coords_r)
        matched_u_single = (d2d_u_single.arcsec < max_sep)
        matched_u_stack = (d2d_u_stack.arcsec < max_sep)

        ### match g band ###
        coords_g_single = SkyCoord(sex_cat_g['ALPHA_J2000'],
                                   sex_cat_g['DELTA_J2000'], unit='deg')
        coords_g_stack = SkyCoord(sex_cat_g_stack['ALPHA_J2000'],
                                  sex_cat_g_stack['DELTA_J2000'], unit='deg')
        idx_g_single, d2d_g_single, d3d = coords_g_single.match_to_catalog_sky(coords_r)
        idx_g_stack, d2d_g_stack, d3d = coords_g_stack.match_to_catalog_sky(coords_r)
        matched_g_single = (d2d_g_single.arcsec < max_sep)
        matched_g_stack = (d2d_g_stack.arcsec < max_sep)

        for i in range(0, num_all):
            prev_time = time.time()
            if i in idx_r_single[matched_r_single]:                                      #### for r band
                idx_match_r_single = np.where(idx_r_single == ref_r['NUMBER'][i])
                r_iso_single = sex_cat_r[galaxy_r][0:len(is_dqm_r_zero)][is_dqm_r_zero][mag_sys][idx_match_r_single][0]
                rerr_iso_single = sex_cat_r[galaxy_r][0:len(is_dqm_r_zero)][is_dqm_r_zero][magerr_sys][idx_match_r_single][0]
                r_auto_single = sex_cat_r[galaxy_r][0:len(is_dqm_r_zero)][is_dqm_r_zero][mag_sys1][idx_match_r_single][0]
                rerr_auto_single = sex_cat_r[galaxy_r][0:len(is_dqm_r_zero)][is_dqm_r_zero][magerr_sys1][idx_match_r_single][0]
            else:
                r_iso_single = 999.0
                rerr_iso_single = 999.0
                r_auto_single = 999.0
                rerr_auto_single = 999.0
            if i in idx_r_stack[matched_r_stack]:
                idx_match_r = np.where(idx_r_stack == ref_r['NUMBER'][i])
                r_iso = sex_cat_r_stack[galaxy_r_stack][0:len(is_dqm_r_zero_stack)][is_dqm_r_zero_stack][matched_r_stack][mag_sys][idx_match_r][0]
                rerr_iso = sex_cat_r_stack[galaxy_r_stack][0:len(is_dqm_r_zero_stack)][is_dqm_r_zero_stack][matched_r_stack][magerr_sys][idx_match_r][0]
                r_auto = sex_cat_r_stack[galaxy_r_stack][0:len(is_dqm_r_zero_stack)][is_dqm_r_zero_stack][matched_r_stack][mag_sys1][idx_match_r][0]
                rerr_auto = sex_cat_r_stack[galaxy_r_stack][0:len(is_dqm_r_zero_stack)][is_dqm_r_zero_stack][matched_r_stack][magerr_sys1][idx_match_r][0]
            else:
                r_iso = 999.0
                rerr_iso = 999.0
                r_auto = 999.0
                rerr_auto = 999.0
            if i in idx_u_single[matched_u_single]:                                      #### for u band
                idx_match_u_single = np.where(idx_u_single[matched_u_single] == ref_r['NUMBER'][i])
                u_iso_single = sex_cat_u[matched_u_single][mag_sys][idx_match_u_single][0]
                uerr_iso_single = sex_cat_u[matched_u_single][magerr_sys][idx_match_u_single][0]
                u_auto_single = sex_cat_u[matched_u_single][mag_sys1][idx_match_u_single][0]
                uerr_auto_single = sex_cat_u[matched_u_single][magerr_sys1][idx_match_u_single][0]
            else:
                u_iso_single = 999.0
                uerr_iso_single = 999.0
                u_auto_single = 999.0
                uerr_auto_single = 999.0
            if i in idx_u_stack[matched_u_stack]:
                idx_match_u = np.where(idx_u_stack[matched_u_stack] == ref_r['NUMBER'][i])
                u_iso = sex_cat_u_stack[matched_u_stack][mag_sys][idx_match_u][0]
                uerr_iso = sex_cat_u_stack[matched_u_stack][magerr_sys][idx_match_u][0]
                u_auto = sex_cat_u_stack[matched_u_stack][mag_sys1][idx_match_u][0]
                uerr_auto = sex_cat_u_stack[matched_u_stack][magerr_sys1][idx_match_u][0]
            else:
                u_iso = 999.0
                uerr_iso = 999.0
                u_auto = 999.0
                uerr_auto = 999.0
            if i in idx_g_stack[matched_g_stack]:                                        #### for g band
                idx_match_g = np.where(idx_g_stack[matched_g_stack] == ref_r['NUMBER'][i])
                g_iso = sex_cat_g_stack[matched_g_stack][mag_sys][idx_match_g][0]
                gerr_iso = sex_cat_g_stack[matched_g_stack][magerr_sys][idx_match_g][0]
                g_auto = sex_cat_g_stack[matched_g_stack][mag_sys1][idx_match_g][0]
                gerr_auto = sex_cat_g_stack[matched_g_stack][magerr_sys1][idx_match_g][0]
            else:
                g_iso = 999.0
                gerr_iso = 999.0
                g_auto = 999.0
                gerr_auto = 999.0
            if i in idx_g_single[matched_g_single]:
                idx_match_g_single = np.where(idx_g_single[matched_g_single] == ref_r['NUMBER'][i])
                g_iso_single = sex_cat_g[matched_g_single][mag_sys][idx_match_g_single][0]
                gerr_iso_single = sex_cat_g[matched_g_single][magerr_sys][idx_match_g_single][0]
                g_auto_single = sex_cat_g[matched_g_single][mag_sys1][idx_match_g_single][0]
                gerr_auto_single = sex_cat_g[matched_g_single][magerr_sys1][idx_match_g_single][0]
            else:
                g_iso_single = 999.0
                gerr_iso_single = 999.0
                g_auto_single = 999.0
                gerr_auto_single = 999.0

            fh.writelines("{:4}".format(i+1) + ' ' +                                    # NUMBER
                "{:12.7f}".format(coords_r[i].ra.value)+' '+                            # RA
                "{:12.7f}".format(coords_r[i].dec.value) + ' ' +                        # DEC
                "{:7.3f}".format(ref_r['A_WORLD'][i]*3600) + ' ' +      # A_WORLD
                "{:7.3f}".format(ref_r['B_WORLD'][i]*3600) + ' ' +      # B_WORLD
                "{:4.2f}".format(ref_r['CLASS_STAR'][i]) + ' ' +        # CLASS_STAR
                "{:7.3f}".format(u_auto_single) + ' ' +                                        # u-band 'MAG_AUTO'
                "{:7.3f}".format(uerr_auto_single) + ' ' +                                     # u-band 'MAGERR_AUTO'
                "{:7.3f}".format(u_iso_single) + ' ' +                                         # u-band 'MAG_ISO'
                "{:7.3f}".format(uerr_iso_single) + ' ' +                                      # u-band 'MAGERR_ISO'
                "{:7.3f}".format(g_auto_single) + ' ' +                                        # g-band 'MAG_AUTO'
                "{:7.3f}".format(gerr_auto_single) + ' ' +                                     # g-band 'MAGERR_AUTO'
                "{:7.3f}".format(g_iso_single) + ' ' +                                         # g-band 'MAG_ISO'
                "{:7.3f}".format(gerr_iso_single) + ' ' +                                      # g-band 'MAGERR_ISO'
                "{:7.3f}".format(r_auto_single) + ' ' +            # r-band 'MAG_AUTO'
                "{:7.3f}".format(rerr_auto_single) + ' ' +         # r-band 'MAGERR_AUTO'
                "{:7.3f}".format(r_iso_single) + ' ' +             # r-band 'MAG_ISO'
                "{:7.3f}".format(rerr_iso_single) + ' ' +         # r-band 'MAGERR_ISO'
                "{:7.3f}".format(u_auto) + ' ' +  # u-band 'MAG_AUTO'
                "{:7.3f}".format(uerr_auto) + ' ' +  # u-band 'MAGERR_AUTO'
                "{:7.3f}".format(u_iso) + ' ' +  # u-band 'MAG_ISO'
                "{:7.3f}".format(uerr_iso) + ' ' +  # u-band 'MAGERR_ISO'
                "{:7.3f}".format(g_auto) + ' ' +  # g-band 'MAG_AUTO'
                "{:7.3f}".format(gerr_auto) + ' ' +  # g-band 'MAGERR_AUTO'
                "{:7.3f}".format(g_iso) + ' ' +  # g-band 'MAG_ISO'
                "{:7.3f}".format(gerr_iso) + ' ' +  # g-band 'MAGERR_ISO'
                "{:7.3f}".format(r_auto) + ' ' +  # r-band 'MAG_AUTO'
                "{:7.3f}".format(rerr_auto) + ' ' +  # r-band 'MAGERR_AUTO'
                "{:7.3f}".format(r_iso) + ' ' +  # r-band 'MAG_ISO'
                "{:7.3f}".format(rerr_iso) + '\n')  # r-band 'MAGERR_ISO'
            if ref_r['FLAGS'][i] > 0:
                color = 'red'
            else:
                color = 'green'
            regFile.writelines("j2000; ellipse({:12.7f}, {:12.7f}, {:7.3f}\", {:7.3f}\", {:7.3f}) # text=\'{}\', "
                               "color={}, dash=1 \n".format(
                coords_r[i].ra.value,
                coords_r[i].dec.value,
                ref_r['A_WORLD'][i]*3600,
                ref_r['B_WORLD'][i]*3600,
                180-ref_r['THETA_WORLD'][i],
                i+1,
                color))

            if i % 10000 == 0:
                print("--- %06.2f minutes ---" % (((time.time() - prev_time))/60.0))
                print("--- %d / %d (write cat) ---" % (i, num_all))
                my_module.print_time()
                prev_time = time.time()

    print("--- {} Done ---".format(clusters[k]))
    my_module.print_time()