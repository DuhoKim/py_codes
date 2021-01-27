from astropy.io import ascii
from astropy.table import vstack
import numpy as np
import time
from astropy.coordinates import SkyCoord
import my_module
import pandas as pd
from astropy.io import fits
from astropy import wcs
from datetime import datetime
from multiprocessing import Pool
import abell_cluster_module as ab

tag = 'deblend'

mw_ext = [[0., 0., 0.],             # corrected already
    [0., 0., 0.],
    [0., 0., 0.],
    [0., 0., 0.],
    [0., 0., 0.],
    [0., 0., 0.],
    [0., 0., 0.]]

am_cat = pd.read_csv(ab.work_dir+'sex/cat/airmass.csv')

def dqm_check_single(n):
    prev_time = start_time
    is_dqm_zero = np.ones(ni, dtype=bool)
    hdu_r = fits.open(ab.work_dir + 'fits/best_single/' + ab.short_dqm_fn[k][2])
    for ii in range(0, ni):
        # read the celestial coordinate
        cel_coord = [[sex_cat_r[galaxy_r]['ALPHA_J2000'][ii+n*ni],
                      sex_cat_r[galaxy_r]['DELTA_J2000'][ii+n*ni]], [0, 0]]

        # a few bright galaxies have saturated center but we want to include those so...
        # if (sex_cat_r[galaxy_r][ab.mag_sys][ii+n*ni] < 20) &  (sex_cat_r[galaxy_r]['CLASS_STAR'][ii+n*ni] < 0.2):
        #     continue
        if (sex_cat_r[galaxy_r][ab.mag_sys][ii+n*ni] < 19) &  (sex_cat_r[galaxy_r]['CLASS_STAR'][ii+n*ni] < 0.1):
            continue

        for jj in range(1, len(hdu_r)):
            # read WCS
            w = wcs.WCS(hdu_r[jj].header)
            pixcrd = w.wcs_world2pix(cel_coord, 1)
            # if w.footprint_contains(sky_coord):
            if (pixcrd[0][0] > 0) & (pixcrd[0][0] < hdu_r[jj].shape[1]) & \
                    (pixcrd[0][1] > 0) & (pixcrd[0][1] < hdu_r[jj].shape[0]):
                dqm_fits = hdu_r[jj].data
                # check the value of DQM
                if k == 0:      # Somehow A754 has value 7 (cosmic ray) in galaxies
                    if dqm_fits[int(pixcrd[0][1])][int(pixcrd[0][0])] % 7:
                        # toggle on saturated bool array
                        is_dqm_zero[ii] = False
                else:           # somehow other clusters have value 128 (?) in galaxies
                    if dqm_fits[int(pixcrd[0][1])][int(pixcrd[0][0])] % 128:
                        # toggle on saturated bool array
                        is_dqm_zero[ii] = False

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
    hdu_r_stack = fits.open(ab.work_dir + 'fits/best_single/' + ab.short_dqm_fn[k][2])
    for ii in range(0, ni):
        # read the celestial coordinate
        cel_coord = [[sex_cat_r_stack[galaxy_r_stack]['ALPHA_J2000'][ii+n*ni],
                      sex_cat_r_stack[galaxy_r_stack]['DELTA_J2000'][ii+n*ni]], [0, 0]]

        for jj in range(1, len(hdu_r_stack)):
            # read WCS
            w = wcs.WCS(hdu_r_stack[jj].header)
            pixcrd = w.wcs_world2pix(cel_coord, 1)
            # if w.footprint_contains(sky_coord):
            if (pixcrd[0][0] > 0) & (pixcrd[0][0] < hdu_r_stack[jj].shape[1]) & \
                    (pixcrd[0][1] > 0) & (pixcrd[0][1] < hdu_r_stack[jj].shape[0]):
                dqm_fits = hdu_r_stack[jj].data
                # check the value of DQM
                if k == 0:      # Somehow A754 has value 7 (cosmic ray) in galaxies
                    if dqm_fits[int(pixcrd[0][1])][int(pixcrd[0][0])] % 7:
                        # toggle on saturated bool array
                        is_dqm_zero[ii] = False
                else:           # somehow other clusters have value 128 (?) in galaxies
                    if dqm_fits[int(pixcrd[0][1])][int(pixcrd[0][0])] % 128:
                        # toggle on saturated bool array
                        is_dqm_zero[ii] = False
        if ii % 1000 == 0:
            print("--- %06.2f minutes ---" % ((time.time() - prev_time) / 60.0))
            print("--- %d / %d - %d (stack r)---" % (ii, n*ni, (n+1)*ni))
            my_module.print_time()
            prev_time = time.time()

    hdu_r_stack.close()
    return is_dqm_zero

for k in range(0, len(ab.clusters)):
# for k in range(0, 1):
    start_time = time.time()
    fn = ab.work_dir+'sex/cat/DECam_merged_SEx_cat_{}_Gal_ext_corrected_{}'.format(
        ab.clusters[k], ab.ver)+'_'+tag

    with open(fn + '_single.txt', 'w') as fh1, open(fn + '_stack.txt', 'w') as fh2, open(fn + '_merged.txt', 'w') as fh, open(fn + '_kron.reg', 'w') as regFile:
        regFile.writelines('# Region file format: DS9 version 4.1\n')
        regFile.writelines('global color=green dashlist=8 3 width=2 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
        regFile.writelines('fk5\n')

        fh.writelines('############################################################################################\n')
        fh.writelines('#   This is a catalog of photometry of galaxies detected in either single-exposure or stac-#\n')
        fh.writelines(
            '# ked exposure of Abell cluster ' + ab.clusters[k] + ' mosaics in u, g, and r band taken by Dark Energy     #\n')
        fh.writelines('# Camera (DECam) mounted at the prime focus of the Blanco 4-m telescope at CTIO over       #\n')
        fh.writelines('# 19th - 21st Aug 2014. The single exposures that I used for photometry are:               #\n')
        if k == 0:
            fh.writelines(
                '#   A2399 u-band : 300s exposure taken at 02:36 in 20-Aug observing date                   #\n')
            fh.writelines(
                '#   A2399 g-band :  60s exposure taken at 02:39 in 19-Aug observing date                   #\n')
            fh.writelines(
                '#   A2399 r-band : 300s exposure taken at 01:42 in 19-Aug observing date                   #\n')
        if k == 1:
            fh.writelines(
                '#   A2670 u-band : 300s exposure taken at 04:09 in 21-Aug observing date                   #\n')
            fh.writelines(
                '#   A2670 g-band :  60s exposure taken at 04:20 in 21-Aug observing date                   #\n')
            fh.writelines(
                '#   A2670 r-band : 300s exposure taken at 03:51 in 19-Aug observing date                   #\n')
        if k == 2:
            fh.writelines(
                '#   A3716 u-band : 300s exposure taken at 23:57 in 21-Aug observing date                   #\n')
            fh.writelines(
                '#   A3716 g-band : 300s exposure taken at 00:30 in 21-Aug observing date                   #\n')
            fh.writelines(
                '#   A3716 r-band : 300s exposure taken at 04:18 in 19-Aug observing date                   #\n')
        fh.writelines('# , which were selected based on seeing and airmass.                                       #\n')
        if k == 0:
            fh.writelines(
                '#   Stacked mosaics with 6300s, 4200s, and 9000s exposure time for u, g, and r,            #\n')
        if k == 1:
            fh.writelines(
                '#   Stacked mosaics with 6000s, 4200s, and 10500s exposure time for u, g, and r,           #\n')
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
        fh.writelines(
            '# Source Extractor parameter \'CLASS_STAR\' < {:4.2f} (close to 0: likely to be galaxies,       #\n'
            .format(ab.class_star_lim))
        fh.writelines(
            '# close to 1: likely to be stars), and central pixel values of \'DQM\' (Data Quality Mask;   #\n')
        fh.writelines('# DECam pipeline product which was useful for getting rid of bright stars) == 0 or 128. The#\n')
        fh.writelines('# u-band and g-band magnitudes are added if there were matched sources within 1\" radii. I  #\n')
        fh.writelines('# put arbitrary number 999.0 if there were no matched source.                              #\n')
        fh.writelines('#                                                                                          #\n')
        fh.writelines('#    If you have any question, please send me an email.                                    #\n')
        fh.writelines('#    Duho Kim [email:dhkim@kasi.re.kr]                                                     #\n')
        fh.writelines('#                                                                                          #\n')
        fh.writelines('#    This Version {} was written on {}                               #\n'
                      .format(ab.ver, datetime.today().strftime('%Y-%m-%d')))
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
        fh.writelines('#   9 MAG_AUTO_g             Kron-like elliptical aperture magnitude in g band    [mag]    #\n')
        fh.writelines('#  10 MAGERR_AUTO_g          RMS error for AUTO magnitude in g band               [mag]    #\n')
        fh.writelines('#  11 MAG_AUTO_r             Kron-like elliptical aperture magnitude in r band    [mag]    #\n')
        fh.writelines('#  12 MAGERR_AUTO_r          RMS error for AUTO magnitude in r band               [mag]    #\n')
        fh.writelines('#  13 FROM_WHERE_u           Source of MAG_AUTO_u     [0:single, 1:stacked, 2:no match]    #\n')
        fh.writelines('#  14 FROM_WHERE_g           Source of MAG_AUTO_g     [0:single, 1:stacked, 2:no match]    #\n')
        fh.writelines('#  15 FROM_WHERE_r           Source of MAG_AUTO_r     [0:single, 1:stacked, 2:no match]    #\n')
        fh.writelines('############################################################################################\n')

        # read SEx cat
        sex_cat_u = ascii.read(ab.short_cat_dir + ab.short_cat_fn[k][0] + '_' + tag + '.cat')
        sex_cat_g = ascii.read(ab.short_cat_dir + ab.short_cat_fn[k][1] + '_' + tag + '.cat')
        sex_cat_r = ascii.read(ab.short_cat_dir + ab.short_cat_fn[k][2] + '_' + tag + '.cat')

        sex_cat_u_stack = ascii.read(ab.work_dir + 'sex/cat/stack/' + ab.clusters[k] + '_usi.cat')
        sex_cat_g_stack = ascii.read(ab.work_dir + 'sex/cat/stack/' + ab.clusters[k] + '_gsi.cat')
        sex_cat_r_stack = ascii.read(ab.work_dir + 'sex/cat/stack/' + ab.clusters[k] + '_rsi.cat')

        # add column for FROM_WHERE parameter
        sex_cat_u.add_column(np.zeros(len(sex_cat_u)), name='FROM_WHERE')
        sex_cat_g.add_column(np.zeros(len(sex_cat_g)), name='FROM_WHERE')
        sex_cat_r.add_column(np.zeros(len(sex_cat_r)), name='FROM_WHERE')

        sex_cat_u_stack.add_column(np.ones(len(sex_cat_u_stack)), name='FROM_WHERE')
        sex_cat_g_stack.add_column(np.ones(len(sex_cat_g_stack)), name='FROM_WHERE')
        sex_cat_r_stack.add_column(np.ones(len(sex_cat_r_stack)), name='FROM_WHERE')

        # standardize magnitudes and Milky Way extinction correction
        sex_cat_u[ab.mag_sys] = sex_cat_u[ab.mag_sys] + ab.short_a[k][0] + ab.short_b[k][0] - mw_ext[k][0]
        sex_cat_g[ab.mag_sys] = sex_cat_g[ab.mag_sys] + ab.short_a[k][1] + ab.short_b[k][1] - mw_ext[k][1]
        sex_cat_r[ab.mag_sys] = sex_cat_r[ab.mag_sys] + ab.short_a[k][2] + ab.short_b[k][2] - mw_ext[k][2]

        sex_cat_u_stack[ab.mag_sys] = sex_cat_u_stack[ab.mag_sys] + ab.stack_a[k][0] - mw_ext[k][0]
        sex_cat_g_stack[ab.mag_sys] = sex_cat_g_stack[ab.mag_sys] + ab.stack_a[k][1] - mw_ext[k][1]
        sex_cat_r_stack[ab.mag_sys] = sex_cat_r_stack[ab.mag_sys] + ab.stack_a[k][2] - mw_ext[k][2]

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

        ### exclude stars and magnitudes less than mag lim and sources outside 0.8deg circle ###
        coords_r_single = SkyCoord(sex_cat_r['ALPHA_J2000'], sex_cat_r['DELTA_J2000'], unit='deg')
        coords_r_stack = SkyCoord(sex_cat_r_stack['ALPHA_J2000'], sex_cat_r_stack['DELTA_J2000'], unit='deg')

        galaxy_r = (sex_cat_r['CLASS_STAR'] < ab.class_star_lim) & (sex_cat_r[ab.mag_sys] < ab.mag_lim) & \
                   (coords_r_single.separation(ab.coords_cl_cen[k]).value < ab.rad_lim)
        galaxy_r_stack = (sex_cat_r_stack['CLASS_STAR'] < ab.class_star_lim) & (sex_cat_r_stack[ab.mag_sys] < ab.mag_lim) & \
                         (coords_r_stack.separation(ab.coords_cl_cen[k]).value < ab.rad_lim)
        # crit_r = sex_cat_r[criteria_r]

        ### exclude saturated sources in single-exposure r band ###
        num_process = 12
        num_gal_r = len(sex_cat_r[galaxy_r])
        ni = int(num_gal_r / num_process)
        with Pool(num_process) as p:
            is_dqm_r_zero = p.map(dqm_check_single, range(num_process))
        is_dqm_r_zero = np.array(is_dqm_r_zero).reshape(-1)

        for i in range(0, sum(is_dqm_r_zero)):
            fh1.writelines("{:4}".format(i + 1) + ' ' +  # NUMBER
                          "{:12.7f}".format(sex_cat_r[galaxy_r][0:len(is_dqm_r_zero)][is_dqm_r_zero][i]['ALPHA_J2000']) + ' ' +  # RA
                          "{:12.7f}".format(sex_cat_r[galaxy_r][0:len(is_dqm_r_zero)][is_dqm_r_zero][i]['DELTA_J2000']) + ' ' +  # DEC
                          "{:7.3f}".format(sex_cat_r[galaxy_r][0:len(is_dqm_r_zero)][is_dqm_r_zero][i]['MAG_AUTO']) + ' ' +  # r
                          "{:7.3f}".format(sex_cat_r[galaxy_r][0:len(is_dqm_r_zero)][is_dqm_r_zero][i]['MAGERR_AUTO']) + '\n')  # rerr

        num_gal_r_stack = len(sex_cat_r_stack[galaxy_r_stack])
        ni = int(num_gal_r_stack / num_process)
        with Pool(num_process) as p:
            is_dqm_r_zero_stack = p.map(dqm_check_stack, range(num_process))
        is_dqm_r_zero_stack = np.array(is_dqm_r_zero_stack).reshape(-1)

        ### make ref_r to match other bands ###
        coords_single = SkyCoord(sex_cat_r[galaxy_r][0:len(is_dqm_r_zero)][is_dqm_r_zero]['ALPHA_J2000'],
                                 sex_cat_r[galaxy_r][0:len(is_dqm_r_zero)][is_dqm_r_zero]['DELTA_J2000'],
                                 unit='deg')
        coords_stack = SkyCoord(
            sex_cat_r_stack[galaxy_r_stack][0:len(is_dqm_r_zero_stack)][is_dqm_r_zero_stack]['ALPHA_J2000'],
            sex_cat_r_stack[galaxy_r_stack][0:len(is_dqm_r_zero_stack)][is_dqm_r_zero_stack]['DELTA_J2000'],
            unit='deg')
        idx_r_stack, d2d, d3d = coords_stack.match_to_catalog_sky(coords_single)
        matched_r_stack = (d2d.arcsec < ab.max_sep)

        # ref_r = vstack([sex_cat_r_stack[galaxy_r_stack][0:len(is_dqm_r_zero_stack)][is_dqm_r_zero_stack],
        #                 sex_cat_r[galaxy_r][0:len(is_dqm_r_zero)][is_dqm_r_zero][~matched_r_single]])
        # num_stack = len(sex_cat_r_stack[galaxy_r_stack][0:len(is_dqm_r_zero_stack)][is_dqm_r_zero_stack])
        # num_all = len(ref_r)

        ref_r = vstack([sex_cat_r[galaxy_r][0:len(is_dqm_r_zero)][is_dqm_r_zero],
                        sex_cat_r_stack[galaxy_r_stack][0:len(is_dqm_r_zero_stack)][is_dqm_r_zero_stack][~matched_r_stack]])
        num_single = len(sex_cat_r[galaxy_r][0:len(is_dqm_r_zero)][is_dqm_r_zero])
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
        ref_r.sort(ab.mag_sys)

        for i in range(0, len(ref_r)):
            # make 'NUMBER' unique
            ref_r['NUMBER'][i] = i

        coords_r = SkyCoord(ref_r['ALPHA_J2000'], ref_r['DELTA_J2000'], unit='deg')

        ### match u band ###
        coords_u_single = SkyCoord(sex_cat_u['ALPHA_J2000'],
                                   sex_cat_u['DELTA_J2000'], unit='deg')
        coords_u_stack = SkyCoord(sex_cat_u_stack['ALPHA_J2000'],
                                  sex_cat_u_stack['DELTA_J2000'], unit='deg')
        idx_u_single, d2d_u_single, d3d = coords_u_single.match_to_catalog_sky(coords_r)
        idx_u_stack, d2d_u_stack, d3d = coords_u_stack.match_to_catalog_sky(coords_r)
        matched_u_single = (d2d_u_single.arcsec < ab.max_sep)
        matched_u_stack = (d2d_u_stack.arcsec < ab.max_sep)

        ### match g band ###
        coords_g_single = SkyCoord(sex_cat_g['ALPHA_J2000'],
                                   sex_cat_g['DELTA_J2000'], unit='deg')
        coords_g_stack = SkyCoord(sex_cat_g_stack['ALPHA_J2000'],
                                  sex_cat_g_stack['DELTA_J2000'], unit='deg')
        idx_g_single, d2d_g_single, d3d = coords_g_single.match_to_catalog_sky(coords_r)
        idx_g_stack, d2d_g_stack, d3d = coords_g_stack.match_to_catalog_sky(coords_r)
        matched_g_single = (d2d_g_single.arcsec < ab.max_sep)
        matched_g_stack = (d2d_g_stack.arcsec < ab.max_sep)

        ### write catalog for stacked mosaics only ###
        idx_u2r, d2d_r2u, d3d = coords_stack.match_to_catalog_sky(coords_u_stack)
        idx_g2r, d2d_r2g, d3d = coords_stack.match_to_catalog_sky(coords_g_stack)
        match_r2u = (d2d_r2u.arcsec < ab.max_sep)
        match_r2g = (d2d_r2g.arcsec < ab.max_sep)

        for i in range(0, sum(is_dqm_r_zero_stack)):
            if match_r2u[i]:
                umag = sex_cat_u_stack[idx_u2r][i][ab.mag_sys]
                uerr = sex_cat_u_stack[idx_u2r][i][ab.magerr_sys]
            else:
                umag = 999.0
                uerr = 999.0
            if match_r2g[i]:
                gmag = sex_cat_g_stack[idx_g2r][i][ab.mag_sys]
                gerr = sex_cat_g_stack[idx_g2r][i][ab.magerr_sys]
            else:
                gmag = 999.0
                gerr = 999.0
            fh2.writelines("{:4}".format(i + 1) + ' ' +  # NUMBER
                            "{:12.7f}".format(sex_cat_r_stack[galaxy_r_stack][0:len(is_dqm_r_zero_stack)]
                                            [is_dqm_r_zero_stack][i]['ALPHA_J2000']) + ' ' +  # RA
                            "{:12.7f}".format(sex_cat_r_stack[galaxy_r_stack][0:len(is_dqm_r_zero_stack)]
                                            [is_dqm_r_zero_stack][i]['DELTA_J2000']) + ' ' +  # DEC
                            "{:7.3f}".format(umag) + ' ' +  # u
                            "{:7.3f}".format(uerr) + ' ' +  # uerr
                            "{:7.3f}".format(gmag) + ' ' +  # g
                            "{:7.3f}".format(gerr) + ' ' +  # gerr
                            "{:7.3f}".format(sex_cat_r_stack[galaxy_r_stack][0:len(is_dqm_r_zero_stack)]
                                           [is_dqm_r_zero_stack][i]['MAG_AUTO']) + ' ' +  # r
                            "{:7.3f}".format(sex_cat_r_stack[galaxy_r_stack][0:len(is_dqm_r_zero_stack)]
                                           [is_dqm_r_zero_stack][i]['MAGERR_AUTO']) + '\n')  # rerr


        for i in range(0, num_all):
            prev_time = time.time()
            if i in idx_u_single[matched_u_single]:  #### for u band
                idx_match_u_single = np.where(idx_u_single[matched_u_single] == ref_r['NUMBER'][i])
                u_auto = sex_cat_u[matched_u_single][ab.mag_sys][idx_match_u_single][0]
                uerr_auto = sex_cat_u[matched_u_single][ab.magerr_sys][idx_match_u_single][0]
                from_where_u = 0
            elif i in idx_u_stack[matched_u_stack]:
                idx_match_u = np.where(idx_u_stack[matched_u_stack] == ref_r['NUMBER'][i])
                u_auto = sex_cat_u_stack[matched_u_stack][ab.mag_sys][idx_match_u][0]
                uerr_auto = sex_cat_u_stack[matched_u_stack][ab.magerr_sys][idx_match_u][0]
                from_where_u = 1
            else:
                u_auto = 999.0
                uerr_auto = 999.0
                from_where_u = 2
            if i in idx_g_stack[matched_g_stack]:  #### for g band
                idx_match_g = np.where(idx_g_stack[matched_g_stack] == ref_r['NUMBER'][i])
                g_auto = sex_cat_g_stack[matched_g_stack][ab.mag_sys][idx_match_g][0]
                gerr_auto = sex_cat_g_stack[matched_g_stack][ab.magerr_sys][idx_match_g][0]
                from_where_g = 0
            elif i in idx_g_single[matched_g_single]:
                idx_match_g_single = np.where(idx_g_single[matched_g_single] == ref_r['NUMBER'][i])
                g_auto = sex_cat_g[matched_g_single][ab.mag_sys][idx_match_g_single][0]
                gerr_auto = sex_cat_g[matched_g_single][ab.magerr_sys][idx_match_g_single][0]
                from_where_g = 1
            else:
                g_auto = 999.0
                gerr_auto = 999.0
                from_where_g = 2

            fh.writelines("{:4}".format(i + 1) + ' ' +  # NUMBER
                          "{:12.7f}".format(coords_r[i].ra.value) + ' ' +  # RA
                          "{:12.7f}".format(coords_r[i].dec.value) + ' ' +  # DEC
                          "{:7.3f}".format(ref_r['A_WORLD'][i] * 3600) + ' ' +  # A_WORLD
                          "{:7.3f}".format(ref_r['B_WORLD'][i] * 3600) + ' ' +  # B_WORLD
                          "{:4.2f}".format(ref_r['CLASS_STAR'][i]) + ' ' +  # CLASS_STAR
                          "{:7.3f}".format(u_auto) + ' ' +  # u-band 'MAG_AUTO'
                          "{:7.3f}".format(uerr_auto) + ' ' +  # u-band 'MAGERR_AUTO'
                          "{:7.3f}".format(g_auto) + ' ' +  # g-band 'MAG_AUTO'
                          "{:7.3f}".format(gerr_auto) + ' ' +  # g-band 'MAGERR_AUTO'
                          "{:7.3f}".format(ref_r['MAG_AUTO'][i]) + ' ' +  # r-band 'MAG_AUTO'
                          "{:7.3f}".format(ref_r['MAGERR_AUTO'][i]) + ' ' +  # r-band 'MAGERR_AUTO'
                          "{:2}".format(from_where_u) + ' ' +  # u-band 'FROM_WHERE'
                          "{:2}".format(from_where_g) + ' ' +  # g-band 'FROM_WHERE'
                          "{:2}".format(int(ref_r['FROM_WHERE'][i])) + '\n')  # r-band 'FROM_WHERE'

            color = 'red' if int(ref_r['FROM_WHERE'][i]) else 'green'   # red for stack green for short
            # color = 'darkorange'

            # regFile.writelines("j2000; ellipse({:12.7f}, {:12.7f}, {:7.3f}\", {:7.3f}\", {:7.3f}) #  "
            #                    "color={} \n".format(
            #     coords_r[i].ra.value,
            #     coords_r[i].dec.value,
            #     ref_r['A_WORLD'][i] * 3600 * ref_r['KRON_RADIUS'][i],
            #     ref_r['B_WORLD'][i] * 3600 * ref_r['KRON_RADIUS'][i],
            #     180 - ref_r['THETA_WORLD'][i],
            #     color))

            regFile.writelines("j2000; ellipse({:12.7f}, {:12.7f}, {:7.3f}\", {:7.3f}\", {:7.3f}) # text=\'{}\' "
                               "color={} \n".format(
                coords_r[i].ra.value,
                coords_r[i].dec.value,
                ref_r['A_WORLD'][i] * 3600 * ref_r['KRON_RADIUS'][i],
                ref_r['B_WORLD'][i] * 3600 * ref_r['KRON_RADIUS'][i],
                180 - ref_r['THETA_WORLD'][i],
                i,
                color))

            if i % 10000 == 0:
                print("--- %06.2f minutes ---" % (((time.time() - prev_time)) / 60.0))
                print("--- %d / %d (write cat) ---" % (i, num_all))
                my_module.print_time()
                prev_time = time.time()

    print(f"--- {ab.clusters[k]} Done ---")
    my_module.print_time()

