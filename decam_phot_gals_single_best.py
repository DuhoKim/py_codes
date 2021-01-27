from astropy.io import ascii
import numpy as np
import time
from astropy.coordinates import SkyCoord
import my_module
import pandas as pd
from astropy.io import fits
from astropy import wcs
from datetime import datetime

work_dir=("/Users/duhokim/work/abell/")

# v1.2
# class_star_lim = 0.5
# fwhm_lim = 6.0    # 1.578"
# flag_lim = 999    # no neighbor no blend no saturate

# v1.3
# class_star_lim = 0.3
# fwhm_lim = 5.0    # 1.315"
# flag_lim = 999    # no neighbor no blend no saturate

# v1.4
# class_star_lim = 0.5
# fwhm_lim = 10.0
# flag_lim = 999    # no neighbor no blend no saturate

# v1.5
# class_star_lim = 0.5
# fwhm_lim = 7.0      # 1.841"
# flag_lim = 999    # no neighbor no blend no saturate

# v1.6
# ver = 'v1.6'
# class_star_lim = 0.5
# fwhm_lim = 5.0
# flag_lim = 999    # no neighbor no blend no saturate

# ver = 'v1.7'
# class_star_lim = 0.7
# fwhm_lim = 3.0      # 1.052"
# flag_lim = 999    # no neighbor no blend no saturate

# ver = 'v1.8'
# class_star_lim = 0.9
# fwhm_lim = 0.0      # 1.052"
# mag_lim = [21, 22, 22]

# ver = 'v1.0'
# class_star_lim = 0.9
# fwhm_lim = 0.0
# mag_lim = [21, 22, 22]

# ver = 'v1.1'
# class_star_lim = 0.9
# fwhm_lim = 0.0
# mag_lim = [99, 99, 99]

ver = 'deblended'
class_star_lim = 0.5
fwhm_lim = 0.0
mag_lim = [99, 99, 99]

max_sep = 1.0   # matching radius limit among each individual exposures
mag_sys = 'MAG_ISO'
magerr_sys = 'MAGERR_ISO'
mag_sys1 = 'MAG_AUTO'
magerr_sys1 = 'MAGERR_AUTO'
ccd_x = 2046
ccd_y = 4094

clusters = ['A2399', 'A2670', 'A3716']
# single_fn = [['A2399_us_300_0236_20aug.cat', 'A2399_gs_60_0239_19aug.cat', 'A2399_rs_300_0142_19aug.cat'],
#              ['A2670_us_300_0409_21aug.cat', 'A2670_gs_60_0420_21aug.cat', 'A2670_rs_300_0351_19aug.cat'],
#              ['A3716_us_300_2357_21aug.cat', 'A3716_gs_300_0030_21aug.cat', 'A3716_rs_300_0418_19aug.cat']]
single_fn = [['A2399_us_300_0236_deblend.cat', 'A2399_gs_60_0239_deblend.cat', 'A2399_rs_300_0142_deblend.cat'],
             ['A2670_us_300_0409_deblend.cat', 'A2670_gs_60_0420_deblend.cat', 'A2670_rs_300_0351_deblend.cat'],
             ['A3716_us_300_2357_deblend.cat', 'A3716_gs_300_0030_deblend.cat', 'A3716_rs_300_0418_deblend.cat']]

single_dqm_fn = [['A2399_ud_300_0236_20aug.fits', 'A2399_gd_60_0239_19aug.fits', 'A2399_rd_300_0142_19aug.fits'],
             ['A2670_ud_300_0409_21aug.fits', 'A2670_gd_60_0420_21aug.fits', 'A2670_rd_300_0351_19aug.fits'],
             ['A3716_ud_300_2357_21aug.fits', 'A3716_gd_300_0030_21aug.fits', 'A3716_rd_300_0418_19aug.fits']]

mag_exp_t_corr_300 = 2.5 * np.log10(300)
mag_exp_t_corr_60 = 2.5 * np.log10(60)

### Standardization parameters      mag = inst_mag [ZP:25 with flux/s] + a + b * AIRMASS
# single_a = [[-1.257 + mag_exp_t_corr_300, 0.381 + mag_exp_t_corr_60, 0.567 + mag_exp_t_corr_300],
#             [-1.796 + mag_exp_t_corr_300, 0.179 + mag_exp_t_corr_60, 0.567 + mag_exp_t_corr_300],
#             [-1.796 + mag_exp_t_corr_300, 0.179 + mag_exp_t_corr_300, 0.567 + mag_exp_t_corr_300]]
#
# single_b = [[-0.592 * 1.26, -0.186 * 1.26, -0.136 * 1.53],
#             [-0.132 * 1.3, -0.036 * 1.26, -0.136 * 1.41],
#             [-0.132 * 1.42, -0.036 * 1.31, -0.136 * 1.09]]

single_a = [[-1.9467 + mag_exp_t_corr_300, 0.432 + mag_exp_t_corr_60, 0.521 + mag_exp_t_corr_300],
            [-1.8171 + mag_exp_t_corr_300, 0.192 + mag_exp_t_corr_60, 0.521 + mag_exp_t_corr_300],
            [-1.8171 + mag_exp_t_corr_300, 0.192 + mag_exp_t_corr_300, 0.521 + mag_exp_t_corr_300]]

single_b = [[0.04658 * 1.26,    -0.215 * 1.26, -0.097 * 1.53],
            [-0.0951 * 1.3,     -0.032 * 1.26, -0.097 * 1.41],
            [-0.0951 * 1.42,    -0.032 * 1.31, -0.097 * 1.09]]

# irsa.ipac.caltech.edu, S and F (2011)
mw_ext = [[0.159, 0.124, 0.086],
               [0.188, 0.146, 0.101],
               [0.157, 0.122, 0.085]]

am_cat = pd.read_csv(work_dir+'sex/cat/airmass.csv')

for k in range(0, len(clusters)):
#for k in range(0, 1):
    start_time = time.time()
    fn = work_dir+'sex/cat/DECam_19_21_aug_2014_single_best_exposure_SEx_cat_' + clusters[k] + \
         '_match_rad_1as_Gal_ext_corrected_'+ver

    with open(fn + '.txt', 'w') as fh, open(fn + '.reg', 'w') as regFile:
        fh.writelines('############################################################################################\n')
        fh.writelines('#   This is a catalog of photometry of galaxies in Abell cluster '+clusters[k]+' mosaics in          #\n')
        fh.writelines('# u, g, and r band taken by Dark Energy Camera (DECam) mounted at the prime focus of the B-#\n')
        fh.writelines('# lanco 4-m telescope at CTIO over 19th - 21st Aug 2014. The single exposures that I used  #\n')
        fh.writelines('# for photometry are:                                                                      #\n')
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
        fh.writelines('# , which was selected based on seeing and airmass.                                        #\n')
        fh.writelines('#                                                                                          #\n')
        fh.writelines('#                                  ***  Standardization  ***                               #\n')
        fh.writelines('#   All 3-bands in 3-days of observing data are standardized using standard star observati-#\n')
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
        fh.writelines('#  I added galaxies into this catalog if they have MAG_AUTO < {} in r band, values of the  #\n'
                      .format(mag_lim[2]))
        fh.writelines('# Source Extractor parameter \'CLASS_STAR\' < {:4.2f} (close to 0: likely to be galaxies,  #\n'
                      .format(class_star_lim))
        fh.writelines('# close to 1: likely to be stars), and central pixel values of \'DQM\' (Data Quality Mask; #\n')
        fh.writelines('# DECam pipeline product which was useful for getting rid of bright stars) == 0 or 128. The#\n')
        fh.writelines('# u-band and g-band magnitudes are added if there were matched sources within 1\" radius w-#\n')
        fh.writelines('# ith their MAG_AUTO < {} and {}, respectively. I put arbitrary numb-er 999.0 if there were#\n'
                      .format(mag_lim[0], mag_lim[1]))
        fh.writelines('# no matched source.                                                                                           #\n')
        fh.writelines('#    If you have any question, please send me an email.                                    #\n')
        fh.writelines('#    Duho Kim [email:dhkim@kasi.re.kr]                                                     #\n')
        fh.writelines('#                                                                                          #\n')
        fh.writelines('#    This Version 1.0 was written on {}                                            #\n'
                      .format(datetime.today().strftime('%d-%m-%Y')))
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
        fh.writelines('############################################################################################\n')

        # read SEx cat
        sex_cat_u = ascii.read(work_dir + 'sex/cat/best_single/' + single_fn[k][0])
        sex_cat_g = ascii.read(work_dir + 'sex/cat/best_single/' + single_fn[k][1])
        sex_cat_r = ascii.read(work_dir + 'sex/cat/best_single/' + single_fn[k][2])

        # standardize magnitudes and Milky Way extinction correction
        sex_cat_u[mag_sys] = sex_cat_u[mag_sys] + single_a[k][0] + single_b[k][0] - mw_ext[k][0]
        sex_cat_g[mag_sys] = sex_cat_g[mag_sys] + single_a[k][1] + single_b[k][1] - mw_ext[k][1]
        sex_cat_r[mag_sys] = sex_cat_r[mag_sys] + single_a[k][2] + single_b[k][2] - mw_ext[k][2]

        sex_cat_u[mag_sys1] = sex_cat_u[mag_sys1] + single_a[k][0] + single_b[k][0] - mw_ext[k][0]
        sex_cat_g[mag_sys1] = sex_cat_g[mag_sys1] + single_a[k][1] + single_b[k][1] - mw_ext[k][1]
        sex_cat_r[mag_sys1] = sex_cat_r[mag_sys1] + single_a[k][2] + single_b[k][2] - mw_ext[k][2]

        # for every catalog item
        for i in range(0, len(sex_cat_u)):
            # make 'NUMBER' unique
            sex_cat_u['NUMBER'][i] = i

        for i in range(0, len(sex_cat_g)):
            # make 'NUMBER' unique
            sex_cat_g['NUMBER'][i] = i

        for i in range(0, len(sex_cat_r)):
            # make 'NUMBER' unique
            sex_cat_r['NUMBER'][i] = i

        criteria_r = (sex_cat_r[mag_sys1] < mag_lim[2]) & (sex_cat_r['CLASS_STAR'] < class_star_lim)  & \
                    (sex_cat_r['FWHM_IMAGE'] > fwhm_lim)

        crit_r = sex_cat_r[criteria_r]

        # read Data Quality Mask (DQM) fits file to check saturation
        hdu_r = fits.open(work_dir + 'fits/best_single/' + single_dqm_fn[k][2])

        # make a boolean array to set FLAGS to SATURATED = 4 (could not directly)
        is_dqm_r_zero = np.ones(len(crit_r), dtype=bool)

        crit_r.sort(mag_sys1)

        for i in range(0, len(crit_r)):
            # read the celestial coordinate
            cel_coord = [[crit_r['ALPHA_J2000'][i], crit_r['DELTA_J2000'][i]], [0, 0]]
            sky_coord = SkyCoord(crit_r['ALPHA_J2000'][i], crit_r['DELTA_J2000'][i], unit='deg')
            # for every CCD
            for j in range(1, 61):
                # read WCS
                w = wcs.WCS(hdu_r[j].header)
                pixcrd = w.wcs_world2pix(cel_coord, 1)
                # if w.footprint_contains(sky_coord):
                if (pixcrd[0][0] > 0) & (pixcrd[0][0] < ccd_x) & (pixcrd[0][1] > 0) & (pixcrd[0][1] < ccd_y):
                    dqm_fits = hdu_r[j].data
                    # check the value of DQM
                    if dqm_fits[int(pixcrd[0][1])][int(pixcrd[0][0])] % 128:
                        # toggle on saturated bool array
                        is_dqm_r_zero[i] = False
            if i % 1000 == 0:
                print("--- %s minutes ---" % (((time.time() - start_time)) / 60.0))
                print("--- %d / %d ---" % (i, len(crit_r)))
                my_module.print_time()

        # r-band is the reference
        num_of_gal = len(crit_r[is_dqm_r_zero])

        coords_r = SkyCoord(crit_r['ALPHA_J2000'][is_dqm_r_zero], crit_r['DELTA_J2000'][is_dqm_r_zero], unit='deg')
        coords_u = SkyCoord(sex_cat_u['ALPHA_J2000'], sex_cat_u['DELTA_J2000'], unit='deg')

        idx_u, d2d, d3d = coords_u.match_to_catalog_sky(coords_r)
        sep_constraint_u = (d2d.arcsec < max_sep) & (sex_cat_u[mag_sys1] < mag_lim[0])
        u_match_to_ref = sex_cat_u[sep_constraint_u]
        r_match_to_u = crit_r[is_dqm_r_zero][idx_u[sep_constraint_u]]

        coords_g = SkyCoord(sex_cat_g['ALPHA_J2000'], sex_cat_g['DELTA_J2000'], unit='deg')

        idx_g, d2d, d3d = coords_g.match_to_catalog_sky(coords_r)
        sep_constraint_g = (d2d.arcsec < max_sep) & (sex_cat_g[mag_sys1] < mag_lim[1])
        g_match_to_ref = sex_cat_g[sep_constraint_g]
        r_match_to_g = crit_r[is_dqm_r_zero][idx_g[sep_constraint_g]]

        for i in range(0, num_of_gal):
            if i in idx_u[sep_constraint_u]:
                idx_match_u = np.where(r_match_to_u['NUMBER'] == crit_r['NUMBER'][is_dqm_r_zero][i])
                u_iso = u_match_to_ref[mag_sys][idx_match_u][0]
                uerr_iso = u_match_to_ref[magerr_sys][idx_match_u][0]
                u_auto = u_match_to_ref[mag_sys1][idx_match_u][0]
                uerr_auto = u_match_to_ref[magerr_sys1][idx_match_u][0]
            else:
                u_iso = 999.0
                uerr_iso = 999.0
                u_auto = 999.0
                uerr_auto = 999.0
            if i in idx_g[sep_constraint_g]:
                idx_match_g = np.where(r_match_to_g['NUMBER'] == crit_r['NUMBER'][is_dqm_r_zero][i])
                g_iso = g_match_to_ref[mag_sys][idx_match_g][0]
                gerr_iso = g_match_to_ref[magerr_sys][idx_match_g][0]
                g_auto = g_match_to_ref[mag_sys1][idx_match_g][0]
                gerr_auto = g_match_to_ref[magerr_sys1][idx_match_g][0]
            else:
                g_iso = 999.0
                gerr_iso = 999.0
                g_auto = 999.0
                gerr_auto = 999.0
            fh.writelines("{:4}".format(i+1) + ' ' +                                    # NUMBER
                "{:12.7f}".format(coords_r[i].ra.value)+' '+                            # RA
                "{:12.7f}".format(coords_r[i].dec.value) + ' ' +                        # DEC
                "{:7.3f}".format(crit_r['A_WORLD'][is_dqm_r_zero][i]*3600) + ' ' +      # A_WORLD
                "{:7.3f}".format(crit_r['B_WORLD'][is_dqm_r_zero][i]*3600) + ' ' +      # B_WORLD
                "{:4.2f}".format(crit_r['CLASS_STAR'][is_dqm_r_zero][i]) + ' ' +        # CLASS_STAR
                "{:7.3f}".format(u_auto) + ' ' +                                        # u-band 'MAG_AUTO'
                "{:7.3f}".format(uerr_auto) + ' ' +                                     # u-band 'MAGERR_AUTO'
                "{:7.3f}".format(u_iso) + ' ' +                                         # u-band 'MAG_ISO'
                "{:7.3f}".format(uerr_iso) + ' ' +                                      # u-band 'MAGERR_ISO'
                "{:7.3f}".format(g_auto) + ' ' +                                        # g-band 'MAG_AUTO'
                "{:7.3f}".format(gerr_auto) + ' ' +                                     # g-band 'MAGERR_AUTO'
                "{:7.3f}".format(g_iso) + ' ' +                                         # g-band 'MAG_ISO'
                "{:7.3f}".format(gerr_iso) + ' ' +                                      # g-band 'MAGERR_ISO'
                "{:7.3f}".format(crit_r[mag_sys1][is_dqm_r_zero][i]) + ' ' +            # r-band 'MAG_AUTO'
                "{:7.3f}".format(crit_r[magerr_sys1][is_dqm_r_zero][i]) + ' ' +         # r-band 'MAGERR_AUTO'
                "{:7.3f}".format(crit_r[mag_sys][is_dqm_r_zero][i]) + ' ' +             # r-band 'MAG_ISO'
                "{:7.3f}".format(crit_r[magerr_sys][is_dqm_r_zero][i]) + '\n')         # r-band 'MAGERR_ISO'
            if crit_r['FLAGS'][is_dqm_r_zero][i] > 0:
                color = 'red'
            else:
                color = 'green'
            regFile.writelines("j2000; ellipse({:12.7f}, {:12.7f}, {:7.3f}\", {:7.3f}\", {:7.3f}) # text=\'{}\', "
                               "color={}, dash=1 \n".format(
                coords_r[i].ra.value,
                coords_r[i].dec.value,
                crit_r['A_WORLD'][is_dqm_r_zero][i]*3600,
                crit_r['B_WORLD'][is_dqm_r_zero][i]*3600,
                180-crit_r['THETA_WORLD'][is_dqm_r_zero][i],
                i+1,
                color))

    print("--- %s minutes ---" % (((time.time() - start_time))/60.0))

    my_module.print_time()