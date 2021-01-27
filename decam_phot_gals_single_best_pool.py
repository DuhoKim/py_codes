from astropy.io import ascii
import numpy as np
import time
from astropy.coordinates import SkyCoord
import my_module
import pandas as pd
from astropy.io import fits
from astropy import wcs
from datetime import datetime
from multiprocessing import Pool

work_dir=("/Users/duhokim/work/abell/")
server_dir=("/Volumes/APPLE SSD/Server/work/sex/best_single/")

ver = '2013-2014_bright'
class_star_lim = 0.5
fwhm_lim = 0.0
mag_lim = [20, 20, 20, 20, 20, 20, 20]

max_sep = 1.0   # matching radius limit among each individual exposures
mag_sys = 'MAG_ISO'
magerr_sys = 'MAGERR_ISO'
mag_sys1 = 'MAG_AUTO'
magerr_sys1 = 'MAGERR_AUTO'
ccd_x = 2046
ccd_y = 4094

clusters = ['A754', 'A2399', 'A2670', 'A3558', 'A3574', 'A3659', 'A3716']

single_fn = [['A754_us_300', 'A754_gs_300', 'A754_rs_300'],
            ['A2399_us_300_0236', 'A2399_gs_60_0239', 'A2399_rs_300_0142'],
            ['A2670_us_300_0409', 'A2670_gs_60_0420', 'A2670_rs_300_0351'],
            ['A3558_us_300', 'A3558_gs_300', 'A3558_rs_300'],
            ['A3574_us_300', 'A3574_gs_300', 'A3574_rs_300'],
            ['A3659_us_300', 'A3659_gs_300', 'A3659_rs_300'],
             ['A3716_us_300_2357', 'A3716_gs_300_0030', 'A3716_rs_300_0418']]

single_dqm_fn = [['A754_ud_300_11apr.fits', 'A754_gd_300_11apr.fits', 'A754_rd_300_10apr.fits'],
                ['A2399_ud_300_0236_20aug.fits', 'A2399_gd_60_0239_19aug.fits', 'A2399_rd_300_0142_19aug.fits'],
                ['A2670_ud_300_0409_21aug.fits', 'A2670_gd_60_0420_21aug.fits', 'A2670_rd_300_0351_19aug.fits'],
                ['A3558_ud_300_11apr.fits', 'A3558_gd_300_11apr.fits', 'A3558_rd_300_10apr.fits'],
                ['A3574_ud_300_11apr.fits', 'A3574_gd_300_11apr.fits', 'A3574_rd_300_10apr.fits'],
                ['A3659_ud_300_11apr.fits', 'A3659_gd_300_11apr.fits', 'A3659_rd_300_11apr.fits'],
                ['A3716_ud_300_2357_21aug.fits', 'A3716_gd_300_0030_21aug.fits', 'A3716_rd_300_0418_19aug.fits']]

mag_exp_t_corr_300 = 2.5 * np.log10(300)
mag_exp_t_corr_60 = 2.5 * np.log10(60)

### Standardization parameters      mag = inst_mag [ZP:25 with flux/s] + a + b * AIRMASS
# manually exclude outliers
# single_a = [[-1.78154 + mag_exp_t_corr_300, 0.4244 + mag_exp_t_corr_60, 0.5308 + mag_exp_t_corr_300],
#             [-1.78052 + mag_exp_t_corr_300, 0.2787 + mag_exp_t_corr_60, 0.5308 + mag_exp_t_corr_300],
#             [-1.78052 + mag_exp_t_corr_300, 0.2787 + mag_exp_t_corr_300, 0.5308 + mag_exp_t_corr_300]]
#
# single_b = [[-0.09498 * 1.26,    -0.2086 * 1.26, -0.1018 * 1.53],
#             [-0.12552 * 1.3,     -0.0792 * 1.26, -0.1018 * 1.41],
#             [-0.12552 * 1.42,    -0.0792 * 1.31, -0.1018 * 1.09]]

# statsmodel + default SEx param
single_a = [[-1.534 + mag_exp_t_corr_300, 0.426 + mag_exp_t_corr_300, 0.498 + mag_exp_t_corr_300],
            [-1.6863 + mag_exp_t_corr_300, 0.3501 + mag_exp_t_corr_60, 0.4448 + mag_exp_t_corr_300],
            [-1.7019 + mag_exp_t_corr_300, 0.2680 + mag_exp_t_corr_60, 0.4448 + mag_exp_t_corr_300],
            [-1.534 + mag_exp_t_corr_300, 0.426 + mag_exp_t_corr_300, 0.498 + mag_exp_t_corr_300],
            [-1.534 + mag_exp_t_corr_300, 0.426 + mag_exp_t_corr_300, 0.498 + mag_exp_t_corr_300],
            [-1.534 + mag_exp_t_corr_300, 0.426 + mag_exp_t_corr_300, 0.513 + mag_exp_t_corr_300],
            [-1.7019 + mag_exp_t_corr_300, 0.2680 + mag_exp_t_corr_300, 0.4448 + mag_exp_t_corr_300]]

single_b = [[-0.355 * 1.19,    -0.166 * 1.07, -0.061 * 1.07],
            [-0.1724 * 1.26,    -0.1631 * 1.26, -0.0544 * 1.53],
            [-0.2034 * 1.3,     -0.0846 * 1.26, -0.0544 * 1.41],
            [-0.355 * 1.04,     -0.166 * 1.0, -0.061 * 1.0],
            [-0.355 * 1.32,     -0.166 * 1.13, -0.061 * 1.23],
            [-0.355 * 1.1,     -0.166 * 1.2, -0.073 * 1.09],
            [-0.2034 * 1.42,    -0.0846 * 1.31, -0.0544 * 1.09]]

# irsa.ipac.caltech.edu, S and F (2011)
mw_ext = [  [0.308, 0.240, 0.166],
            [0.159, 0.124, 0.086],
            [0.188, 0.146, 0.101],
            [0.212, 0.165, 0.114],
            [0.247, 0.192, 0.133],
            [0.488, 0.380, 0.263],
            [0.157, 0.122, 0.085]]

am_cat = pd.read_csv(work_dir+'sex/cat/airmass.csv')

def dqm_check_single(n):
    prev_time = start_time
    is_dqm_zero = np.ones(ni, dtype=bool)
    hdu_r = fits.open(work_dir + 'fits/best_single/' + single_dqm_fn[k][2])
    for ii in range(0, ni):
        # read the celestial coordinate
        cel_coord = [[sex_cat_r[criteria_r]['ALPHA_J2000'][ii+n*ni],
                      sex_cat_r[criteria_r]['DELTA_J2000'][ii+n*ni]], [0, 0]]

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


# for k in range(0, len(clusters)):
for k in range(0, 1):
    start_time = time.time()
    fn = work_dir+'sex/cat/DECam_short_exposure_SEx_cat_' + clusters[k] + \
         '_match_rad_1as_Gal_ext_corrected_'+ver

    with open(fn + '.txt', 'w') as fh, open(fn + '.reg', 'w') as regFile:
        fh.writelines('############################################################################################\n')
        fh.writelines('#   This is a catalog of photometry of galaxies in Abell cluster '+clusters[k]+' mosaics in          #\n')
        fh.writelines('# u, g, and r band taken by Dark Energy Camera (DECam) mounted at the prime focus of the B-#\n')
        fh.writelines('# lanco 4-m telescope at CTIO. The single exposures that I used for photometry are:        #\n')
        if k == 0:
            fh.writelines('#   A754  u-band : 300s exposure taken at 02:21 in 11-Apr 2013                             #\n')
            fh.writelines('#   A754  g-band : 300s exposure taken at 00:21 in 11-Apr 2013                             #\n')
            fh.writelines('#   A754  r-band : 300s exposure taken at 01:00 in 10-Apr 2013                             #\n')
        if k == 1:
            fh.writelines('#   A2399 u-band : 300s exposure taken at 02:36 in 20-Aug 2014                             #\n')
            fh.writelines('#   A2399 g-band :  60s exposure taken at 02:39 in 19-Aug 2014                             #\n')
            fh.writelines('#   A2399 r-band : 300s exposure taken at 01:42 in 19-Aug 2014                             #\n')
        if k == 2:
            fh.writelines('#   A2670 u-band : 300s exposure taken at 04:09 in 21-Aug 2014                             #\n')
            fh.writelines('#   A2670 g-band :  60s exposure taken at 04:20 in 21-Aug 2014                             #\n')
            fh.writelines('#   A2670 r-band : 300s exposure taken at 03:51 in 19-Aug 2014                             #\n')
        if k == 3:
            fh.writelines('#   A3558 u-band : 300s exposure taken at 06:06 in 11-Apr 2013                             #\n')
            fh.writelines('#   A3558 g-band : 300s exposure taken at 04:53 in 11-Apr 2013                             #\n')
            fh.writelines('#   A3558 r-band : 300s exposure taken at 04:36 in 10-Apr 2013                             #\n')
        if k == 4:
            fh.writelines('#   A3574 u-band : 300s exposure taken at 08:19 in 11-Apr 2013                             #\n')
            fh.writelines('#   A3574 g-band : 300s exposure taken at 07:19 in 11-Apr 2013                             #\n')
            fh.writelines('#   A3574 r-band : 300s exposure taken at 07:59 in 10-Apr 2013                             #\n')
        if k == 5:
            fh.writelines('#   A3659 u-band : 300s exposure taken at 09:29 in 11-Apr 2013                             #\n')
            fh.writelines('#   A3659 g-band : 300s exposure taken at 08:47 in 11-Apr 2013                             #\n')
            fh.writelines('#   A3659 r-band : 300s exposure taken at 09:35 in 11-Apr 2013                             #\n')
        if k == 6:
            fh.writelines('#   A3716 u-band : 300s exposure taken at 23:57 in 21-Aug 2014                             #\n')
            fh.writelines('#   A3716 g-band : 300s exposure taken at 00:30 in 21-Aug 2014                             #\n')
            fh.writelines('#   A3716 r-band : 300s exposure taken at 04:18 in 19-Aug 2014                             #\n')
        fh.writelines('# , which was selected based on seeing and airmass.                                        #\n')
        fh.writelines('#                                                                                          #\n')
        fh.writelines('#                                  ***  Standardization  ***                               #\n')
        fh.writelines('#   All 3-bands in each observing data are standardized using standard star observations.  #\n')
        fh.writelines('#   All FITS files are available on NOAO archive. I used Source Extractor to measure       #\n')
        fh.writelines('# magnitudes of standard stars and galaxies. I used MAG_APER at 14.8" diameter for standard#\n')
        fh.writelines('# star photometry to convert instrument magnitudes (ZP=25) to AB magnitudes correcting     #\n')
        fh.writelines('# for atmospheric extinction by fitting the airmass term.                                  #\n')
        fh.writelines('#                                                                                          #\n')
        fh.writelines('#                                    ***  Photometry  ***                                  #\n')
        fh.writelines('#   The MAG_AUTO and MAG_ISO of galaxies from the Source Extractor are added. I calibrated #\n')
        fh.writelines('# both magnitudes using the standardization equation that I derived above. The Galactic ex-#\n')
        fh.writelines('# tinction was corrected adopting the values from irsa.ipac.caltech.edu, S and F (2011).   #\n')
        fh.writelines('#                                                                                          #\n')
        fh.writelines('#                                  ***  Sample Selection  ***                              #\n')
        fh.writelines('#  I added galaxies into this catalog if they have MAG_AUTO < {} in r band, values of the  #\n'
                      .format(mag_lim[2]))
        fh.writelines('# Source Extractor parameter \'CLASS_STAR\' < {:4.2f} (close to 0: likely to be galaxies,       #\n'
                      .format(class_star_lim))
        fh.writelines('# close to 1: likely to be stars), and central pixel values of \'DQM\' (Data Quality Mask;   #\n')
        fh.writelines('# DECam pipeline product which was useful for getting rid of bright stars) == 0 or 128. The#\n')
        fh.writelines('# u-band and g-band magnitudes are added if there were matched sources within 1\" radius w- #\n')
        fh.writelines('# ith their MAG_AUTO < {} and {}, respectively. I put arbitrary numb-er 999.0 if there were#\n'
                      .format(mag_lim[0], mag_lim[1]))
        fh.writelines('# no matched source.                                                                       #\n')
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
        # sex_cat_u = ascii.read(server_dir + single_fn[k][0] + '_deblend.cat')
        # sex_cat_g = ascii.read(server_dir + single_fn[k][1] + '_deblend.cat')
        # sex_cat_r = ascii.read(server_dir + single_fn[k][2] + '_deblend.cat')
        sex_cat_u = ascii.read(server_dir + single_fn[k][0] + '_default.cat')
        sex_cat_g = ascii.read(server_dir + single_fn[k][1] + '_default.cat')
        sex_cat_r = ascii.read(server_dir + single_fn[k][2] + '_default.cat')

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

        num_process = 12
        num_gal_r = len(sex_cat_r[criteria_r])
        ni = int(num_gal_r / num_process)
        with Pool(num_process) as p:
            is_dqm_r_zero = p.map(dqm_check_single, range(num_process))
        is_dqm_r_zero = np.array(is_dqm_r_zero).reshape(-1)

        # r-band is the reference
        num_of_gal = len(crit_r[0:len(is_dqm_r_zero)][is_dqm_r_zero])

        coords_r = SkyCoord(crit_r['ALPHA_J2000'][0:len(is_dqm_r_zero)][is_dqm_r_zero], crit_r['DELTA_J2000'][0:len(is_dqm_r_zero)][is_dqm_r_zero], unit='deg')
        coords_u = SkyCoord(sex_cat_u['ALPHA_J2000'], sex_cat_u['DELTA_J2000'], unit='deg')

        idx_u, d2d, d3d = coords_u.match_to_catalog_sky(coords_r)
        sep_constraint_u = (d2d.arcsec < max_sep) & (sex_cat_u[mag_sys1] < mag_lim[0])
        u_match_to_ref = sex_cat_u[sep_constraint_u]
        r_match_to_u = crit_r[0:len(is_dqm_r_zero)][is_dqm_r_zero][idx_u[sep_constraint_u]]

        coords_g = SkyCoord(sex_cat_g['ALPHA_J2000'], sex_cat_g['DELTA_J2000'], unit='deg')

        idx_g, d2d, d3d = coords_g.match_to_catalog_sky(coords_r)
        sep_constraint_g = (d2d.arcsec < max_sep) & (sex_cat_g[mag_sys1] < mag_lim[1])
        g_match_to_ref = sex_cat_g[sep_constraint_g]
        r_match_to_g = crit_r[0:len(is_dqm_r_zero)][is_dqm_r_zero][idx_g[sep_constraint_g]]

        for i in range(0, num_of_gal):
            if i in idx_u[sep_constraint_u]:
                idx_match_u = np.where(r_match_to_u['NUMBER'] == crit_r['NUMBER'][0:len(is_dqm_r_zero)][is_dqm_r_zero][i])
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
                idx_match_g = np.where(r_match_to_g['NUMBER'] == crit_r['NUMBER'][0:len(is_dqm_r_zero)][is_dqm_r_zero][i])
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
                "{:7.3f}".format(crit_r['A_WORLD'][0:len(is_dqm_r_zero)][is_dqm_r_zero][i]*3600) + ' ' +      # A_WORLD
                "{:7.3f}".format(crit_r['B_WORLD'][0:len(is_dqm_r_zero)][is_dqm_r_zero][i]*3600) + ' ' +      # B_WORLD
                "{:4.2f}".format(crit_r['CLASS_STAR'][0:len(is_dqm_r_zero)][is_dqm_r_zero][i]) + ' ' +        # CLASS_STAR
                "{:7.3f}".format(u_auto) + ' ' +                                        # u-band 'MAG_AUTO'
                "{:7.3f}".format(uerr_auto) + ' ' +                                     # u-band 'MAGERR_AUTO'
                "{:7.3f}".format(u_iso) + ' ' +                                         # u-band 'MAG_ISO'
                "{:7.3f}".format(uerr_iso) + ' ' +                                      # u-band 'MAGERR_ISO'
                "{:7.3f}".format(g_auto) + ' ' +                                        # g-band 'MAG_AUTO'
                "{:7.3f}".format(gerr_auto) + ' ' +                                     # g-band 'MAGERR_AUTO'
                "{:7.3f}".format(g_iso) + ' ' +                                         # g-band 'MAG_ISO'
                "{:7.3f}".format(gerr_iso) + ' ' +                                      # g-band 'MAGERR_ISO'
                "{:7.3f}".format(crit_r[mag_sys1][0:len(is_dqm_r_zero)][is_dqm_r_zero][i]) + ' ' +            # r-band 'MAG_AUTO'
                "{:7.3f}".format(crit_r[magerr_sys1][0:len(is_dqm_r_zero)][is_dqm_r_zero][i]) + ' ' +         # r-band 'MAGERR_AUTO'
                "{:7.3f}".format(crit_r[mag_sys][0:len(is_dqm_r_zero)][is_dqm_r_zero][i]) + ' ' +             # r-band 'MAG_ISO'
                "{:7.3f}".format(crit_r[magerr_sys][0:len(is_dqm_r_zero)][is_dqm_r_zero][i]) + '\n')         # r-band 'MAGERR_ISO'
            if crit_r['FLAGS'][0:len(is_dqm_r_zero)][is_dqm_r_zero][i] > 0:
                color = 'red'
            else:
                color = 'green'
            # regFile.writelines("j2000; ellipse({:12.7f}, {:12.7f}, {:7.3f}\", {:7.3f}\", {:7.3f}) # text=\'{}\', "
            #                    "color={} \n".format(
            #     coords_r[i].ra.value,
            #     coords_r[i].dec.value,
            #     crit_r['A_WORLD'][0:len(is_dqm_r_zero)][is_dqm_r_zero][i]*3600,
            #     crit_r['B_WORLD'][0:len(is_dqm_r_zero)][is_dqm_r_zero][i]*3600,
            #     180-crit_r['THETA_WORLD'][0:len(is_dqm_r_zero)][is_dqm_r_zero][i],
            #     i+1,
            #     color))
            regFile.writelines("j2000; ellipse({:12.7f}, {:12.7f}, {:7.3f}\", {:7.3f}\", {:7.3f}) #  "
                               "color={} \n".format(
                coords_r[i].ra.value,
                coords_r[i].dec.value,
                crit_r['A_WORLD'][0:len(is_dqm_r_zero)][is_dqm_r_zero][i]*3600,
                crit_r['B_WORLD'][0:len(is_dqm_r_zero)][is_dqm_r_zero][i]*3600,
                180-crit_r['THETA_WORLD'][0:len(is_dqm_r_zero)][is_dqm_r_zero][i],
                color))

    print("--- %s minutes ---" % (((time.time() - start_time))/60.0))

    my_module.print_time()