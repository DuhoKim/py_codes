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
import importlib
importlib.reload(ab)

def dqm_check(n):
    prev_time = start_time
    is_dqm_zero = np.ones(ni, dtype=bool)
    hdu_r_single = fits.open(ab.work_dir + 'fits/best_single/' + ab.short_dqm_fn[k][2])
    hdu_r_stack = fits.open(ab.work_dir + 'fits/stacked/' + ab.stack_dqm_fn[k][2])

    for ii in range(0, ni):
        # read the celestial coordinate
        cel_coord = [[merged_r[bright_central_galaxies]['ALPHA_J2000'][ii+n*ni],
                      merged_r[bright_central_galaxies]['DELTA_J2000'][ii+n*ni]], [0, 0]]

        for jj in range(1, len(hdu_r_single)):
            # read WCS
            w = wcs.WCS(hdu_r_single[jj].header)
            pixcrd = w.wcs_world2pix(cel_coord, 1)
            # if w.footprint_contains(sky_coord):
            if (pixcrd[0][0] > ab.dqm_boundary) & (pixcrd[0][0] < hdu_r_single[jj].shape[1] - ab.dqm_boundary) & \
                    (pixcrd[0][1] > ab.dqm_boundary) & (pixcrd[0][1] < hdu_r_single[jj].shape[0] - ab.dqm_boundary):
                dqm_fits = hdu_r_single[jj].data
                # check the value of DQM
                if k == 0:      # Somehow A754 has value 7 (cosmic ray) in galaxies
                    if dqm_fits[int(pixcrd[0][1])][int(pixcrd[0][0])] % 7:
                        # toggle on saturated bool array
                        is_dqm_zero[ii] = False
                else:           # somehow other clusters have value 128 (?) in galaxies
                    if dqm_fits[int(pixcrd[0][1])][int(pixcrd[0][0])] % 128:
                        # toggle on saturated bool array
                        is_dqm_zero[ii] = False

        for jj in range(1, len(hdu_r_stack)):
            # read WCS
            w = wcs.WCS(hdu_r_stack[jj].header)
            pixcrd = w.wcs_world2pix(cel_coord, 1)
            # if w.footprint_contains(sky_coord):
            if (pixcrd[0][0] > 0) & (pixcrd[0][0] < hdu_r_stack[jj].shape[1]) & \
                    (pixcrd[0][1] > 0) & (pixcrd[0][1] < hdu_r_stack[jj].shape[0]):
                dqm_fits = hdu_r_stack[jj].data
                # check the value of DQM
                if dqm_fits[int(pixcrd[0][1])][int(pixcrd[0][0])]:
                    # toggle on saturated bool array
                    is_dqm_zero[ii] = False

        if ii % 1000 == 0:
            print("--- %06.2f minutes ---" % ((time.time() - prev_time) / 60.0))
            print("--- %d / %d - %d (DQM check) ---" % (ii, n*ni, (n+1)*ni))
            prev_time = time.time()

    hdu_r_single.close()
    hdu_r_stack.close()
    return is_dqm_zero

for k in range(0, len(ab.clusters)):
# for k in range(6, 7):
    print(f"--- {ab.clusters[k]} {ab.ver} initiate ---")
    my_module.print_time()

    start_time = time.time()
    fn = ab.sex_dir+f'DECam_merged_SEx_cat_{ab.clusters[k]}_Gal_ext_corrected_{ab.ver}'

    psf_sex_dir = ab.work_dir + f'sex/run/{ab.clusters[k]}_rsi_'

    with open(fn + '.txt', 'w') as fh, open(fn + '_kron_w_text.reg', 'w') as regFile, open(fn + '_kron.reg', 'w') as regFile2, open(fn + '_kron16.reg', 'w') as regFile16, open(fn + '_kron18.reg', 'w') as regFile18:
        regFile.writelines('# Region file format: DS9 version 4.1\n')
        regFile.writelines('global color=green dashlist=8 3 width=2 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
        regFile.writelines('fk5\n')

        regFile2.writelines('# Region file format: DS9 version 4.1\n')
        regFile2.writelines('global color=green dashlist=8 3 width=2 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
        regFile2.writelines('fk5\n')

        regFile16.writelines('# Region file format: DS9 version 4.1\n')
        regFile16.writelines(
            'global color=green dashlist=8 3 width=2 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
        regFile16.writelines('fk5\n')

        regFile18.writelines('# Region file format: DS9 version 4.1\n')
        regFile18.writelines(
            'global color=green dashlist=8 3 width=2 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
        regFile18.writelines('fk5\n')

        fh.writelines('############################################################################################\n')
        fh.writelines('#   This is a catalog of galaxy sources detected primarily in short-exposure and addition- #\n')
        fh.writelines(f'# ally from stacked Abell cluster {ab.clusters[k]} mosaics in u, g, and r band taken by Dark Energy   #\n')
        fh.writelines('# Camera (DECam) mounted at the prime focus of the Blanco 4-m telescope at CTIO. All FITS  #\n')
        fh.writelines('# files are available on NOAO archive.                                                     #\n')
        fh.writelines('#   Short exposures that I used for photometry are:                                        #\n')
        fh.writelines(f'#   {ab.clusters[k]} u-band : {str(ab.short_exp_times[k][0])}s exposure taken in {ab.short_obs_datas[k][0]}                                       #\n')
        fh.writelines(f'#   {ab.clusters[k]} g-band : {str(ab.short_exp_times[k][1])}s exposure taken in {ab.short_obs_datas[k][1]}                                       #\n')
        fh.writelines(f'#   {ab.clusters[k]} r-band : {str(ab.short_exp_times[k][2])}s exposure taken in {ab.short_obs_datas[k][2]}.                                      #\n')
        fh.writelines('#                                                                                          #\n')
        fh.writelines(f'#   Stacked mosaics with {str(ab.stack_exp_times[k][0])}s, {str(ab.stack_exp_times[k][1])}s, and {str(ab.stack_exp_times[k][2])}s exposure time for u, g, and r, respectiv- #\n')
        fh.writelines('# ely.                                                                                     #\n')
        fh.writelines('#                           ***  Standard Star Photometry  ***                             #\n')
        fh.writelines('#   All data are standardized using standard-star observation taken at same night as the s-#\n')
        fh.writelines('# cience observations. I used Source Extractor to measure magnitudes of standard stars and #\n')
        fh.writelines('# galaxies. I used MAG_APER at 14.8" diameter for measuring standard star magnitudes and   #\n')
        fh.writelines('# compared those with Southern Standard Stars or SDSS catalogue, so that I could fit the   #\n')
        fh.writelines('# zeropoints (ZPs) which are corrected for airmass term at each night.                     #\n')
        fh.writelines('#                                                                                          #\n')
        fh.writelines('#                              ***  Galaxy Photometry  ***                                 #\n')
        fh.writelines('#   Using the ZPs, I measured galaxy magnitudes (MAG_AUTO) from the SExtractor. Galactic   #\n')
        fh.writelines('# extinctions were corrected using the values from http://irsa.ipac.caltech.edu.           #\n')
        fh.writelines('#                                                                                          #\n')
        fh.writelines('#                               ***  Sample Selection  ***                                 #\n')
        fh.writelines(f'#  This catalogue include galaxies with SExtractor parameter \'CLASS_STAR\' < {ab.class_star_lim} (close to  #\n')
        fh.writelines(f'# 0: likely to be galaxies, close to 1: likely to be stars), r-band mag < {ab.mag_lim}, galactic #\n')
        fh.writelines('# central stacked DQM (Data Quality Mask; DECam pipeline product) pixel values are zero,   #\n')
        fh.writelines('# and detected in at least two bands to avoid the false detection in the single band.      #\n')
        fh.writelines('# The u-band and g-band magnitudes are added if there were matched sources within 1\" radii #\n')
        fh.writelines('#   I put arbitrary number 999.0 if there were no matched source.                          #\n')
        fh.writelines('#                                                                                          #\n')
        fh.writelines('#    If you have any question, please send me an email.                                    #\n')
        fh.writelines('#    Duho Kim [email:dhkim@kasi.re.kr]                                                     #\n')
        fh.writelines('#                                                                                          #\n')
        fh.writelines(f"#    This Version {ab.ver} was written on {datetime.today().strftime('%Y-%m-%d')}                               #\n")
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
        sex_cat_u = ascii.read(ab.short_cat_dir + ab.short_cat_fn[k][0] + '_deblend.cat')
        sex_cat_g = ascii.read(ab.short_cat_dir + ab.short_cat_fn[k][1] + '_deblend.cat')
        sex_cat_r = ascii.read(ab.short_cat_dir + ab.short_cat_fn[k][2] + '_deblend.cat')

        sex_cat_u_stack = ascii.read(ab.work_dir + 'sex/cat/stack/' + ab.clusters[k] + '_usi.cat')
        sex_cat_g_stack = ascii.read(ab.work_dir + 'sex/cat/stack/' + ab.clusters[k] + '_gsi.cat')
        sex_cat_r_stack = ascii.read(ab.work_dir + 'sex/cat/stack/' + ab.clusters[k] + '_rsi.cat')

        # read SEx+PSF in r cat for CLASS_STAR
        sex_cat1 = ascii.read(psf_sex_dir + '1_psf.cat')
        sex_cat2 = ascii.read(psf_sex_dir + '2_psf.cat')
        sex_cat3 = ascii.read(psf_sex_dir + '3_psf.cat')
        sex_cat4 = ascii.read(psf_sex_dir + '4_psf.cat')
        sex_cat5 = ascii.read(psf_sex_dir + '5_psf.cat')
        sex_cat6 = ascii.read(psf_sex_dir + '6_psf.cat')
        sex_cat7 = ascii.read(psf_sex_dir + '7_psf.cat')
        sex_cat8 = ascii.read(psf_sex_dir + '8_psf.cat')
        if k != 0:
            sex_cat9 = ascii.read(psf_sex_dir + '9_psf.cat')
            spr_sex_cat = vstack([sex_cat1, sex_cat2, sex_cat3, sex_cat4, sex_cat5, sex_cat6, sex_cat7, sex_cat8, sex_cat9])
        else:
            spr_sex_cat = vstack([sex_cat1, sex_cat2, sex_cat3, sex_cat4, sex_cat5, sex_cat6, sex_cat7, sex_cat8])

        # add column for FROM_WHERE parameter
        sex_cat_u.add_column(np.zeros(len(sex_cat_u)), name='FROM_WHERE')
        sex_cat_g.add_column(np.zeros(len(sex_cat_g)), name='FROM_WHERE')
        sex_cat_r.add_column(np.zeros(len(sex_cat_r)), name='FROM_WHERE')

        sex_cat_u_stack.add_column(np.ones(len(sex_cat_u_stack)), name='FROM_WHERE')
        sex_cat_g_stack.add_column(np.ones(len(sex_cat_g_stack)), name='FROM_WHERE')
        sex_cat_r_stack.add_column(np.ones(len(sex_cat_r_stack)), name='FROM_WHERE')

        # standardize magnitudes and Milky Way extinction correction
        sex_cat_u[ab.mag_sys] = sex_cat_u[ab.mag_sys] + ab.short_a[k][0] + ab.short_b[k][0] - ab.mw_ext[k][0]
        sex_cat_g[ab.mag_sys] = sex_cat_g[ab.mag_sys] + ab.short_a[k][1] + ab.short_b[k][1] - ab.mw_ext[k][1]
        sex_cat_r[ab.mag_sys] = sex_cat_r[ab.mag_sys] + ab.short_a[k][2] + ab.short_b[k][2] - ab.mw_ext[k][2]

        # ZP already include MW extinction correction
        sex_cat_u_stack[ab.mag_sys] = sex_cat_u_stack[ab.mag_sys] + ab.stack_a[k][0] #- ab.mw_ext[k][0]
        sex_cat_g_stack[ab.mag_sys] = sex_cat_g_stack[ab.mag_sys] + ab.stack_a[k][1] #- ab.mw_ext[k][1]
        sex_cat_r_stack[ab.mag_sys] = sex_cat_r_stack[ab.mag_sys] + ab.stack_a[k][2] #- ab.mw_ext[k][2]

        # match stack SEx sources to short SEx to avoid duplication
        coords_u_single = SkyCoord(sex_cat_u['ALPHA_J2000'], sex_cat_u['DELTA_J2000'], unit='deg')
        coords_u_stack = SkyCoord(sex_cat_u_stack['ALPHA_J2000'], sex_cat_u_stack['DELTA_J2000'], unit='deg')
        idx_u_stack, d2d, d3d = coords_u_stack.match_to_catalog_sky(coords_u_single)
        matched_u_stack = (d2d.arcsec < ab.max_sep)

        coords_g_single = SkyCoord(sex_cat_g['ALPHA_J2000'], sex_cat_g['DELTA_J2000'], unit='deg')
        coords_g_stack = SkyCoord(sex_cat_g_stack['ALPHA_J2000'], sex_cat_g_stack['DELTA_J2000'], unit='deg')
        idx_g_stack, d2d, d3d = coords_g_stack.match_to_catalog_sky(coords_g_single)
        matched_g_stack = (d2d.arcsec < ab.max_sep)

        coords_r_single = SkyCoord(sex_cat_r['ALPHA_J2000'], sex_cat_r['DELTA_J2000'], unit='deg')
        coords_r_stack = SkyCoord(sex_cat_r_stack['ALPHA_J2000'], sex_cat_r_stack['DELTA_J2000'], unit='deg')
        idx_r_stack, d2d, d3d = coords_r_stack.match_to_catalog_sky(coords_r_single)
        matched_r_stack = (d2d.arcsec < ab.max_sep)

        # merge stack & short /  base is short to include deblended bright sources from short exposure
        merged_u = vstack([sex_cat_u, sex_cat_u_stack[~matched_u_stack]])
        merged_g = vstack([sex_cat_g, sex_cat_g_stack[~matched_g_stack]])
        merged_r = vstack([sex_cat_r, sex_cat_r_stack[~matched_r_stack]])

        ### exclude stars | magnitudes less than mag lim | outer | saturated ###
        coords_merged_r = SkyCoord(merged_r['ALPHA_J2000'], merged_r['DELTA_J2000'], unit='deg')
        coords_psf = SkyCoord(spr_sex_cat['ALPHA_J2000'], spr_sex_cat['DELTA_J2000'], unit='deg')
        idx_psf2merged, d2d, d3d = coords_merged_r.match_to_catalog_sky(coords_psf)
        merged_r['CLASS_STAR'] = spr_sex_cat[idx_psf2merged]['CLASS_STAR']
        # matched_merged2psf = (d2d.arcsec < ab.max_sep)

        bright_central_galaxies = (merged_r['CLASS_STAR'] < ab.class_star_lim) & (merged_r[ab.mag_sys] < ab.mag_lim)

        num_process = 12
        num_gal_r = sum(bright_central_galaxies)
        ni = int(num_gal_r / num_process)
        with Pool(num_process) as p:
            is_dqm_r_zero = p.map(dqm_check, range(num_process))
        is_dqm_r_zero = np.array(is_dqm_r_zero).reshape(-1)

        # finalize reference r catalogue
        ref_r = merged_r[bright_central_galaxies][0:len(is_dqm_r_zero)][is_dqm_r_zero]
        ref_r.sort(ab.mag_sys)
        coords_r = SkyCoord(ref_r['ALPHA_J2000'], ref_r['DELTA_J2000'], unit='deg')

        ### match r to u and g band ###
        coords_u = SkyCoord(merged_u['ALPHA_J2000'], merged_u['DELTA_J2000'], unit='deg')
        idx_u2r, d2d_r2u, d3d = coords_r.match_to_catalog_sky(coords_u)
        match_r2u = (d2d_r2u.arcsec < ab.max_sep)

        coords_g = SkyCoord(merged_g['ALPHA_J2000'], merged_g['DELTA_J2000'], unit='deg')
        idx_g2r, d2d_r2g, d3d = coords_r.match_to_catalog_sky(coords_g)
        match_r2g = (d2d_r2g.arcsec < ab.max_sep)

        j = 0   # ID number

        for i in range(0, len(ref_r)):
            prev_time = time.time()
            if match_r2u[i]:  #### for u band
                umag = merged_u[idx_u2r][i][ab.mag_sys]
                uerr = merged_u[idx_u2r][i][ab.magerr_sys]
                from_where_u = merged_u[idx_u2r][i]['FROM_WHERE']
            else:
                umag = 999.0
                uerr = 999.0
                from_where_u = 2
            if match_r2g[i]:  #### for g band
                gmag = merged_g[idx_g2r][i][ab.mag_sys]
                gerr = merged_g[idx_g2r][i][ab.magerr_sys]
                from_where_g = merged_g[idx_g2r][i]['FROM_WHERE']
            else:
                gmag = 999.0
                gerr = 999.0
                from_where_g = 2

            if from_where_u == 2 and from_where_g == 2:
                continue
            else:
                j += 1

            fh.writelines("{:4}".format(j) + ' ' +  # NUMBER
                          "{:12.7f}".format(coords_r[i].ra.value) + ' ' +  # RA
                          "{:12.7f}".format(coords_r[i].dec.value) + ' ' +  # DEC
                          "{:7.3f}".format(ref_r['A_WORLD'][i] * 3600) + ' ' +  # A_WORLD
                          "{:7.3f}".format(ref_r['B_WORLD'][i] * 3600) + ' ' +  # B_WORLD
                          "{:4.2f}".format(ref_r['CLASS_STAR'][i]) + ' ' +  # CLASS_STAR
                          "{:7.3f}".format(umag) + ' ' +  # u-band 'MAG_AUTO'
                          "{:7.3f}".format(uerr) + ' ' +  # u-band 'MAGERR_AUTO'
                          "{:7.3f}".format(gmag) + ' ' +  # g-band 'MAG_AUTO'
                          "{:7.3f}".format(gerr) + ' ' +  # g-band 'MAGERR_AUTO'
                          "{:7.3f}".format(ref_r['MAG_AUTO'][i]) + ' ' +  # r-band 'MAG_AUTO'
                          "{:7.3f}".format(ref_r['MAGERR_AUTO'][i]) + ' ' +  # r-band 'MAGERR_AUTO'
                          "{:2}".format(int(from_where_u)) + ' ' +  # u-band 'FROM_WHERE'
                          "{:2}".format(int(from_where_g)) + ' ' +  # g-band 'FROM_WHERE'
                          "{:2}".format(int(ref_r['FROM_WHERE'][i])) + '\n')  # r-band 'FROM_WHERE'

            color = 'red' if int(ref_r['FROM_WHERE'][i]) else 'green'   # red for stack green for short

            regFile.writelines("j2000; ellipse({:12.7f}, {:12.7f}, {:7.3f}\", {:7.3f}\", {:7.3f}) # text=\'{}\' "
                               "color={} \n".format(
                coords_r[i].ra.value,
                coords_r[i].dec.value,
                ref_r['A_WORLD'][i] * 3600 * ref_r['KRON_RADIUS'][i],
                ref_r['B_WORLD'][i] * 3600 * ref_r['KRON_RADIUS'][i],
                180 - ref_r['THETA_WORLD'][i],
                j,
                color))

            regFile2.writelines("j2000; ellipse({:12.7f}, {:12.7f}, {:7.3f}\", {:7.3f}\", {:7.3f})  \n".format(
                coords_r[i].ra.value,
                coords_r[i].dec.value,
                ref_r['A_WORLD'][i] * 3600 * ref_r['KRON_RADIUS'][i],
                ref_r['B_WORLD'][i] * 3600 * ref_r['KRON_RADIUS'][i],
                180 - ref_r['THETA_WORLD'][i]))

            if ref_r['MAG_AUTO'][i] < 16:
                regFile16.writelines("j2000; ellipse({:12.7f}, {:12.7f}, {:7.3f}\", {:7.3f}\", {:7.3f}) \n".format(
                    coords_r[i].ra.value,
                    coords_r[i].dec.value,
                    ref_r['A_WORLD'][i] * 3600 * ref_r['KRON_RADIUS'][i],
                    ref_r['B_WORLD'][i] * 3600 * ref_r['KRON_RADIUS'][i],
                    180 - ref_r['THETA_WORLD'][i]))

            if ref_r['MAG_AUTO'][i] < 18:
                regFile18.writelines("j2000; ellipse({:12.7f}, {:12.7f}, {:7.3f}\", {:7.3f}\", {:7.3f}) \n".format(
                    coords_r[i].ra.value,
                    coords_r[i].dec.value,
                    ref_r['A_WORLD'][i] * 3600 * ref_r['KRON_RADIUS'][i],
                    ref_r['B_WORLD'][i] * 3600 * ref_r['KRON_RADIUS'][i],
                    180 - ref_r['THETA_WORLD'][i]))

            if i % 10000 == 0:
                print("--- %06.2f minutes ---" % (((time.time() - prev_time)) / 60.0))
                print("--- %d / %d (write cat) ---" % (i, len(ref_r)))
                my_module.print_time()
                prev_time = time.time()

    print(f"--- {ab.clusters[k]} Done ---")
    my_module.print_time()

