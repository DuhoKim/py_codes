from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
import time
import Funcs
from astropy.table import Table, vstack

plot_dir=("/Users/dkim108/Documents/work/plot/")
cat_dir=("/Users/dkim108/Documents/work/cat/SouthernStandardStars/")

dates = ['19aug', '20aug', '21aug']

band = 'u'
exp_t = 100.0
mag_col = 'col4'
err_col = 'col5'

for m in range(2, len(dates)):
    start_time = time.time()
    if m == 0:
        ss_date = '19aug'
        sex_dir=("/Users/dkim108/Documents/work/sex/19aug/")
        ss_name = ['220100', 'E8-A', 'LSE_259', '190000']
        ss_legend_name = ['220100 (1.96)', 'E8-A (1.21)', 'LSE_259 (1.75)', '190000 (1.23)']
        ss_cat_name = ['220100-300000', 'E8_a', 'LSE_259', '190000-295600']
        ss_airmass = [1.96, 1.21, 1.75, 1.23]
        ss_seeing = [1.04, 0.95, 1.03, 0.92]

        max_sep = 1.0
        class_star_lim = 0.97
        flag_lim = 0

        a_guess = -1.641 + np.linspace(-0.1, 0.1, 100)
        b_guess = -0.304 + np.linspace(-0.1, 0.1, 100)

    elif m == 1:
        ss_date = '20aug'
        sex_dir=("/Users/dkim108/Documents/work/sex/20aug/")
        ss_name = ['E5-A', 'E6-A', 'LSE_259_1_11', '220100_2_44', '220100_1_01', 'LSE_259_1_71']
        ss_legend_name = ['E5-A (1.7)', 'E6-A (1.12)', 'LSE_259 (1.11)', '220100 (2.44)', '220100 (1.01)', 'LSE_259 (1.71)']
        ss_cat_name = ['E5_a', 'E6_a', 'LSE_259', '220100-300000', '220100-300000', 'LSE_259']
        ss_airmass = [1.7, 1.12, 1.11, 2.44, 1.01, 1.71]
        ss_seeing = [1.75, 1.49, 1.31, 1.77, 3.5, 3.8]

        max_sep = 3.0
        class_star_lim = 0.5
        flag_lim = 15

        a_guess = -1.257 + np.linspace(-0.01, 0.01, 100)
        b_guess = -0.592 + np.linspace(-0.01, 0.01, 100)

    elif m == 2:
        ss_date = '21aug'
        sex_dir=("/Users/dkim108/Documents/work/sex/21aug/")
        ss_name = ['E5-A', 'E6-A', 'LSE_259_1_11', 'LSE_259_1_8']
        ss_legend_name = ['E5-A (1.65)', 'E6-A (1.11)', 'LSE_259 (1.11)', 'LSE_259 (1.8)']
        ss_cat_name = ['E5_a', 'E6_a', 'LSE_259', 'LSE_259']
        ss_airmass = [1.65, 1.11, 1.11, 1.8]
        ss_seeing = [1.46, 1.03, 1.07, 1.58]

        max_sep = 1.0
        class_star_lim = 0.9
        flag_lim = 3

        a_guess = -1.795 + np.linspace(-0.01, 0.01, 100)
        b_guess = -0.133 + np.linspace(-0.01, 0.01, 100)

    sample_numbers = np.zeros(len(ss_name))
    for i in range(0, len(ss_name)):
        sex_result = ascii.read(sex_dir + ss_name[i] + '_' + band + 'a.cat')
        sex_result_good = (sex_result['MAG_ISO'] < 90) & (sex_result['FLAGS'] <= flag_lim) & \
                          (sex_result['CLASS_STAR'] > class_star_lim)
        sex_coords = SkyCoord(sex_result['ALPHA_J2000'][sex_result_good],
                              sex_result['DELTA_J2000'][sex_result_good],
                              unit='deg')

        ss_cat = ascii.read(cat_dir + ss_cat_name[i] + '.dat.trimmed', format='no_header')
        cat_coords = SkyCoord(ss_cat['col2'], ss_cat['col3'], unit=(u.hourangle, u.deg))

        idx, d2d, d3d = sex_coords.match_to_catalog_sky(cat_coords)
        sep_constraint = (d2d.arcsec < max_sep) & (ss_cat[idx][err_col] > 0)
        sex_matches = sex_result[sex_result_good][sep_constraint]
        cat_matches = ss_cat[idx[sep_constraint]]

        sample_numbers[i] = len(ss_cat[idx[sep_constraint]])

        mag5 = 25 - 2.5 * np.log10(sex_matches['FLUX_APER_5'] / exp_t)
        mag6 = 25 - 2.5 * np.log10(sex_matches['FLUX_APER_6'] / exp_t)

        sex_matches['MAG_ISO'] = mag5 * 0.186 + mag6 * 0.814 - 0.05 * ss_seeing[i] + 0.04   # temporarily save mag
        sex_matches['PETRO_RADIUS'] = mag5 * 0.0 + ss_airmass[i]                             # temporarily save seeing
        sex_matches['MAGERR_ISO'] = - 2.5 * np.log10(
            (sex_matches['FLUX_APER_6'] - sex_matches['FLUXERR_APER_6']) / exp_t) \
                                    + 2.5 * np.log10(
            (sex_matches['FLUX_APER_6'] + sex_matches['FLUXERR_APER_6']) / exp_t)
        if i == 0:
            sex_stack = sex_matches
            cat_stack = cat_matches
        else:
            sex_stack = vstack([sex_stack, sex_matches])
            cat_stack = vstack([cat_stack, cat_matches])

    plt.close('all')
    plt.figure(figsize=(5,5))
    plt.title('SS '+ss_date+' '+band+'-band standardized DECam vs trimmed cat')

    minchi = 1e9
    a_best = 0
    b_best = 0
    for j in a_guess:
        for k in b_guess:
            chi = 0
            for l in range(0, len(sex_stack)):
                chi += np.abs(cat_stack[mag_col][l] - sex_stack['MAG_ISO'][l] - j - k * sex_stack['PETRO_RADIUS'][l]) \
                       / np.sqrt(cat_stack[err_col][l]**2 + sex_stack['MAGERR_ISO'][l]**2)
            if chi < minchi:
                minchi = chi
                a_best = j
                b_best = k

    for i in range(0, len(ss_name)):
        sex_result = ascii.read(sex_dir+ss_name[i]+'_'+band+'a.cat')
        sex_result_good = (sex_result['MAG_ISO'] < 90) & (sex_result['FLAGS'] <= flag_lim) & \
                          (sex_result['CLASS_STAR'] > class_star_lim)
        sex_coords = SkyCoord(sex_result['ALPHA_J2000'][sex_result_good],
                              sex_result['DELTA_J2000'][sex_result_good],
                              unit='deg')

        ss_cat = ascii.read(cat_dir+ss_cat_name[i]+'.dat.trimmed', format='no_header')

        cat_coords = SkyCoord(ss_cat['col2'], ss_cat['col3'], unit=(u.hourangle, u.deg))

        idx, d2d, d3d = sex_coords.match_to_catalog_sky(cat_coords)
        sep_constraint = (d2d.arcsec < max_sep) & (ss_cat[idx][err_col] > 0)
        sex_matches = sex_result[sex_result_good][sep_constraint]
        cat_matches = ss_cat[idx[sep_constraint]]

        mag5 = 25 - 2.5 * np.log10(sex_matches['FLUX_APER_5'] / exp_t)
        mag6 = 25 - 2.5 * np.log10(sex_matches['FLUX_APER_6'] / exp_t)
        mag14_8 = mag5 * 0.186 + mag6 * 0.814 - 0.05 * ss_seeing[i] + 0.04  # Aperture correction
        mag6err = - 2.5 * np.log10((sex_matches['FLUX_APER_6'] - sex_matches['FLUXERR_APER_7']) / exp_t) \
        + 2.5 * np.log10((sex_matches['FLUX_APER_6'] + sex_matches['FLUXERR_APER_6']) / exp_t)

        # col4: u_mag, col7: g_mag, col10: r_mag
        plt.scatter(cat_matches[mag_col], mag14_8 +a_best +b_best*ss_airmass[i], alpha=0.5, s=2,
                    label=ss_legend_name[i]+' ({:d})'.format(int(sample_numbers[i])))
        # plt.errorbar(cat_matches['col4'], sex_matches['MAG_AUTO']-ss_u_zp[i]+25.0, xerr=cat_matches['col5'],
        #               yerr=sex_matches['MAGERR_AUTO'], fmt='.', capthick=2, alpha=0.1)

    plt.plot([10, 22], [10, 22], '--', alpha=0.1)
    plt.xlim([10, 22])
    plt.ylim([10, 22])
    plt.ylabel('DECam MAG_APER standardized')
    plt.xlabel('Southern Star catalog')
    plt.gca().invert_xaxis()
    plt.gca().invert_yaxis()
    plt.text(21, 11, band+' = '+band+'_inst (ZP=25) + {:5.3f} + {:5.3f} * X'.format(a_best, b_best))
    plt.text(21, 12, r'$\chi^2$/# = {:7.3f}'.format(minchi/np.sum(sample_numbers)))
    plt.text(21, 13, '# of sample = {:4d}'.format(int(np.sum(sample_numbers))))
    plt.text(20.5, 14, 'match < {:4.2f}\"'.format(max_sep))
    plt.text(20.5, 15, 'FLAGS <= {:2d}'.format(int(flag_lim)))
    plt.text(20.5, 16, 'CLASS_STAR > {:4.2f}'.format(class_star_lim))
    plt.legend(title='Standard Star Field (Airmass)(#)', loc='lower right', prop={'size':7})
    plt.savefig(plot_dir+ss_date+'_DECam_'+band+'_band_standardized_vs_trimmed_aper_corr6.pdf')
    plt.show(block=False)

    print(a_guess)
    print(b_guess)

    print("--- %s minutes ---" % ((time.time() - start_time)/60.0))

    Funcs.print_time()
