from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
import time
import my_module
from astropy.table import Table, vstack

plot_dir=("/Users/duhokim/work/abell/plot/")
cat_dir=("/Users/duhokim/old_macbook_documents/work/cat/SouthernStandardStars/remove_asterisk/")
check_dir=("/Volumes/APPLE SSD/Server/work/sex/ss_deblend/")

dates = ['19aug', '20aug', '21aug']

band = 'r'
exp_t = 10
mag_col = 'col10'
err_col = 'col11'

for m in range(0, len(dates)):
    start_time = time.time()
    if m == 0:
        ss_date = '19aug'
        sex_dir=("/Users/duhokim/work/abell/sex/cat/ss_deblend/19aug/")
        ss_name = ['220100', 'E8-A', 'LSE_259', '190000']
        ss_legend_name = ['220100 (1.96)', 'E8-A (1.21)', 'LSE_259 (1.75)', '190000 (1.23)']
        ss_cat_name = ['220100-300000', 'E8_a', 'LSE_259', '190000-295600']
        ss_airmass = [1.96, 1.21, 1.75, 1.23]
        ss_seeing = [1.04, 0.95, 1.03, 0.92]

        max_sep = 1.0
        class_star_lim = 0.97
        flag_lim = 0
        near_nei_dis_lim = 1.0  # in aperture radii which is 7.43"
        excl_num = [[], [], [], [3078, 8282]]

        a_guess = 0.5308 + np.linspace(-0.0001, 0.0001, 11)
        b_guess = -0.1018 + np.linspace(-0.0001, 0.0001, 11)

    elif m == 1:
        ss_date = '20aug'
        sex_dir=("/Users/duhokim/work/abell/sex/cat/ss_deblend/20aug/")
        ss_name = ['E5-A', 'E6-A', 'LSE_259_1_11', '220100_2_44', '220100_1_01', 'LSE_259_1_71']
        ss_legend_name = ['E5-A (1.7)', 'E6-A (1.12)', 'LSE_259 (1.11)', '220100 (2.44)', '220100 (1.01)', 'LSE_259 (1.71']
        ss_cat_name = ['E5_a', 'E6_a', 'LSE_259', '220100-300000', '220100-300000', 'LSE_259']
        ss_airmass = [1.7, 1.12, 1.11, 2.44, 1.01, 1.71]
        ss_seeing = [1.75, 1.49, 1.31, 1.77, 3.5, 3.8]

        max_sep = 1.0
        class_star_lim = 0.8
        flag_lim = 0
        near_nei_dis_lim = 1.0
        excl_num = [[], [1101], [], [], [], []]

        a_guess = 0.5272 + np.linspace(-0.0001, 0.0001, 11)
        b_guess = -0.0847 + np.linspace(-0.0001, 0.0001, 11)

    elif m == 2:
        ss_date = '21aug'
        sex_dir=("/Users/duhokim/work/abell/sex/cat/ss_deblend/21aug/")
        ss_name = ['E5-A', 'E6-A', 'LSE_259_1_11', 'LSE_259_1_8']
        ss_legend_name = ['E5-A (1.65)', 'E6-A (1.11)', 'LSE_259 (1.11)', 'LSE_259 (1.8)']
        ss_cat_name = ['E5_a', 'E6_a', 'LSE_259', 'LSE_259']
        ss_airmass = [1.65, 1.11, 1.11, 1.8]
        ss_seeing = [1.46, 1.03, 1.07, 1.58]

        max_sep = 1.0
        class_star_lim = 0.85
        flag_lim = 0
        near_nei_dis_lim = 1.0
        excl_num = [[], [], [], [], [], []]

        a_guess = 0.4817 + np.linspace(-0.0001, 0.0001, 11)
        b_guess = -0.0457 + np.linspace(-0.0001, 0.0001, 11)

    sample_numbers = np.zeros(len(ss_name))
    for i in range(0, len(ss_name)):
        sex_result = ascii.read(sex_dir + ss_name[i] + '_' + band + 'a.cat')
        good_one = []
        for j in sex_result['NUMBER']:
            good_one.append(j not in excl_num[i])
        sex_result_good = (sex_result['MAG_ISO'] < 90) & (sex_result['FLAGS'] <= flag_lim) & \
                          (sex_result['CLASS_STAR'] > class_star_lim) & good_one
        sex_coords = SkyCoord(sex_result['ALPHA_J2000'][sex_result_good],
                              sex_result['DELTA_J2000'][sex_result_good],
                              unit='deg')

        ss_cat = ascii.read(cat_dir + ss_cat_name[i] + '.dat.trimmed', format='no_header')
        cat_coords = SkyCoord(ss_cat['col2'], ss_cat['col3'], unit=(u.hourangle, u.deg))

        idx, d2d, d3d = sex_coords.match_to_catalog_sky(cat_coords)
        sep_constraint = (d2d.arcsec < max_sep) & (ss_cat[idx][err_col] > 0) & (ss_cat[idx]['col20'] > near_nei_dis_lim)
        sex_matches = sex_result[sex_result_good][sep_constraint]
        cat_matches = ss_cat[idx[sep_constraint]]

        sample_numbers[i] = len(ss_cat[idx[sep_constraint]])

        mag5 = 25 - 2.5 * np.log10(sex_matches['FLUX_APER_5'] / exp_t)
        mag6 = 25 - 2.5 * np.log10(sex_matches['FLUX_APER_6'] / exp_t)

        sex_matches['MAG_ISO'] = mag5 * 0.186 + mag6 * 0.814 - 0.05 * ss_seeing[i] + 0.04   # temporarily save mag after aperture correction
        sex_matches['PETRO_RADIUS'] = ss_airmass[i]                             # temporarily save airmass
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
    plt.title('SS '+ss_date+' '+band+' standard star fit SEx deblending conf.')

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
        good_one = []
        for j in sex_result['NUMBER']:
            good_one.append(j not in excl_num[i])
        sex_result_good = (sex_result['MAG_ISO'] < 90) & (sex_result['FLAGS'] <= flag_lim) & \
                          (sex_result['CLASS_STAR'] > class_star_lim) & good_one
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
        mag6err = - 2.5 * np.log10((sex_matches['FLUX_APER_6'] - sex_matches['FLUXERR_APER_6']) / exp_t) \
        + 2.5 * np.log10((sex_matches['FLUX_APER_6'] + sex_matches['FLUXERR_APER_6']) / exp_t)

        # col4: u_mag, col7: g_mag, col10: r_mag
        plt.scatter(cat_matches[mag_col], mag14_8 +a_best +b_best*ss_airmass[i], alpha=0.5, s=2,
                    label=ss_legend_name[i]+' ({:d})'.format(int(sample_numbers[i])))
        # plt.errorbar(cat_matches['col4'], sex_matches['MAG_AUTO']-ss_u_zp[i]+25.0, xerr=cat_matches['col5'],
        #               yerr=sex_matches['MAGERR_AUTO'], fmt='.', capthick=2, alpha=0.1)

        with open(check_dir + ss_date + '/' + ss_name[i] + '_' + band + '.reg', 'w') as reg:
            reg.writelines("# Region file format: DS9 version 4.1\n global color=green dashlist=8 3 width=4 \
            font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 \
            source=1\n fk5\n")
            for j in range(0, len(sex_matches)):
                if abs(cat_matches[mag_col][j] - (mag14_8[j] +a_best +b_best*ss_airmass[i])) > 1.0:
                    color = 'red'
                    print('{} {}'.format(ss_name[i], sex_matches['NUMBER'][j]))
                    reg.writelines("circle({}, {}, {}\") # text=\'{}\',color={}\n".format(
                        cat_matches['col2'][j],
                        cat_matches['col3'][j],
                        4.0,
                        'SS',
                        color))
                    reg.writelines("circle({:12.7f}, {:12.7f}, {:7.3f}\") # text=\'{}\', color={}\n".format(
                        sex_matches['ALPHA_J2000'][j],
                        sex_matches['DELTA_J2000'][j],
                        4.0,
                        'SEx {:4.1f}'.format(d2d[sep_constraint][j].arcsec),
                        color))
                else:
                    color = 'green'
                    reg.writelines("circle({}, {}, {}\") # text=\'{}\',color={}\n".format(
                        cat_matches['col2'][j],
                        cat_matches['col3'][j],
                        4.0,
                        '',
                        color))

    plt.plot([10, 22], [10, 22], '--', alpha=0.1)
    plt.xlim([10, 22])
    plt.ylim([10, 22])
    plt.ylabel('DECam MAG_APER at 14.8\" (aper corrected)')
    plt.xlabel('Southern Star catalog')
    plt.gca().invert_xaxis()
    plt.gca().invert_yaxis()
    plt.text(21.5, 11,
             band + ' = ' + band + r"$_{{inst}}$ (ZP=25) + " + '{:7.4f} + {:7.4f} * X (airmass)'.format(a_best, b_best))
    plt.text(21.5, 12, r'$\chi^2$/# = {:7.3f}'.format(minchi/np.sum(sample_numbers)))
    plt.text(21.5, 13, '# of sample = {:4d}'.format(int(np.sum(sample_numbers))))
    plt.text(20.8, 13.5, 'match < {:4.2f}\"'.format(max_sep))
    plt.text(20.8, 14, 'FLAGS <= {:2d}'.format(int(flag_lim)))
    plt.text(20.8, 14.5, 'CLASS_STAR > {:4.2f}'.format(class_star_lim))
    plt.text(20.8, 15, 'nei_dist > {:3.1f} x Aper'.format(near_nei_dis_lim))
    plt.text(20.8, 15.5, 'Aper : 7.43\" for SS')
    plt.legend(title='Standard Star Field (Airmass)(#)', loc='lower right', prop={'size':7})
    plt.savefig(plot_dir+ss_date+'_DECam_'+band+'_band_standardized_vs_trimmed_aper_corr_flag0.png')
    plt.show(block=False)

    print(a_guess)
    print(b_guess)

    print("--- %s minutes ---" % ((time.time() - start_time)/60.0))

    my_module.print_time()