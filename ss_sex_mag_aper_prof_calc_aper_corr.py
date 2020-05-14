from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
import Funcs

plot_dir=("/Users/dkim108/Documents/work/plot/")
cat_dir=("/Users/dkim108/Documents/work/cat/SouthernStandardStars/")

bands = ['u', 'g', 'r']
ss_date = ['19aug', '20aug', '21aug']
rad_in_arcsec = np.multiply( [4, 8, 12, 20, 30, 40, 60, 100], 0.263)

num = 30

fig, ax = plt.subplots()

for m in range(0, len(bands)):
    band = bands[m]
    if band == 'u':  # col4: u_mag, col7: g_mag, col10: r_mag
        mag_col = 'col4'
        err_col = 'col5'
    elif band == 'g':
        mag_col = 'col7'
        err_col = 'col8'
    elif band == 'r':
        mag_col = 'col10'
        err_col = 'col11'

    aper_corr_mean = np.zeros(14)   # mean of aperture correction at the aperture 14.8" for 14 SS fields
    aper_corr_std = np.zeros(14)    # std
    seeing = np.zeros(14)           # seeing
    k = 0
    for j in range(0, len(ss_date)):
        if ss_date[j] == '19aug':
            sex_dir=("/Users/dkim108/Documents/work/sex/19aug/")
            ss_name = ['220100', 'E8-A', 'LSE_259', '190000']
            ss_legend_name = ['220100 (1.96)(1.0)', 'E8-A (1.21)(1.0)', 'LSE_259 (1.75)(1.0)', '190000 (1.23)(0.9)']
            ss_cat_name = ['220100-300000', 'E8_a', 'LSE_259', '190000-295600']
            ss_airmass = [1.96, 1.21, 1.75, 1.23]
            ss_seeing = [1.04, 0.95, 1.03, 0.92]
            ss_col = ['indianred', 'tan', 'olivedrab', 'mediumpurple']

        elif ss_date[j] == '20aug':
            sex_dir=("/Users/dkim108/Documents/work/sex/20aug/")
            ss_name = ['E5-A', 'E6-A', 'LSE_259_1_11', '220100_2_44', '220100_1_01', 'LSE_259_1_71']
            ss_legend_name = ['E5-A (1.7)(1.8)', 'E6-A (1.12)(1.5)', 'LSE_259 (1.11)(1.3)', '220100 (2.44)(1.8)', '220100 (1.01)(3.5)', 'LSE_259 (1.71)(4.0)']
            ss_cat_name = ['E5_a', 'E6_a', 'LSE_259', '220100-300000', '220100-300000', 'LSE_259']
            ss_airmass = [1.7, 1.12, 1.11, 2.44, 1.01, 1.71]
            ss_seeing = [1.75, 1.49, 1.31, 1.77, 3.5, 3.8]
            ss_col = ['darkred', 'darkslategray', 'navy', 'purple', 'darkolivegreen', 'darkorange']

        elif ss_date[j] == '21aug':
            sex_dir=("/Users/dkim108/Documents/work/sex/21aug/")
            ss_name = ['E5-A', 'E6-A', 'LSE_259_1_11', 'LSE_259_1_8']
            ss_legend_name = ['E5-A  1.65 1.5 ', 'E6-A  1.11 1.1 ', 'LSE_259  1.11 1.05 ', 'LSE_259  1.8 1.6 ']
            ss_cat_name = ['E5_a', 'E6_a', 'LSE_259', 'LSE_259']
            ss_airmass = [1.65, 1.11, 1.11, 1.8]
            ss_seeing = [1.46, 1.03, 1.07, 1.58]
            ss_col = ['darkred', 'darkslategray', 'navy', 'purple']


        for i in range(0, len(ss_name)):
            sex_result = ascii.read(sex_dir+ss_name[i]+'_'+band+'a.cat')
            sex_result_good = (sex_result['MAG_ISO'] < 90) & (sex_result['FLAGS'] == 0) & \
                              (sex_result['CLASS_STAR'] > 0.95)
            star_cat_all = sex_result[sex_result_good]
            star_cat_all.sort('MAG_ISO')
            star_cat = star_cat_all[0:num]


            aper_corr_temp = np.zeros(len(star_cat))

            for l in range(0, len(star_cat)):
                flux_annul = [star_cat['FLUX_APER'][l],
                        star_cat['FLUX_APER_1'][l] - star_cat['FLUX_APER'][l],
                        star_cat['FLUX_APER_2'][l] - star_cat['FLUX_APER_1'][l],
                        star_cat['FLUX_APER_3'][l] - star_cat['FLUX_APER_2'][l],
                        star_cat['FLUX_APER_4'][l] - star_cat['FLUX_APER_3'][l],
                        star_cat['FLUX_APER_5'][l] - star_cat['FLUX_APER_4'][l],
                        star_cat['FLUX_APER_6'][l] - star_cat['FLUX_APER_5'][l],
                        star_cat['FLUX_APER_7'][l] - star_cat['FLUX_APER_6'][l]]

                flux_tot = [star_cat['FLUX_APER'][l],
                              star_cat['FLUX_APER_1'][l],
                              star_cat['FLUX_APER_2'][l],
                              star_cat['FLUX_APER_3'][l],
                              star_cat['FLUX_APER_4'][l],
                              star_cat['FLUX_APER_5'][l],
                              star_cat['FLUX_APER_6'][l],
                              star_cat['FLUX_APER_7'][l]]

                mags = 2.5 * np.log10(flux_tot/max(flux_tot))
                aper_corr_temp[l] = mags[5] * 0.186 + mags[6] * 0.814   # interpolated mag at 14.8"
            aper_corr_mean[k] = np.nanmean(aper_corr_temp)
            aper_corr_std[k] = np.nanstd(aper_corr_temp)
            seeing[k] = ss_seeing[i]
            k += 1

    ax.errorbar(seeing, aper_corr_mean, yerr=aper_corr_std, fmt='o', alpha = 0.5, label=bands[m])

    a_guess = -0.05 + np.linspace(-0.005, 0.005, 5)
    b_guess = 0.04 + np.linspace(-0.005, 0.005, 5)

    minchi = 1e6
    for n in range(0, len(a_guess)):
        for o in range(0, len(b_guess)):
            chi = 0
            for l in range(0, len(seeing)):
                chi += np.abs(a_guess[n] * seeing[l] - aper_corr_mean[l] + b_guess[o]) \
                        / np.sqrt(a_guess[n]**2 + 1) \
                        / aper_corr_std[l]
            if chi < minchi:
                minchi = chi
                a_best = a_guess[n]
                b_best = b_guess[o]

    ax.plot([0, 4], [b_best, a_best * 4 + b_best], alpha=0.5,
            label=bands[m]+' fit:{:4.2f}*seeing+{:4.2f}'.format(a_best, b_best))

ax.plot([0, 4], [0, 0], linestyle=':', alpha=0.5)
ax.set_xlim([0, 4])
ax.set_ylim([-0.1, 0.05])
ax.set_xlabel('Seeing \"')
ax.set_ylabel('Aperture correction at 14.8\" [mag]')
ax.legend(title='Brightest sample # = '+str(num)+' for each')
plt.savefig(plot_dir + 'DECam_Aper_Corr_at_14_8_brightest_'+str(num)+'.pdf')
plt.close()