from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
from astroquery.sdss import SDSS
from astropy.coordinates import SkyCoord
from astropy.coordinates import FK4, Angle, Latitude, Longitude
from astropy import coordinates as coords
import os
import astropy.units as u
from astropy.table import Table, vstack
# import matplotlib
# matplotlib.rcParams['text.usetex'] = True

sex_dir=("/Users/dkim108/Documents/work/sex/19aug/")
plot_dir=("/Users/dkim108/Documents/work/plot/")
cat_dir=("/Users/dkim108/Documents/work/cat/SouthernStandardStars/")

ss_name = ['220100', 'E8-A', 'LSE_259', '190000']
ss_legend_name = ['220100 (1.96)', 'E8-A (1.21)', 'LSE_259 (1.75)', '190000 (1.23)']
ss_cat_name = ['220100-300000', 'E8_a', 'LSE_259', '190000-295600']
ss_airmass = [1.96, 1.21, 1.75, 1.23]

band = 'u'

mag_col_u = 'col4'
err_col_u = 'col5'
mag_col_g = 'col7'

max_sep = 1.0

plt.figure(figsize=(5,5))
plt.title('SS 19 Aug '+band+'-band standardized DECam vs trimmed cat')

a_guess = 3.349 + np.linspace(-0.005, 0.005, 7)
b_guess = -0.318 + np.linspace(-0.005, 0.005, 7)
c_guess = 0.010 + np.linspace(-0.005, 0.005, 7)

minchi = 1e9
c_best = 0
sample_numbers = np.zeros(len(ss_name))
for j in a_guess:
    for k in b_guess:
        for m in c_guess:
            chi = 0
            for i in range(0, len(ss_name)):
                sex_result = ascii.read(sex_dir+ss_name[i]+'_'+band+'wi.cat')
                sex_result_good = (sex_result['MAG_ISO'] < 90) & (sex_result['FLAGS'] == 0) & (sex_result['CLASS_STAR'] > 0.98)
                sex_coords = SkyCoord(sex_result['ALPHA_J2000'][sex_result_good], sex_result['DELTA_J2000'][sex_result_good],
                                      unit='deg')

                ss_cat = ascii.read(cat_dir+ss_cat_name[i]+'.dat.trimmed', format='no_header')

                cat_coords = SkyCoord(ss_cat['col2'], ss_cat['col3'], unit=(u.hourangle, u.deg))

                # col4: u_mag, col7: g_mag, col10: r_mag
                idx, d2d, d3d = sex_coords.match_to_catalog_sky(cat_coords)
                sep_constraint = (d2d.arcsec < max_sep) & (ss_cat[idx][err_col_u] > 0)
                sex_matches = sex_result[sex_result_good][sep_constraint]
                cat_matches = ss_cat[idx[sep_constraint]]

                for l in range(0, len(sex_matches)):
                    chi += np.abs(cat_matches[mag_col_u][l] - sex_matches['MAG_ISO'][l] - j - k * ss_airmass[i] \
                            - m * ( cat_matches[mag_col_u][l] - cat_matches[mag_col_g][l] ) ) \
                            / np.sqrt(cat_matches[err_col_u][l]**2 + sex_matches['MAGERR_ISO'][l]**2)
                    if (j == a_guess[0]) and (k == b_guess[0]) and (m == c_guess[0]):
                        sample_numbers[i] += 1
            if chi < minchi:
                minchi = chi
                a_best = j
                b_best = k
                c_best = m

for i in range(0, len(ss_name)):
    sex_result = ascii.read(sex_dir+ss_name[i]+'_'+band+'wi.cat')
    sex_result_good = (sex_result['MAG_ISO'] < 90) & (sex_result['FLAGS'] == 0) & (sex_result['CLASS_STAR'] > 0.98)
    sex_coords = SkyCoord(sex_result['ALPHA_J2000'][sex_result_good], sex_result['DELTA_J2000'][sex_result_good],
                          unit='deg')

    ss_cat = ascii.read(cat_dir+ss_cat_name[i]+'.dat.trimmed', format='no_header')

    cat_coords = SkyCoord(ss_cat['col2'], ss_cat['col3'], unit=(u.hourangle, u.deg))

    idx, d2d, d3d = sex_coords.match_to_catalog_sky(cat_coords)
    sep_constraint = (d2d.arcsec < max_sep) & (ss_cat[idx][err_col_u] > 0)
    sex_matches = sex_result[sex_result_good][sep_constraint]
    cat_matches = ss_cat[idx[sep_constraint]]

    # col4: u_mag, col7: g_mag, col10: r_mag
    plt.scatter(cat_matches[mag_col_u], sex_matches['MAG_ISO'] + a_best + b_best * ss_airmass[i] + c_best * \
                (cat_matches[mag_col_u] - cat_matches[mag_col_g]), alpha=0.5, s=2,
                label=ss_legend_name[i]+'({:3d})'.format(int(sample_numbers[i])))
    # plt.errorbar(cat_matches['col4'], sex_matches['MAG_AUTO']-ss_u_zp[i]+25.0, xerr=cat_matches['col5'],
    #               yerr=sex_matches['MAGERR_AUTO'], fmt='.', capthick=2, alpha=0.1)

plt.plot([10, 25], [10, 25], '--', alpha=0.1)
plt.xlim([10, 25])
plt.ylim([10, 25])
plt.ylabel('DECam MAG_ISO standardized')
plt.xlabel('Southern Star catalog')
plt.gca().invert_xaxis()
plt.gca().invert_yaxis()
plt.text(24, 11, band+' = '+band+'_inst (ZP=25) + {:5.3f} + {:5.3f} * X + {:5.3f} * (u-g)'.format(a_best, b_best, c_best))
plt.text(24, 12, r'$\chi^2$/# = {:7.3f}'.format(minchi/np.sum(sample_numbers)))
plt.text(24, 13, '# of sample = {:4d}'.format(int(np.sum(sample_numbers))))
plt.text(23.5, 14, 'match < 1\"')
plt.text(23.5, 15, 'FLAGS = 0')
plt.text(23.5, 16, 'CLASS_STAR > 0.98')
plt.legend(title='Standard Star Field (Airmass)(#)', loc='lower right', prop={'size':7})
plt.savefig(plot_dir+'19Aug_DECam_'+band+'_band_standardized_vs_trimmed_color.pdf')
plt.show()

print(a_guess)
print(b_guess)
print(c_guess)
