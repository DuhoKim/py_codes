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

sex_dir=("/Users/dkim108/Documents/work/sex/19aug/")
plot_dir=("/Users/dkim108/Documents/work/plot/")
cat_dir=("/Users/dkim108/Documents/work/cat/SouthernStandardStars/")

ss_name = ['220100', 'E8-A', 'LSE_259', '190000']
ss_legend_name = ['220100 (1.96)', 'E8-A (1.21)', 'LSE_259 (1.75)', '190000 (1.23)']
ss_cat_name = ['220100-300000', 'E8_a', 'LSE_259', '190000-295600']
ss_airmass = [1.96, 1.21, 1.75, 1.23]
ss_u_zp = [24.5029, 25.252, 24.2819, 24.9523]
ss_g_zp = [26.4458, 26.7149, 26.3407, 26.5773]
ss_r_zp = [28.6607, 28.4248, 28.1103, 28.8181]

u_sam_num = 256
g_sam_num = 273
r_sam_num = 272

band = "r"
max_sep = 1.0

plt.figure(figsize=(5,5))
plt.title('SS 19 Aug '+band+'-band standardized DECam vs trimmed cat')

a_guess = 3.045 + np.linspace(-0.05, 0.05, 10)
b_guess = -0.165 + np.linspace(-0.005, 0.005, 10)

minchi = 1e9
a_best = 0
b_best = 0
for j in a_guess:
    for k in b_guess:
        chi=0
        #sample_number=0
        for i in range(0, len(ss_name)):
            sex_result = ascii.read(sex_dir+ss_name[i]+'_'+band+'.cat')
            #sex_result_good = sex_result['MAG_AUTO'] < 90
            sex_coords = SkyCoord(sex_result['ALPHA_J2000'], sex_result['DELTA_J2000'], unit='deg')

            ss_cat = ascii.read(cat_dir+ss_cat_name[i]+'.dat.trimmed', format='no_header')

            # cat_coords = coords.SkyCoord(ss_cat['col2'], ss_cat['col3'], unit=(u.deg, u.deg))
            cat_coords = SkyCoord(ss_cat['col2'], ss_cat['col3'], unit=(u.hourangle, u.deg))

            # col4: u_mag, col7: g_mag, col10: r_mag
            idx, d2d, d3d = sex_coords.match_to_catalog_sky(cat_coords)
            sep_constraint = (d2d.arcsec < max_sep) & (ss_cat[idx]['col11'] > 0)
            sex_matches = sex_result[sep_constraint]
            cat_matches = ss_cat[idx[sep_constraint]]

            for l in range(0, len(sex_matches)):
                chi += np.abs(cat_matches['col10'][l] - (sex_matches['MAG_AUTO'][l]-ss_r_zp[i]+25.0) - j
                              - k * ss_airmass[i]) / np.sqrt(cat_matches['col11'][l]**2 + sex_matches['MAGERR_AUTO'][l]**2)
                #sample_number += 1
        if chi < minchi:
            minchi = chi
            a_best = j
            b_best = k

for i in range(0, len(ss_name)):
    sex_result = ascii.read(sex_dir+ss_name[i]+'_'+band+'.cat')
    #sex_result_good = sex_result['MAG_AUTO'] < 90
    sex_coords = SkyCoord(sex_result['ALPHA_J2000'], sex_result['DELTA_J2000'], unit='deg')

    ss_cat = ascii.read(cat_dir+ss_cat_name[i]+'.dat.trimmed', format='no_header')

    # cat_coords = coords.SkyCoord(ss_cat['col2'], ss_cat['col3'], unit=(u.deg, u.deg))
    cat_coords = SkyCoord(ss_cat['col2'], ss_cat['col3'], unit=(u.hourangle, u.deg))

    idx, d2d, d3d = sex_coords.match_to_catalog_sky(cat_coords)
    sep_constraint = (d2d.arcsec < max_sep) & (ss_cat[idx]['col11'] > 0)
    sex_matches = sex_result[sep_constraint]
    cat_matches = ss_cat[idx[sep_constraint]]

    # col4: u_mag, col7: g_mag, col10: r_mag
    plt.scatter(cat_matches['col10'], sex_matches['MAG_AUTO']-ss_r_zp[i]+25.0 +a_best +b_best*ss_airmass[i], alpha=0.5, s=2, label=ss_legend_name[i])
    # plt.errorbar(cat_matches['col4'], sex_matches['MAG_AUTO']-ss_u_zp[i]+25.0, xerr=cat_matches['col5'], yerr=sex_matches['MAGERR_AUTO'],
    #              fmt='.', capthick=2, alpha=0.1)

plt.plot([10, 25], [10, 25], '--', alpha=0.1)
plt.xlim([10, 25])
plt.ylim([10, 25])
plt.ylabel('DECam MAG_AUTO standardized')
plt.xlabel('Southern Star catalog')
plt.gca().invert_xaxis()
plt.gca().invert_yaxis()
plt.text(24, 11,'r = r_inst + '+str(a_best)+' + '+str(b_best)+' * AM')
plt.text(24, 12,'chi^2 = '+str(minchi/r_sam_num))
#plt.legend(title='Standard Star Field (Airmass)')
plt.savefig(plot_dir+'19Aug_DECam_'+band+'_band_standardized_vs_trimmed_cat.pdf')
plt.show()
