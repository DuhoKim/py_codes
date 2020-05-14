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

band = 'r'

mag_col_u = 'col4'
err_col_u = 'col5'
a_best_u = 3.349
b_best_u = -0.318

mag_col_g = 'col7'
err_col_g = 'col8'
a_best_g = 2.857
b_best_g = -0.183

mag_col_r = 'col10'
err_col_r = 'col11'
a_best_r = 3.024
b_best_r = -0.117

max_sep = 1.0

plt.figure(figsize=(5,5))
plt.title('SS 19 Aug '+band+'-band standardized DECam r vs g-r')


for i in range(0, len(ss_name)):
    sex_result = ascii.read(sex_dir+ss_name[i]+'_'+band+'wi.cat')
    sex_result_good = (sex_result['MAG_ISO'] < 90) & (sex_result['FLAGS'] == 0) & (sex_result['CLASS_STAR'] > 0.98)
    sex_coords = SkyCoord(sex_result['ALPHA_J2000'][sex_result_good], sex_result['DELTA_J2000'][sex_result_good],
                          unit='deg')

    ss_cat = ascii.read(cat_dir+ss_cat_name[i]+'.dat.trimmed', format='no_header')

    cat_coords = SkyCoord(ss_cat['col2'], ss_cat['col3'], unit=(u.hourangle, u.deg))

    idx, d2d, d3d = sex_coords.match_to_catalog_sky(cat_coords)
    sep_constraint = (d2d.arcsec < max_sep) & (ss_cat[idx][err_col_r] > 0)
    sex_matches = sex_result[sex_result_good][sep_constraint]
    cat_matches = ss_cat[idx[sep_constraint]]

    # col4: u_mag, col7: g_mag, col10: r_mag
    plt.scatter(cat_matches[mag_col_g] - cat_matches[mag_col_r], sex_matches['MAG_ISO'] +a_best_r +b_best_r*ss_airmass[i],  alpha=0.5, s=2,
                label=ss_legend_name[i])
    # plt.errorbar(cat_matches['col4'], sex_matches['MAG_AUTO']-ss_u_zp[i]+25.0, xerr=cat_matches['col5'],
    #               yerr=sex_matches['MAGERR_AUTO'], fmt='.', capthick=2, alpha=0.1)

#plt.plot([10, 25], [0, 25], '--', alpha=0.1)
plt.xlim([0, 2])
plt.ylim([10, 25])
plt.ylabel('DECam MAG_ISO standardized (r)')
plt.xlabel('Southern Star catalog (g-r)')
plt.gca().invert_xaxis()
plt.gca().invert_yaxis()
plt.text(1.9, 11, band+' = '+band+'_inst (ZP=25) + {:5.3f} + {:5.3f} * X'.format(a_best_r, b_best_r))
#plt.text(-2, 12, r'$\chi^2$/# = {:7.3f}'.format(minchi/np.sum(sample_numbers)))
#plt.text(-2, 13, '# of sample = {:4d}'.format(int(np.sum(sample_numbers))))
plt.text(1.9, 12, 'match < 1\"')
plt.text(1.9, 13, 'FLAGS = 0')
plt.text(1.9, 14, 'CLASS_STAR > 0.98')
plt.legend(title='Standard Star Field (Airmass)(#)', loc='lower right', prop={'size':7})
plt.savefig(plot_dir+'19Aug_DECam_'+band+'_band_standardized_color.pdf')
plt.show()
