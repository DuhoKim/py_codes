from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
#from astroquery.sdss import SDSS
from astropy.coordinates import SkyCoord
from astropy import wcs
from astropy.io import fits
from astropy.coordinates import FK4, Angle, Latitude, Longitude
from astropy import coordinates as coords
import os
import astropy.units as u
from astropy.table import Table, vstack
# import matplotlib
# matplotlib.rcParams['text.usetex'] = True

sex_dir=("/Users/dkim108/Documents/work/sex/19aug/")
plot_dir=("/Users/dkim108/Documents/work/plot/")
cat_dir=("/Users/dkim108/Documents/work/DECam/ss/")

ss_name = ['220100', 'E8-A', 'LSE_259', '190000']
#ss_name = ['220100']
ss_legend_name = ['220100 (1.96)', 'E8-A (1.21)', 'LSE_259 (1.75)', '190000 (1.23)']
ss_cat_name = ['220100-300000', 'E8_a', 'LSE_259', '190000-295600']
ss_airmass = [1.96, 1.21, 1.75, 1.23]
ss_col = ['darkorange', 'yellowgreen', 'aqua', 'orchid']

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

fig = plt.figure(figsize=(6, 6))

plt.title('SS 19 Aug '+band+'-band MAG_ISO vs IRAF phot')

h1 = fig.add_axes([0.1, 0.2, 0.7, 0.7])
h2 = fig.add_axes([0.1, 0.1, 0.7, 0.1], sharex=h1)
h3 = fig.add_axes([0.8, 0.2, 0.1, 0.7], sharey=h1)

for i in range(0, len(ss_name)):
    sex_result = ascii.read(sex_dir+ss_name[i]+'_'+band+'wi.cat')
    sex_result_good = (sex_result['MAG_ISO'] < 90) & (sex_result['FLAGS'] == 0) & (sex_result['CLASS_STAR'] > 0.98)
    sex_coords = SkyCoord(sex_result['ALPHA_J2000'][sex_result_good], sex_result['DELTA_J2000'][sex_result_good],
                          unit='deg')

    hdulist = fits.open(cat_dir+'19aug/'+ss_name[i]+'_'+band+'.fits')
    for j in range(1, 61):
        w = wcs.WCS(hdulist[j].header)

        ss_cat = ascii.read(cat_dir+'19aug_phot/'+ss_name[i]+'_'+band+'.mag.'+str(j))
        pixcrd = [[ss_cat['XINIT'][k], ss_cat['YINIT'][k]] for k in range(len(ss_cat))]
        wldcrd = w.wcs_pix2world(pixcrd, 0) # have 0-based (Numpy-like) coordinates
        cat_coords = SkyCoord([k[0] for k in wldcrd], [k[1] for k in wldcrd], unit='deg')

        idx, d2d, d3d = sex_coords.match_to_catalog_sky(cat_coords)
        sep_constraint = (d2d.arcsec < max_sep)
        sex_matches = sex_result[sex_result_good][sep_constraint]
        cat_matches = ss_cat[idx[sep_constraint]]

        # col4: u_mag, col7: g_mag, col10: r_mag
        if j == 1:
            h1.scatter(cat_matches['MAG']+2.5, sex_matches['MAG_ISO']+2.5,  alpha=0.5, s=2, color=ss_col[i],
                                                                        label=ss_legend_name[i])
        else:
            h1.scatter(cat_matches['MAG'] + 2.5, sex_matches['MAG_ISO'] + 2.5, alpha=0.5, s=2, color=ss_col[i])
        h2.scatter(cat_matches['MAG'] + 2.5, sex_matches['MAG_ISO'] - cat_matches['MAG'], alpha=0.5, s=2,
                    color=ss_col[i])
        h3.scatter(sex_matches['MAG_ISO'] - cat_matches['MAG'], sex_matches['MAG_ISO'] + 2.5, alpha=0.5, s=2,
                   color=ss_col[i])

        # plt.errorbar(cat_matches['col4'], sex_matches['MAG_AUTO']-ss_u_zp[i]+25.0, xerr=cat_matches['col5'],
        #               yerr=sex_matches['MAGERR_AUTO'], fmt='.', capthick=2, alpha=0.1)

h1.plot([10, 20], [10, 20], '--', alpha=0.1)
h2.plot([10, 20], [0, 0], '--', alpha=0.1)
h3.plot([0, 0], [10, 20], '--', alpha=0.1)
h1.set_xlim([20, 10])
h1.set_ylim([20, 10])
h1.set_ylabel('MAG_ISO '+band+'-band mag')
h2.set_ylabel('MAG_ISO - IRAF PHOT')
h2.set_xlabel('IRAF phot task'+band+'-band mag')
h3.set_yticklabels([])
# h1.gca().invert_xaxis()
# h1.gca().invert_yaxis()
#h1.legend(title='Standard Star Field (Airmass)(#)', loc='lower right', prop={'size':7})
plt.savefig(plot_dir+'19Aug_DECam_'+band+'_band_MAG_ISO_vs_aper_corr_phot.pdf')
plt.show()

