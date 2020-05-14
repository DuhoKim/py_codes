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
import pandas as pd
from os import listdir

sex_dir=("/Users/dkim108/Documents/work/sex/gals/")
plot_dir=("/Users/dkim108/Documents/work/plot/")
cat_dir=("/Users/dkim108/Documents/work/cat/")

cluster = 'A2399'
date = '19aug'

band = 'r'
date == '19aug'

max_sep = 1.0

am_cat = pd.read_csv('/Users/dkim108/Documents/work/cat/airmass.csv')

sdss_galaxy_cat_all = ascii.read(cat_dir + 'sdss_galaxy_' + cluster + '.csv')
sdss_galaxy_cat = sdss_galaxy_cat_all[:200]
cat_coords = coords.SkyCoord(sdss_galaxy_cat['ra'], sdss_galaxy_cat['dec'], unit=(u.deg, u.deg))

f = listdir(sex_dir + cluster)
exps = [s for s in f if '_'+band+'i_' in s and date in s]

a_guess = 1 + np.linspace(-0.5, 0.5, 5)
b_guess = -0.15 + np.linspace(-0.1, 0.1, 5)
chi = np.zeros((len(a_guess), len(b_guess)))

for i in range(0, len(exps)):
    time = exps[i].split('_')[3]
    time_with_colon = time[:2] + ':' + time[2:]
    am_match = am_cat.loc[(am_cat['date'] == 19) & (am_cat['time'] == time_with_colon)]
    X = am_match.iloc[0, 4]
    t = am_match.iloc[0, 3]

    exp_indi = ascii.read(sex_dir + cluster + '/' + exps[i])
    exp_indi_good = exp_indi['MAG_ISO'] < 90
    exp_coords = SkyCoord(exp_indi['ALPHA_J2000'][exp_indi_good], exp_indi['DELTA_J2000'][exp_indi_good], unit='deg')

    idx, d2d, d3d = exp_coords.match_to_catalog_sky(cat_coords)
    sep_constraint = d2d.arcsec < max_sep
    exp_matches = exp_indi[exp_indi_good][sep_constraint]
    sdss_matches = sdss_galaxy_cat[idx[sep_constraint]]

    for j in range(0, len(a_guess)):
        for k in range(0, len(b_guess)):
            for l in range(0, len(exp_matches)):
                chi[j, k] += (exp_matches['MAG_ISO'][l] + 2.5 * np.log10(t) + a_guess[j] + b_guess[k] * X - sdss_matches['petroMag_r'][l])**2 \
                               / (exp_matches['MAGERR_ISO'][l]**2 + sdss_matches['petroMagErr_r'][l]**2)

min_a_ind, min_b_ind = np.where(chi == np.min(chi))

plt.figure(figsize=(5,5))
plt.title(cluster+' '+date+' '+band+'-band DECam vs SDSS Galaxy cat')

for i in range(0, len(exps)):
    time = exps[i].split('_')[3]
    time_with_colon = time[:2] + ':' + time[2:]
    am_match = am_cat.loc[(am_cat['date'] == 19) & (am_cat['time'] == time_with_colon)]
    X = am_match.iloc[0, 4]

    exp_indi = ascii.read(sex_dir + cluster + '/' + exps[i])
    exp_indi_good = exp_indi['MAG_ISO'] < 90
    exp_coords = SkyCoord(exp_indi['ALPHA_J2000'][exp_indi_good], exp_indi['DELTA_J2000'][exp_indi_good], unit='deg')

    idx, d2d, d3d = exp_coords.match_to_catalog_sky(cat_coords)
    sep_constraint = d2d.arcsec < max_sep
    exp_matches = exp_indi[exp_indi_good][sep_constraint]
    sdss_matches = sdss_galaxy_cat[idx[sep_constraint]]

    plt.scatter(sdss_matches['petroMag_r'], exp_matches['MAG_ISO'] + a_guess[min_a_ind][0] + b_guess[min_b_ind][0] * X,
        alpha=0.2, s=1, label='{:4.2f} {:3d}'.format(X, am_match.iloc[0, 3]))

plt.plot([5, 25], [5, 25], '--', alpha=0.1)
plt.xlim([5, 25])
plt.ylim([5, 25])
plt.ylabel('DECam MAG_ISO standardized')
plt.xlabel('SDSS Galaxy catalog petroMag_r')
plt.gca().invert_xaxis()
plt.gca().invert_yaxis()
plt.text(24, 6, band+' = '+band+'_inst (ZP=25) + {:5.3f} + {:5.3f} * X'.format(a_guess[min_a_ind][0], b_guess[min_b_ind][0]))
plt.text(24, 7, r'$\chi^2$/# = {:7.3f}'.format(np.min(chi)/len(exps)))
plt.text(24, 8, '# of exposures = {:4d}'.format(len(exps)))
plt.text(24, 9, 'match < {:4.2f}\"'.format(max_sep))
#plt.text(18.5, 13, 'FLAGS < {:2d}'.format(int(flag_lim)))
#plt.text(18.5, 14, 'CLASS_STAR > {:4.2f}'.format(class_star_lim))
plt.legend(title='Airmass', loc='lower right', prop={'size':7})
plt.savefig(plot_dir+date+'_DECam_'+band+'_band_standardized_vs_SDSS_Galaxy_catalog.pdf')
plt.show(block=False)

print(a_guess)
print(b_guess)