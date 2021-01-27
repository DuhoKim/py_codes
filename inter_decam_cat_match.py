from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
from astroquery.sdss import SDSS
from astropy.coordinates import SkyCoord
from astropy import coordinates as coords
import os
import astropy.units as u
from astropy.table import Table, vstack
import my_module

data_dir=("/Users/dkim108/Documents/work/sex/")
plot_dir=("/Users/dkim108/Documents/work/plot/")
cat_dir=("/Users/dkim108/Documents/work/cat/")

cluster_name = "A2670"
band = "g"
pos = coords.SkyCoord(358.557, -10.419, unit="deg")
max_sep = 1.0

sex_result_i = ascii.read(data_dir+"a2670_g_nthresh64_back64.cat")
sex_result_j = ascii.read(cat_dir+"decam_"+cluster_name+"_"+band+"_osj.cat")
#sex_result_good = (sex_result['FLAGS'] < 4) & (sex_result['MAG_AUTO'] < 90)
sex_result_i_good = sex_result_i['MAG_AUTO'] < 90
sex_result_j_good = sex_result_j['MAG_AUTO'] < 90
sex_coords_i = SkyCoord(sex_result_i['ALPHA_J2000'][sex_result_i_good], sex_result_i['DELTA_J2000'][sex_result_i_good], unit='deg')
sex_coords_j = SkyCoord(sex_result_j['ALPHA_J2000'][sex_result_j_good], sex_result_j['DELTA_J2000'][sex_result_j_good], unit='deg')

#ss_cat_tot = ascii.read(cat_dir+'sdss_star_'+cluster_name+'.csv')

#cat_coords = coords.SkyCoord(ss_cat_tot['ra'], ss_cat_tot['dec'], unit=(u.deg, u.deg))
idx, d2d, d3d = sex_coords_i.match_to_catalog_sky(sex_coords_j)
sep_constraint = d2d.arcsec < max_sep
sex_matches_i = sex_result_i[sex_result_i_good][sep_constraint]
sex_matches_j = sex_result_j[sex_result_j_good][idx[sep_constraint]]


plt.figure(figsize=(5,5))
plt.title(cluster_name+' '+band+'-band mag osi vs osj')
plt.scatter(sex_matches_i['MAG_AUTO'], sex_matches_j['MAG_AUTO'], alpha=0.5, s=2)
plt.plot([0, 30], [0,30], '--', alpha=0.1)
plt.xlim([5, 30])
plt.ylim([5, 30])
plt.xlabel('MAG_AUTO of osi')
plt.ylabel('MAG_AUTO of osj')
plt.gca().invert_xaxis()
plt.gca().invert_yaxis()
plt.savefig(plot_dir+cluster_name+'_'+band+'_mag_SEx_osi_vs_osj.pdf')
plt.show()
