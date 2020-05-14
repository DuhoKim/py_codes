from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
from astroquery.sdss import SDSS
from astropy.coordinates import SkyCoord
from astropy import coordinates as coords
import os
import astropy.units as u
from astropy.table import Table, vstack
import Funcs

data_dir=("/Users/dkim108/Documents/work/sex/")
plot_dir=("/Users/dkim108/Documents/work/plot/")
cat_dir=("/Users/dkim108/Documents/work/cat/")

cluster_name = "A2670"
band = "r"
pos = coords.SkyCoord(358.557, -10.419, unit="deg")
max_sep = 1.0

sex_result_decam = ascii.read(data_dir+cluster_name+"_"+band+"_nthresh64_back64.cat")
sex_result_sdss = ascii.read(cat_dir+"sdss_"+cluster_name+"_"+band+".cat")
#sex_result_good = (sex_result['FLAGS'] < 4) & (sex_result['MAG_AUTO'] < 90)
sex_result_decam_good = sex_result_decam['MAG_AUTO'] < 90
sex_result_sdss_good = sex_result_sdss['MAG_AUTO'] < 90
sex_coords_decam = SkyCoord(sex_result_decam['ALPHA_J2000'][sex_result_decam_good], sex_result_decam['DELTA_J2000'][sex_result_decam_good], unit='deg')
sex_coords_sdss = SkyCoord(sex_result_sdss['ALPHA_J2000'][sex_result_sdss_good], sex_result_sdss['DELTA_J2000'][sex_result_sdss_good], unit='deg')

#sdss_galaxy_cat = ascii.read(cat_dir+'sdss_galaxy_'+cluster_name+'.csv')
#cat_coords = coords.SkyCoord(sdss_galaxy_cat['ra'], sdss_galaxy_cat['dec'], unit=(u.deg, u.deg))

# idx, d2d, d3d = sex_coords_sdss.match_to_catalog_sky(cat_coords)
# sep_constraint = d2d.arcsec < max_sep
# sex_matches_sdss = sex_result_sdss[sex_result_sdss_good][sep_constraint]
# sex_matches_cat = sdss_galaxy_cat[idx[sep_constraint]]

idx, d2d, d3d = sex_coords_sdss.match_to_catalog_sky(sex_coords_decam)
sep_constraint = d2d.arcsec < max_sep
sex_matches_sdss = sex_result_sdss[sex_result_sdss_good][sep_constraint]
sex_matches_decam = sex_result_decam[sex_result_decam_good][idx[sep_constraint]]

plt.figure(figsize=(5,5))
plt.title(cluster_name+' '+band+'-band MAG_AUTO SDSS vs DECam')
plt.scatter(sex_matches_sdss['MAG_AUTO'], sex_matches_decam['MAG_AUTO'], alpha=0.5, s=2)
plt.plot([0, 30], [0,30], '--', alpha=0.1)
plt.xlim([5, 30])
plt.ylim([5, 30])
plt.xlabel('SDSS MAG_AUTO')
plt.ylabel('DECam MAG_AUTO')
plt.gca().invert_xaxis()
plt.gca().invert_yaxis()
plt.savefig(plot_dir+cluster_name+'_'+band+'_SDSS_MAG_AUTO_vs_DECam_MAG_AUTO.pdf')
plt.show()
