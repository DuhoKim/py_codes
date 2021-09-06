import numpy as np
from astropy import units as u
from astropy.io import ascii
from astropy.coordinates import SkyCoord

sheen = ascii.read("/Users/duhokim/work/abell/cat/sheen_2012_table_6.txt")
gala = ascii.read('/Users/duhokim/work/abell/galapagos/A2670_test/A2670_back_values.txt')

coords_sheen = SkyCoord(sheen['col2'], sheen['col3'], unit=(u.hourangle, u.deg))
coords_gala = SkyCoord(gala['ra'], gala['dec'], unit='deg')

idx_s2g, d2d_g2s, d3d = coords_gala.match_to_catalog_sky(coords_sheen)
match_g2s = (d2d_g2s.arcsec < 1)

for i in range(0, sum(match_g2s)):
    if i in idx_s2g[match_g2s]:
        print(f"{sheen['col1'][idx_s2g[match_g2s][i]]} {gala[match_g2s][i]['id']}")