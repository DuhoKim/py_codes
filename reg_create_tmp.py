from astropy.io import ascii
import abell_cluster_module as ab
from astropy.coordinates import SkyCoord
from scipy import constants as const
import my_module as mm
from astropy import units as u
import importlib
importlib.reload(mm)
importlib.reload(ab)

work_dir=("/Users/duhokim/work/abell/cat/")
short_cat_dir=("/Volumes/APPLE SSD/Server/work/sex/best_single/")
out_dir=("/Users/duhokim/work/abell/sex/reg/")

fns = ['n32_m05', 'n32_m0005', 'n64_m005', 'n16_m005', 'default', 'deblend']
color = ['blue', 'cyan', 'magenta', 'yellow', 'green', 'red']
# white
# black
# red
# green
# blue
# cyan
# magenta
# yellow

###  A754  ####
fn = "/Users/duhokim/work/abell/catalogue/matched_flag0_5sig_dqm/"
with open(fn + 'A2670_r_lower.reg', 'w') as reg:
    cat = ascii.read(fn + 'A2670_r_lower.cat')

    coords_cat = SkyCoord(cat['ALPHA_J2000_1'], cat['DELTA_J2000_1'], unit='deg')

    for i in range(0, len(cat)):
        reg.writelines(f"j2000; circle({cat['ALPHA_J2000_1'][i]}, {cat['DELTA_J2000_1'][i]} 10\") # width=3, color=red \n")

with open(fn + 'A2670_r_upper.reg', 'w') as reg:
    cat = ascii.read(fn + 'A2670_r_upper.cat')

    coords_cat = SkyCoord(cat['ALPHA_J2000_1'], cat['DELTA_J2000_1'], unit='deg')

    for i in range(0, len(cat)):
        reg.writelines(
            f"j2000; circle({cat['ALPHA_J2000_1'][i]}, {cat['DELTA_J2000_1'][i]} 10\") # width=3, color=blue \n")
