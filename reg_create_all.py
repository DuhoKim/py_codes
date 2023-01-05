from astropy.io import ascii
import abell_cluster_module as ab
from astropy.coordinates import SkyCoord
from scipy import constants as const
import my_module as mm
import numpy as np
from astropy import units as u
import importlib
importlib.reload(mm)
importlib.reload(ab)

work_dir="/Users/duhokim/work/abell/spec/"

clusters = ['A754', 'A2399', 'A3558', 'A3716']

for cluster in clusters:
    with open(work_dir + f'{cluster}_OWM_ra_dec_z.txt', 'w') as merged:
        cnt = 0
        omega = ascii.read(work_dir + f'WINGS/OmegaWINGS_Moretti+2017_{cluster}.txt')
        coords_omega = SkyCoord(omega['col4'], omega['col5'], unit='deg')

        # only members
        omem = (omega['col10'] == 1)

        for i in range(0, len(omega[omem])):
            cnt = cnt + 1
            z = omega['col6'][omem][i]
            merged.writelines(f"{cnt} {coords_omega[omem][i].ra.value} "
                              f"{coords_omega[omem][i].dec.value} {z}  \n")
        print(f'{cluster} n_o={len(omem)}, n_o_mem={sum(omem)}, tot={cnt}')