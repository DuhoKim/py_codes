import abell_cluster_module as ab
from astropy.io import ascii
from astropy.coordinates import SkyCoord
import numpy as np

# for cluster in ab.clusters:
#     fn = ab.sex_dir+f'DECam_merged_SEx_cat_{cluster}_Gal_ext_corrected_{ab.ver}_kron'
#     with open(fn+'.reg') as reg, open(fn+'_notext.reg', 'w') as new_rg:
#         lines = reg.readlines()
#         for line in lines:
#             if 'text' in line:
#                 split = line.split("\'")
#                 new_line = (split[0]+split[2]).replace('text=', '')
#                 new_rg.writelines(new_line)
#             else:
#                 new_rg.writelines(line)

# # sift through bright sources to exclude garbage SEx detections
# for cluster in ab.clusters:
#     fn = ab.sex_dir+f'DECam_merged_SEx_cat_{cluster}_Gal_ext_corrected_{ab.ver}'
#     with open(fn+'_kron.reg') as reg, open(fn+'_20rmag.reg', 'w') as new_rg:
#         lines = reg.readlines()
#         cat = ascii.read(fn+'.txt')
#         for i in range(0, len(cat)):
#             if cat['MAG_AUTO_r'][i] < 20:
#                 line = lines[i]
#                 new_rg.writelines(line)

# add 'circle' in front
fit_dir = f'{ab.work_dir}fits/extracted_mock/'
sex_dir = f'{ab.work_dir}sex/run_stack_mock/'
# for cluster in ab.clusters:
# with open(f'{sex_dir}{cluster}_rsi_1_mock_in.reg', 'w') as reg, \
#     open(f'{sex_dir}{cluster}_rsi_1_mock_in.cat', 'r') as cat:
with open(f'{fit_dir}A754_rsi_1_mock_in_2028.reg', 'w') as reg:
    #open(f'{fit_dir}A754_rsi_1_mock_out_2028.reg', 'w') as reg_out, \
    cat_in = ascii.read(f'{fit_dir}A754_rsi_1_mock_in_2028.cat')
    cat_out = ascii.read(f'{sex_dir}A754_rsi_1_mock_out_2028.cat')
    # open(f'{fit_dir}A754_rsi_1_mock_in_2028.cat', 'r') as cat_in, \
    # open(f'{sex_dir}A754_rsi_1_mock_out_2028.cat', 'r') as cat_out:
    # reg.writelines('# Region file format: DS9 version 4.1\n')
    # reg.writelines(
    #     'global color=green dashlist=8 3 width=2 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    # reg.writelines('fk5\n')

    matched_in = np.ones(len(cat_in), dtype=bool)

    for i in range(0, len(cat_in)):
        for j in range(0, len(cat_out)):
            matched_in[i] = (abs(cat_in['x'][i] - cat_out['X_IMAGE'][j]) < 2) & \
                            (abs(cat_in['y'][i] - cat_out['Y_IMAGE'][j]) < 2)
            if matched_in[i]:
                break

    # in_coords = SkyCoord(cat_in['x'], cat_in['y'], unit='deg')
    # out_coords = SkyCoord(cat_out['ALPHA_J2000'], cat_out['DELTA_J2000'], unit='deg')
    #
    # idx, d2d, d3d = in_coords.match_to_catalog_sky(out_coords)
    # sep_constraint = d2d.arcsec < 1.0

    for i in range(len(cat_in)):
        if matched_in[i]:
            col = 'green'
        else:
            col = 'red'
        reg.writelines(f'circle {cat_in["x"][i]} {cat_in["y"][i]} 15 # color={col} dash=1 width=4 \n')
    #
    # sex_matches = sex_cat[sep_constraint]
    #
    # lines = cat_in.readlines()
    # new_lines = ['global color=red dashlist=8 3 width=4 \n']
    # for line in lines:
    #     if not line[0] == '#':
    #         new_lines.append(f'circle {line}')
    # reg_in.writelines(new_lines)
    #
    # lines = cat_out.readlines()
    # new_lines = ['global color=green dashlist=8 3 width=4 \n']
    # for line in lines:
    #     if not line[0] == '#':
    #         new_lines.append(f'circle {line}')
    # reg_out.writelines(new_lines)