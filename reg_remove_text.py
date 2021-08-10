import abell_cluster_module as ab
from astropy.io import ascii

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


# sift through bright sources to exclude garbage SEx detections
for cluster in ab.clusters:
    fn = ab.sex_dir+f'DECam_merged_SEx_cat_{cluster}_Gal_ext_corrected_{ab.ver}'
    with open(fn+'_kron.reg') as reg, open(fn+'_20rmag.reg', 'w') as new_rg:
        lines = reg.readlines()
        cat = ascii.read(fn+'.txt')
        for i in range(0, len(cat)):
            if cat['MAG_AUTO_r'][i] < 20:
                line = lines[i]
                new_rg.writelines(line)
