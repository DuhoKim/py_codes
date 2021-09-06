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

## from SDSS
# fn = "/Users/duhokim/work/abell/spec/SDSS/A2670"
# with open(fn+'.reg', 'w') as reg:
#     cat = ascii.read(fn+'.csv')
#     for j in range(0, len(cat)):
#         if cat['z'][j] > 0.05 and cat['z'][j] < 0.1:
#             color = 'green'
#         else:
#             color = 'red'
#         reg.writelines(f"j2000; circle({cat['ra'][j]}, {cat['dec'][j]}, 10\") # color={color}\n")

# fn = "/Users/duhokim/work/abell/spec/"
# with open(fn+'A2670_spec.radec', 'w') as mem:
#     cnt = 0
#     cat = ascii.read(fn + 'SDSS/A2670.csv')
#     for i in range(0, len(cat)):
#         if cat['z'][i] > 0.066 and cat['z'][i] < 0.086:
#             cnt = cnt + 1
#             # mem.writelines(f"{cnt} {cat['ra'][i]} {cat['dec'][i]} {cat['z'][i]} {cat['zErr'][i]} \n")
#             mem.writelines(f"{cat['ra'][i]} {cat['dec'][i]} \n")

## from Visir catalog
# fn = "/Users/duhokim/work/abell/spec/"
# with open(fn+'A3558_spec_Shapley.txt', 'w') as spec:
#     cat = ascii.read(fn+'.dat')
#     for j in range(0, len(cat)):
#         coord_sha = SkyCoord(f"{cat['col3'][j]}:{cat['col4'][j]}:{cat['col5'][j]} "
#                              f"{cat['col6'][j]}:{cat['col7'][j]}:{cat['col8'][j]}",
#                             unit = (u.hourangle, u.deg))
#         reg.writelines(f"{coord_sha.ra.value} {coord_sha.dec.value} \n")
# fn = "/Users/duhokim/work/abell/spec/WINGS/"
# with open(fn + 'OmegaWINGS_Moretti+2017_A754.reg', 'w') as reg:
#     cat = ascii.read(fn + 'OmegaWINGS_Moretti+2017_A754.txt')
#     for j in range(0, len(cat)):
#         if cat['col10'][j]:
#             reg.writelines(f"j2000; circle({cat['col4'][j]}, {cat['col5'][j]} 10\") # width=3, color=black \n")

### KYDISC + Shapley
# fn = "/Users/duhokim/work/abell/spec/"
# with open(fn+'A3574_spec.radec', 'w') as merged:
#     cnt = 0
#     cat = ascii.read(fn+'KYDISC/Oh+2019_A3574.txt')
#     for j in range(0, len(cat)):
#         cnt = cnt + 1
#         # coord_sha = SkyCoord(f"{cat['col2'][j]}:{cat['col3'][j]}:{cat['col4'][j]} "
#         #                      f"{cat['col5'][j]}:{cat['col6'][j]}:{cat['col7'][j]}",
#         #                     unit = (u.hourangle, u.deg))
#         # merged.writelines(f"{cnt} {cat['col10'][j]} {cat['col11'][j]} {cat['col12'][j]} {cat['col13'][j]} \n")
#         merged.writelines(f"{cat['col10'][j]} {cat['col11'][j]} \n")
#
#     sha = ascii.read(fn+'Shapley/catalog.dat')
#     for i in range(0, len(sha)):
#         z = sha['col18'][i] / const.c * 1e3
#         ze = sha['col19'][i] / const.c * 1e3
#         if z > 0.01 and z < 0.02:
#             cnt = cnt + 1
#             coord_sha = SkyCoord(f"{sha['col3'][i]}:{sha['col4'][i]}:{sha['col5'][i]} "
#                                  f"{sha['col6'][i]}:{sha['col7'][i]}:{sha['col8'][i]}",
#                                 unit = (u.hourangle, u.deg))
#             # merged.writelines(f"{cnt} {coord_sha.ra.value} {coord_sha.dec.value} {z} {ze} \n")
#             merged.writelines(f"{coord_sha.ra.value} {coord_sha.dec.value} \n")


### KYDISC
# fn = "/Users/duhokim/work/abell/spec/"
# with open(fn+'A3659_spec.radec', 'w') as merged:
#     cnt = 0
#     cat = ascii.read(fn+'KYDISC/Oh+2019_A3659.txt')
#     for j in range(0, len(cat)):
#         cnt = cnt + 1
#         # coord_sha = SkyCoord(f"{cat['col2'][j]}:{cat['col3'][j]}:{cat['col4'][j]} "
#         #                      f"{cat['col5'][j]}:{cat['col6'][j]}:{cat['col7'][j]}",
#         #                     unit = (u.hourangle, u.deg))
#         # merged.writelines(f"{cnt} {cat['col10'][j]} {cat['col11'][j]} {cat['col12'][j]} {cat['col13'][j]} \n")
#         merged.writelines(f"{cat['col10'][j]} {cat['col11'][j]} \n")

## merge WINGS+OmegaWINGS
# fn = "/Users/duhokim/work/abell/spec/"
# with open(fn + 'A754_spec.radec', 'w') as merged:
#     cnt = 0
#     wings = ascii.read(fn + 'WINGS/WINGS_Cava+2009_A754.txt')
#     for i in range(0, len(wings)):
#         if wings['col8'][i]:
#             cnt = cnt + 1
#             z = wings['col5'][i]/const.c*1e3
#             coord = SkyCoord(wings['col2'][i], wings['col3'][i], unit='deg')
#             dist = mm.clcendist(ab.redshifts[0], z, ab.coords_cl[0], coord)
#             # merged.writelines(f"{cnt} {coord.ra.value} {coord.dec.value} {z} {wings['col6'][i]/const.c*1e3} {dist} \n")
#             merged.writelines(
#                 f"{coord.ra.value} {coord.dec.value} \n")
#     omega = ascii.read(fn + 'WINGS/OmegaWINGS_Moretti+2017_A754.txt')
#     for i in range(0, len(omega)):
#         if omega['col10'][i]:
#             cnt = cnt + 1
#             z = omega['col6'][i]
#             coord = SkyCoord(omega['col4'][i], omega['col5'][i], unit='deg')
#             dist = mm.clcendist(ab.redshifts[0], z, ab.coords_cl[0], coord)
#             # merged.writelines(f"{cnt} {coord.ra.value} {coord.dec.value} {z} {omega['col7'][i]} {dist} \n")
#             merged.writelines(f"{coord.ra.value} {coord.dec.value} \n")

# fn = "/Users/duhokim/work/abell/spec/"
# with open(fn + 'A2399_spec.radec', 'w') as merged:
#     cnt = 0
#     wings = ascii.read(fn + 'WINGS/WINGS_Cava+2009_A2399.txt')
#     for i in range(0, len(wings)):
#         if wings['col8'][i]:
#             cnt = cnt + 1
#             z = wings['col5'][i]/const.c*1e3
#             coord = SkyCoord(wings['col2'][i], wings['col3'][i], unit='deg')
#             dist = mm.clcendist(ab.redshifts[0], z, ab.coords_cl[0], coord)
#             # merged.writelines(f"{cnt} {coord.ra.value} {coord.dec.value} {z} {wings['col6'][i]/const.c*1e3} {dist} \n")
#             merged.writelines(
#                 f"{coord.ra.value} {coord.dec.value} \n")
#     omega = ascii.read(fn + 'WINGS/OmegaWINGS_Moretti+2017_A2399.txt')
#     for i in range(0, len(omega)):
#         if omega['col10'][i]:
#             cnt = cnt + 1
#             z = omega['col6'][i]
#             coord = SkyCoord(omega['col4'][i], omega['col5'][i], unit='deg')
#             dist = mm.clcendist(ab.redshifts[0], z, ab.coords_cl[0], coord)
#             # merged.writelines(f"{cnt} {coord.ra.value} {coord.dec.value} {z} {omega['col7'][i]} {dist} \n")
#             merged.writelines(f"{coord.ra.value} {coord.dec.value} \n")

## merge OmegaWINGS+Shapley
# fn = "/Users/duhokim/work/abell/spec/"
# with open(fn + 'A3558_spec.radec', 'w') as merged:
#     cnt = 0
#     omega = ascii.read(fn + 'WINGS/OmegaWINGS_Moretti+2017_A3558.txt')
#     for i in range(0, len(omega)):
#         if omega['col10'][i]:
#             cnt = cnt + 1
#             # merged.writelines(f"{cnt} {omega['col4'][i]} {omega['col5'][i]} {omega['col6'][i]} {omega['col7'][i]} \n")
#             merged.writelines(f"{omega['col4'][i]} {omega['col5'][i]} \n")
#
#     sha = ascii.read(fn+'Shapley/catalog.dat')
#
#     for i in range(0, len(sha)):
#         z = sha['col18'][i] / const.c * 1e3
#         ze = sha['col19'][i] / const.c * 1e3
#         if z > 0.032 and z < 0.06:
#             cnt = cnt + 1
#             coord_sha = SkyCoord(f"{sha['col3'][i]}:{sha['col4'][i]}:{sha['col5'][i]} "
#                                  f"{sha['col6'][i]}:{sha['col7'][i]}:{sha['col8'][i]}",
#                                 unit = (u.hourangle, u.deg))
#             # merged.writelines(f"{cnt} {coord_sha.ra.value} {coord_sha.dec.value} {z} {ze} \n")
#             merged.writelines(f"{coord_sha.ra.value} {coord_sha.dec.value} \n")

###  OmegaWINGS
fn = "/Users/duhokim/work/abell/spec/"
with open(fn + 'A3716_spec.radec', 'w') as merged:
    cnt = 0
    omega = ascii.read(fn + 'WINGS/OmegaWINGS_Moretti+2017_A3716.txt')
    for i in range(0, len(omega)):
        if omega['col10'][i]:
            cnt = cnt + 1
            # merged.writelines(f"{cnt} {omega['col4'][i]} {omega['col5'][i]} {omega['col6'][i]} {omega['col7'][i]} \n")
            merged.writelines(f"{omega['col4'][i]} {omega['col5'][i]} \n")

### from Topcat catalog
# fn = "/Users/duhokim/work/abell/sex/cat/A3558_g-r"
# with open(fn+'.reg', 'w') as reg:
#     cat = ascii.read(fn+'_blue')
#     color = 'blue'
#     for j in range(0, len(cat)):
#         reg.writelines(f"j2000; circle({cat['col2'][j]}, {cat['col3'][j]}, "
#                        f"10\") # width=3, color={color} \n")
#     cat = ascii.read(fn+'_red')
#     color = 'red'
#     for j in range(0, len(cat)):
#         reg.writelines(f"j2000; circle({cat['col2'][j]}, {cat['col3'][j]}, "
#                        f"10\") # width=3, color={color} \n")

### from SEx catalog
# for k in range(0, len(ab.clusters)):
#     with open(ab.short_cat_dir + ab.short_cat_fn[k][2] + '_deblend.reg', 'w') as reg:
#         cat = ascii.read(ab.short_cat_dir + ab.short_cat_fn[k][2] + '_deblend.cat')
#         for j in range(0, len(cat)):
#             reg.writelines(f"j2000; ellipse({cat[j]['ALPHA_J2000']}, {cat[j]['DELTA_J2000']}, "
#                            f"{cat[j]['A_WORLD'] * 3600 * cat[j]['KRON_RADIUS']}\","
#                            f"{cat[j]['B_WORLD'] * 3600 * cat[j]['KRON_RADIUS']}\","
#                            f"{180 - cat[j]['THETA_WORLD']}) # text=\'{cat[j]['NUMBER']}\' \n")
# kron = 3.5
# for k in range(0, len(ab.clusters)):
#     with open(ab.work_dir + 'sex/reg/' + ab.clusters[k] + '_rmag_16.reg', 'w') as reg16,    \
#         open(ab.work_dir + 'sex/reg/' + ab.clusters[k] + '_rmag_18.reg', 'w') as reg18, \
#         open(ab.work_dir + 'sex/reg/' + ab.clusters[k] + '_rmag_25.reg', 'w') as reg25:
#         cat = ascii.read(ab.sex_dir+f'DECam_merged_SEx_cat_{ab.clusters[k]}_Gal_ext_corrected_{ab.ver}_dqm_edit.txt')
#         for j in range(0, len(cat)):
#             if cat[j]['MAG_AUTO_r'] < 16:
#                 reg16.writelines(f"j2000; circle({cat[j]['ALPHA_J2000']}, {cat[j]['DELTA_J2000']}, "
#                                f"{cat[j]['A_WORLD'] * kron}\") \n")
#                 reg18.writelines(f"j2000; circle({cat[j]['ALPHA_J2000']}, {cat[j]['DELTA_J2000']}, "
#                                  f"{cat[j]['A_WORLD'] * kron}\") \n")
#                 reg25.writelines(f"j2000; circle({cat[j]['ALPHA_J2000']}, {cat[j]['DELTA_J2000']}, "
#                                  f"{cat[j]['A_WORLD'] * kron}\") \n")
#             elif cat[j]['MAG_AUTO_r'] < 18:
#                 reg18.writelines(f"j2000; circle({cat[j]['ALPHA_J2000']}, {cat[j]['DELTA_J2000']}, "
#                                f"{cat[j]['A_WORLD'] * kron}\") \n")
#                 reg25.writelines(f"j2000; circle({cat[j]['ALPHA_J2000']}, {cat[j]['DELTA_J2000']}, "
#                                  f"{cat[j]['A_WORLD'] * kron}\") \n")
#             else:
#                 reg25.writelines(f"j2000; circle({cat[j]['ALPHA_J2000']}, {cat[j]['DELTA_J2000']}, "
#                                  f"{cat[j]['A_WORLD'] * kron}\") \n")



# with open(short_cat_dir + 'A3716_rs_300_0418_deblend.reg', 'w') as new_rg:
#     new_rg.writelines('# Region file format: DS9 version 4.1\n')
#     new_rg.writelines(f'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 \
#     highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
#     new_rg.writelines('fk5\n')
#
#     cat = ascii.read(short_cat_dir + 'A3716_rs_300_0418_deblend.cat')
#     for j in range(0, len(cat)):
#         new_rg.writelines("j2000; ellipse({:12.7f}, {:12.7f}, {:7.3f}\", {:7.3f}\", {:7.3f}) # "
#                                "color={} \n".format(
#             cat[j]['ALPHA_J2000'],
#             cat[j]['DELTA_J2000'],
#             cat[j]['A_WORLD'] * 3600 * cat[j]['KRON_RADIUS'],
#             cat[j]['B_WORLD'] * 3600 * cat[j]['KRON_RADIUS'],
#             180 - cat[j]['THETA_WORLD'],
#             'green'))



# with open(out_dir + 'A754_satu.reg', 'w') as new_rg:
#     new_rg.writelines('# Region file format: DS9 version 4.1\n')
#     new_rg.writelines(f'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 \
#     highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
#     new_rg.writelines('fk5\n')
#
#     cat = ascii.read("/Users/duhokim/work/abell/sex/reg/A754_satu.txt")
#     for j in range(0, len(cat)):
#         new_rg.writelines("circle({:12.7f}, {:12.7f}, {:7.3f}\") # text=\'{}\' \n".format(
#             cat[j]['col2'],
#             cat[j]['col3'],
#             cat[j]['col4'] *  5,
#             cat[j]['col1']))

# for i in range(0, len(fns)):
#     with open(out_dir + 'A754_' + fns[i] + '.reg', 'w') as new_rg:
#         new_rg.writelines('# Region file format: DS9 version 4.1\n')
#         new_rg.writelines(f'global color={color[i]} dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 \
#         highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
#         new_rg.writelines('fk5\n')
#
#         cat = ascii.read(short_cat_dir+'A754_rs_300_'+fns[i]+'.cat')
#         for j in range(0, len(cat)):
#             new_rg.writelines("ellipse({:12.7f}, {:12.7f}, {:7.3f}\", {:7.3f}\", {:7.3f}) # \n".format(
#                 cat[j]['ALPHA_J2000'],
#                 cat[j]['DELTA_J2000'],
#                 cat[j]['A_WORLD'] * 3600 * cat[j]['KRON_RADIUS'],
#                 cat[j]['B_WORLD'] * 3600 * cat[j]['KRON_RADIUS'],
#                 180 - cat[j]['THETA_WORLD']))