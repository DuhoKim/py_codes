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

work_dir=("/Users/duhokim/work/abell/cat/")
short_cat_dir=("/Volumes/APPLE SSD/Server/work/sex/best_single/")
out_dir=("/Users/duhokim/work/abell/sex/reg/")

z_ran = 0.05
r_ran = 1.5 # deg

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

## KYDISC + Shapley
fn = "/Users/duhokim/work/abell/spec/"
with open(fn+f'A3574_spec_ra_dec_z_zran_{z_ran}_rran_{r_ran}.txt', 'w') as merged:
    k = 4
    cnt = 0
    #   Shapley Supercluster
    sha = ascii.read(fn+'Shapley/catalog.dat')
    coords_sha_str = []
    for i in range(len(sha)):
        coords_sha_str.append(f"{sha['col3'][i]}:{sha['col4'][i]}:{sha['col5'][i]} "
                              f"{sha['col6'][i]}:{sha['col7'][i]}:{sha['col8'][i]}")
    coords_sha = SkyCoord(coords_sha_str, unit=(u.hourangle, u.deg))

    sha_z = sha['col18'] / const.c * 1e3
    sha_ze = sha['col19'] / const.c * 1e3

    smem = (coords_sha.separation(ab.coords_cl_cen[k]).value < r_ran) & \
           (sha_z > ab.redshifts[4] - z_ran) & (sha_z < ab.redshifts[4] + z_ran)
    sha = sha[smem]

    for i in range(0, len(sha)):
        cnt = cnt + 1
        coord_sha = SkyCoord(f"{sha['col3'][i]}:{sha['col4'][i]}:{sha['col5'][i]} "
                             f"{sha['col6'][i]}:{sha['col7'][i]}:{sha['col8'][i]}",
                            unit = (u.hourangle, u.deg))
        merged.writelines(f"{cnt} {coord_sha.ra.value} {coord_sha.dec.value} {sha_z[smem][i]} {sha_ze[smem][i]} \n")
        # merged.writelines(f"{coord_sha.ra.value} {coord_sha.dec.value} \n")
    snum = cnt

    cat = ascii.read(fn + 'KYDISC/Oh+2019_mrt.txt')
    coords_kyd = SkyCoord(cat['RAdeg'], cat['DEdeg'], unit='deg')
    kmem = (coords_kyd.separation(ab.coords_cl_cen[k]).value < r_ran) & \
           (cat['z'] > ab.redshifts[4] - z_ran) & (cat['z'] < ab.redshifts[4] + z_ran)
    cat = cat[kmem]

    for j in range(0, len(cat)):
        if cat['ID'][j].split('_')[0] == 'A3574':
            coord_kyd = SkyCoord(cat['RAdeg'][j], cat['DEdeg'][j], unit='deg')
            if not np.any(coords_kyd.separation(ab.coords_cl_cen[k]).value < r_ran):
                cnt = cnt + 1
                merged.writelines(f"{cnt} {cat['RAdeg'][j]} {cat['DEdeg'][j]} {cat['z'][j]} {cat['e_z'][j]} \n")

    knum = cnt - snum
    print(f'A3574 n_k_mem={knum}, n_s_mem={snum}, tot={cnt}')


### KYDISC
fn = "/Users/duhokim/work/abell/spec/"
with open(fn+f'A3659_spec_ra_dec_z_zran_{z_ran}_rran_{r_ran}.txt', 'w') as merged:
    k = 5
    cnt = 0
    cat = ascii.read(fn+'KYDISC/Oh+2019_mrt.txt')
    coords_kyd = SkyCoord(cat['RAdeg'], cat['DEdeg'], unit='deg')
    kmem = (coords_kyd.separation(ab.coords_cl_cen[k]).value < r_ran) & \
           (cat['z'] > ab.redshifts[k] - z_ran) & (cat['z'] < ab.redshifts[k] + z_ran)
    cat = cat[kmem]
    for j in range(0, len(cat)):
        if cat['ID'][j].split('_')[0] == 'A3659':
            cnt = cnt + 1
            # coord_sha = SkyCoord(f"{cat['col2'][j]}:{cat['col3'][j]}:{cat['col4'][j]} "
            #                      f"{cat['col5'][j]}:{cat['col6'][j]}:{cat['col7'][j]}",
            #                     unit = (u.hourangle, u.deg))
            merged.writelines(f"{cnt} {cat['RAdeg'][j]} {cat['DEdeg'][j]} {cat['z'][j]} {cat['e_z'][j]} \n")
            # merged.writelines(f"{cat['col10'][j]} {cat['col11'][j]} \n")
    print(f'A3659 n_k_mem={cnt}, tot={cnt}')

###  A754  ####
fn = "/Users/duhokim/work/abell/spec/"
with open(fn + f'A754_spec_ra_dec_z_zran_{z_ran}_rran_{r_ran}.txt', 'w') as merged:
    k = 0
    cnt = 0
    wings = ascii.read(fn + 'WINGS/WINGS_Cava+2009_A754.txt')
    omega = ascii.read(fn + 'WINGS/OmegaWINGS_Moretti+2017_A754.txt')

    coords_wings = SkyCoord(wings['col2'], wings['col3'], unit='deg')
    coords_omega = SkyCoord(omega['col4'], omega['col5'], unit='deg')

    # only members
    wmem = (wings['col8'] == 1) & (coords_wings.separation(ab.coords_cl_cen[k]).value < r_ran) & \
           (wings['col5'] / const.c * 1e3 > ab.redshifts[k] - z_ran) & \
           (wings['col5'] / const.c * 1e3 < ab.redshifts[k] + z_ran)
    omem = (omega['col10'] == 1)& (coords_omega.separation(ab.coords_cl_cen[k]).value < r_ran) & \
           (omega['col6'] > ab.redshifts[k] - z_ran) & \
           (omega['col6'] < ab.redshifts[k] + z_ran)

    idx_o2w, d2d_w2o, d3d = coords_wings[wmem].match_to_catalog_sky(coords_omega[omem])
    not_match_w2o = (d2d_w2o.arcsec > 1)

    for i in range(0, len(wings[wmem][not_match_w2o])):
        cnt = cnt + 1
        z = wings['col5'][wmem][not_match_w2o][i]/const.c*1e3
        merged.writelines(f"{cnt} {coords_wings[wmem][not_match_w2o][i].ra.value} "
                          f"{coords_wings[wmem][not_match_w2o][i].dec.value} "
                          f"{z} {wings['col6'][wmem][not_match_w2o][i]/const.c*1e3} \n")

    for i in range(0, len(omega[omem])):
        cnt = cnt + 1
        z = omega['col6'][omem][i]
        merged.writelines(f"{cnt} {coords_omega[omem][i].ra.value} "
                          f"{coords_omega[omem][i].dec.value} {z} "
                          f"{omega['col7'][omem][i]} \n")
    print(f'A754 n_w={len(wmem)}, n_w_mem={sum(wmem)}, n_o={len(omem)}, n_o_mem={sum(omem)}, tot={cnt}')


###  A2399 for comparison with Ana's ####
# fn = "/Users/duhokim/work/abell/spec/"
# with open(fn + 'A2399_spec_ra_dec_z_ana_dr12_2nd.txt', 'w') as merged:
#     cnt = 0
#     wings = ascii.read(fn + 'WINGS/WINGS_Cava+2009_A2399.txt')
#     omega = ascii.read(fn + 'WINGS/OmegaWINGS_Moretti+2017_A2399.txt')
#     sdss = ascii.read(fn + 'SDSS/A2399_dr12.csv')
#
#     # Ana's cl_cen
#     coord_clcen = SkyCoord(329.37258, -7.79569, unit='deg')
#     coords_wings = SkyCoord(wings['col2'], wings['col3'], unit='deg')
#     coords_omega = SkyCoord(omega['col4'], omega['col5'], unit='deg')
#     coords_sdss = SkyCoord(sdss['ra'], sdss['dec'], unit='deg')
#     wings_cen_dist = coords_wings.separation(coord_clcen)
#     omega_cen_dist = coords_omega.separation(coord_clcen)
#     sdss_cen_dist = coords_sdss.separation(coord_clcen)
#
#     # only members
#     wmem = (wings['col5']/const.c*1e3 > 0.05) & (wings['col5']/const.c*1e3 < 0.07) & \
#            (wings_cen_dist.arcmin < 37.89) & (wings['col8'] == 1)
#     omem = (omega['col6'] > 0.05) & (omega['col6'] < 0.07) & \
#            (omega_cen_dist.arcmin < 37.89) & (omega['col10'] == 1)
#     smem = (sdss['z'] > 0.05) & (sdss['z'] < 0.07) & (sdss_cen_dist.arcmin < 37.89) & \
#            (sdss['zWarning'] == 0) & (sdss['clean'] == 1)
#
#     idx_o2w, d2d_w2o, d3d = coords_wings[wmem].match_to_catalog_sky(coords_omega[omem])
#     not_match_w2o = (d2d_w2o.arcsec > 1)
#
#     idx_w2s, d2d_s2w, d3d = coords_sdss[smem].match_to_catalog_sky(coords_wings)
#     idx_o2s, d2d_s2o, d3d = coords_sdss[smem].match_to_catalog_sky(coords_omega)
#     not_match_s2ow = (d2d_s2w.arcsec > 1.0) & (d2d_s2o.arcsec > 1.0)
#
#     for i in range(0, len(wings[wmem][not_match_w2o])):
#         cnt = cnt + 1
#         z = wings['col5'][wmem][not_match_w2o][i]/const.c*1e3
#         # dist = mm.clcendist(ab.redshifts[0], z, ab.coords_cl[0], coord)
#         merged.writelines(f"{cnt} {coords_wings[wmem][not_match_w2o][i].ra.value} "
#                           f"{coords_wings[wmem][not_match_w2o][i].dec.value} "
#                           f"{z} {wings['col6'][wmem][not_match_w2o][i]/const.c*1e3} \n")
#         # merged.writelines(f"{coord.ra.value} {coord.dec.value} \n")
#
#     for i in range(0, len(omega[omem])):
#         cnt = cnt + 1
#         z = omega['col6'][omem][i]
#         # dist = mm.clcendist(ab.redshifts[0], z, ab.coords_cl[0], coord)
#         merged.writelines(f"{cnt} {coords_omega[omem][i].ra.value} "
#                           f"{coords_omega[omem][i].dec.value} {z} "
#                           f"{omega['col7'][omem][i]} \n")
#         # merged.writelines(f"{coord.ra.value} {coord.dec.value} \n")
#
#     for i in range(0, len(sdss[smem][not_match_s2ow])):
#         cnt = cnt + 1
#         merged.writelines(f"{cnt} {coords_sdss[smem][not_match_s2ow][i].ra.value} "
#                           f"{coords_sdss[smem][not_match_s2ow][i].dec.value} "
#                           f"{sdss['z'][smem][not_match_s2ow][i]} "
#                           f"{sdss['zErr'][smem][not_match_s2ow][i]} \n")
#
#     print(f'A2399_Ana n_w={len(wmem)}, n_w_mem={sum(wmem)}, n_o={len(omem)}, n_o_mem={sum(omem)}, '
#           f'n_s={len(smem)}, n_s_mem={sum(smem)}, tot={cnt}')

####  A2399  ####
with open(fn + f'A2399_spec_ra_dec_z_zran_{z_ran}_rran_{r_ran}.txt', 'w') as merged:
    k = 1
    cnt = 0
    wings = ascii.read(fn + 'WINGS/WINGS_Cava+2009_A2399.txt')
    omega = ascii.read(fn + 'WINGS/OmegaWINGS_Moretti+2017_A2399.txt')
    sdss = ascii.read(fn + 'SDSS/A2399_dr12.csv')

    coords_wings = SkyCoord(wings['col2'], wings['col3'], unit='deg')
    coords_omega = SkyCoord(omega['col4'], omega['col5'], unit='deg')
    coords_sdss = SkyCoord(sdss['ra'], sdss['dec'], unit='deg')

    # only members
    wmem = (wings['col8'] == 1) & (coords_wings.separation(ab.coords_cl_cen[k]).value < r_ran) & \
           (wings['col5'] / const.c * 1e3 > ab.redshifts[k] - z_ran) & \
           (wings['col5'] / const.c * 1e3 < ab.redshifts[k] + z_ran)
    omem = (omega['col10'] == 1) & (coords_omega.separation(ab.coords_cl_cen[k]).value < r_ran) & \
           (omega['col6'] > ab.redshifts[k] - z_ran) & \
           (omega['col6'] < ab.redshifts[k] + z_ran)
    smem = (coords_sdss.separation(ab.coords_cl_cen[k]).value < r_ran) & \
           (sdss['z'] > ab.redshifts[k] - z_ran) & (sdss['z'] < ab.redshifts[k] + z_ran) & \
           (sdss['zWarning'] == 0) & (sdss['clean'] == 1)

    idx_o2w, d2d_w2o, d3d = coords_wings[wmem].match_to_catalog_sky(coords_omega[omem])
    not_match_w2o = (d2d_w2o.arcsec > 1)

    idx_w2s, d2d_s2w, d3d = coords_sdss[smem].match_to_catalog_sky(coords_wings)
    idx_o2s, d2d_s2o, d3d = coords_sdss[smem].match_to_catalog_sky(coords_omega)
    not_match_s2ow = (d2d_s2w.arcsec > 1.0) & (d2d_s2o.arcsec > 1.0)

    for i in range(0, len(wings[wmem][not_match_w2o])):
        cnt = cnt + 1
        z = wings['col5'][wmem][not_match_w2o][i]/const.c*1e3
        # dist = mm.clcendist(ab.redshifts[0], z, ab.coords_cl[0], coord)
        merged.writelines(f"{cnt} {coords_wings[wmem][not_match_w2o][i].ra.value} "
                          f"{coords_wings[wmem][not_match_w2o][i].dec.value} "
                          f"{z} {wings['col6'][wmem][not_match_w2o][i]/const.c*1e3} \n")
        # merged.writelines(f"{coord.ra.value} {coord.dec.value} \n")

    for i in range(0, len(omega[omem])):
        cnt = cnt + 1
        z = omega['col6'][omem][i]
        # dist = mm.clcendist(ab.redshifts[0], z, ab.coords_cl[0], coord)
        merged.writelines(f"{cnt} {coords_omega[omem][i].ra.value} "
                          f"{coords_omega[omem][i].dec.value} {z} "
                          f"{omega['col7'][omem][i]} \n")
        # merged.writelines(f"{coord.ra.value} {coord.dec.value} \n")

    for i in range(0, len(sdss[smem][not_match_s2ow])):
        cnt = cnt + 1
        merged.writelines(f"{cnt} {coords_sdss[smem][not_match_s2ow][i].ra.value} "
                          f"{coords_sdss[smem][not_match_s2ow][i].dec.value} "
                          f"{sdss['z'][smem][not_match_s2ow][i]} "
                          f"{sdss['zErr'][smem][not_match_s2ow][i]} \n")

    print(f'A2399 n_w={len(wmem)}, n_w_mem={sum(wmem)}, n_o={len(omem)}, n_o_mem={sum(omem)}, '
          f'n_s={len(smem)}, n_s_mem={sum(smem)}, tot={cnt}')

###  A2670  ####
fn = "/Users/duhokim/work/abell/spec/"
with open(fn + f'A2670_spec_ra_dec_z_zran_{z_ran}_rran_{r_ran}.txt', 'w') as merged:
    k = 2
    cnt = 0
    sdss = ascii.read(fn + 'SDSS/A2670_dr12.csv')

    # cl_cen
    coords_sdss = SkyCoord(sdss['ra'], sdss['dec'], unit='deg')

    # only members
    smem = (coords_sdss.separation(ab.coords_cl_cen[k]).value < r_ran) & (sdss['zWarning'] == 0) & (sdss['clean'] == 1) & \
           (sdss['z'] > ab.redshifts[k] - z_ran) & (sdss['z'] < ab.redshifts[k] + z_ran)

    for i in range(0, len(sdss[smem])):
        cnt = cnt + 1
        merged.writelines(f"{cnt} {coords_sdss[smem][i].ra.value} "
                          f"{coords_sdss[smem][i].dec.value} "
                          f"{sdss['z'][smem][i]} "
                          f"{sdss['zErr'][smem][i]} \n")

    print(f'A2670 n_s={len(smem)}, n_s_mem={sum(smem)}, tot={cnt}')

####   A3558   #####
fn = "/Users/duhokim/work/abell/spec/"
with open(fn + f'A3558_spec_ra_dec_z_zran_{z_ran}_rran_{r_ran}.txt', 'w') as merged,   \
    open(fn + f'A3558_spec_red_blue_zran_{z_ran}_rran_{r_ran}.reg', 'w') as reg:
    k = 3
    cnt = 0
    omega = ascii.read(fn + 'WINGS/OmegaWINGS_Moretti+2017_A3558.txt')
    sha = ascii.read(fn + 'Shapley/catalog.dat')

    # cl_cen
    coords_omega = SkyCoord(omega['col4'], omega['col5'], unit='deg')
    coords_sha_str = []
    for i in range(len(sha)):
        coords_sha_str.append(f"{sha['col3'][i]}:{sha['col4'][i]}:{sha['col5'][i]} "
                         f"{sha['col6'][i]}:{sha['col7'][i]}:{sha['col8'][i]}")
    coords_sha = SkyCoord(coords_sha_str, unit=(u.hourangle, u.deg))

    # only members
    omem = (omega['col10'] == 1) & (coords_omega.separation(ab.coords_cl_cen[k]).value < r_ran) & \
           (omega['col6'] > ab.redshifts[k] - z_ran) & \
           (omega['col6'] < ab.redshifts[k] + z_ran)


    sha_z = sha['col18'] / const.c * 1e3
    sha_ze = sha['col19'] / const.c * 1e3
    smem = (coords_sha.separation(ab.coords_cl_cen[k]).value < r_ran) & \
           (sha_z > ab.redshifts[k] - z_ran) & \
           (sha_z < ab.redshifts[k] + z_ran)

    # cross-match
    idx_o2s, d2d_s2o, d3d = coords_sha[smem].match_to_catalog_sky(coords_omega[omem])
    not_match_s2o = (d2d_s2o.arcsec > 1)

    for i in range(0, len(omega[omem])):
        cnt = cnt + 1
        merged.writelines(f"{cnt} {coords_omega[omem][i].ra.value} "
                          f"{coords_omega[omem][i].dec.value} {omega['col6'][omem][i]} "
                          f"{omega['col7'][omem][i]} \n")
        if omega['col6'][omem][i] > 0.048:
            reg.writelines(f"j2000; circle({coords_omega[omem][i].ra.value}, {coords_omega[omem][i].dec.value}, "
                                                  f"10\") # width=3, color='red' \n")
        else:
            reg.writelines(f"j2000; circle({coords_omega[omem][i].ra.value}, {coords_omega[omem][i].dec.value}, "
                           f"10\") # width=3, color='blue' \n")

    for i in range(0, len(sha[smem][not_match_s2o])):
        cnt = cnt + 1
        merged.writelines(f"{cnt} {coords_sha[smem][not_match_s2o][i].ra.value} "
                          f"{coords_sha[smem][not_match_s2o][i].dec.value} "
                          f"{sha_z[smem][not_match_s2o][i]} "
                          f"{sha_ze[smem][not_match_s2o][i]} \n")
        if sha_z[smem][not_match_s2o][i] > 0.048:
            reg.writelines(f"j2000; circle({coords_sha[smem][not_match_s2o][i].ra.value}, "
                           f"{coords_sha[smem][not_match_s2o][i].dec.value}, "
                           f"10\") # width=3, color='red' \n")
        else:
            reg.writelines(f"j2000; circle({coords_sha[smem][not_match_s2o][i].ra.value}, "
                           f"{coords_sha[smem][not_match_s2o][i].dec.value}, "
                           f"10\") # width=3, color='blue' \n")
    print(f'A3558 n_o={len(omem)}, n_o_mem={sum(omem)}, n_s={len(smem)}, n_s_mem={sum(smem)}, tot={cnt}')

####  A3716  ######
fn = "/Users/duhokim/work/abell/spec/"
with open(fn + f'A3716_spec_ra_dec_z_zran_{z_ran}_rran_{r_ran}.txt', 'w') as merged:
    k = 6
    cnt = 0
    omega = ascii.read(fn + 'WINGS/OmegaWINGS_Moretti+2017_A3716.txt')
    coords_omega = SkyCoord(omega['col4'], omega['col5'], unit='deg')
    # only members
    omem = (omega['col10'] == 1) & (coords_omega.separation(ab.coords_cl_cen[k]).value < r_ran) & \
           (omega['col6'] > ab.redshifts[k] - z_ran) & \
           (omega['col6'] < ab.redshifts[k] + z_ran)

    for i in range(0, len(omega[omem])):
        cnt = cnt + 1
        merged.writelines(f"{cnt} {coords_omega[omem][i].ra.value} "
                          f"{coords_omega[omem][i].dec.value} {omega['col6'][omem][i]} "
                          f"{omega['col7'][omem][i]} \n")

    print(f'A3716 n_o={len(omem)}, n_o_mem={sum(omem)}, tot={cnt}')

# fn = "/Users/duhokim/work/abell/spec/SDSS/A2670"
# with open(fn+'.reg', 'w') as reg:
#     cat = ascii.read(fn+'.csv')
#     for j in range(0, len(cat)):
#         if cat['z'][j] > 0.05 and cat['z'][j] < 0.1:
#             color = 'green'
#         else:
#             color = 'red'
#         reg.writelines(f"j2000; circle({cat['ra'][j]}, {cat['dec'][j]}, 10\") # color={color}\n")

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

# ###  OmegaWINGS
# fn = "/Users/duhokim/work/abell/spec/"
# with open(fn + 'A3716_spec.radec', 'w') as merged:
#     cnt = 0
#     omega = ascii.read(fn + 'WINGS/OmegaWINGS_Moretti+2017_A3716.txt')
#     for i in range(0, len(omega)):
#         if omega['col10'][i]:
#             cnt = cnt + 1
#             # merged.writelines(f"{cnt} {omega['col4'][i]} {omega['col5'][i]} {omega['col6'][i]} {omega['col7'][i]} \n")
#             merged.writelines(f"{omega['col4'][i]} {omega['col5'][i]} \n")

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