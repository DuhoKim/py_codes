from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import abell_cluster_module as ab
import statsmodels.api as sm
from astropy.coordinates import SkyCoord  # High-level coordinates
import astropy.units as u

ver = '25rmag'

size = 15
alpha = 0.5

size2 = 30
alpha2 = 0.15

xran = [13, 19]
g_yran = [0, 1.5]
u_yran = [1, 4]


petals_coords = [SkyCoord('09:08:53', '-09:10:57', unit=(u.hourangle, u.deg)),
                 SkyCoord('09:11:12', '-09:40:48', unit=(u.hourangle, u.deg)),
                 SkyCoord('09:06:17', '-09:42:06', unit=(u.hourangle, u.deg)),
                 SkyCoord('09:08:57', '-09:08:28', unit=(u.hourangle, u.deg))]

cols = ['black', 'red', 'purple', 'darkorange', 'olive', 'darkgreen', 'teal']

params = ['MEAN_MAG_DIFF', 'DATE', 'EXP_TIME', 'SEEING', 'AIRMASS', 'TOTAL_MATCH_NUM']
num_of_param = len(params)

am_cat = pd.read_csv('/Users/duhokim/work/abell/sex/cat/airmass.csv')

mmt_dir = '/Users/duhokim/work/abell/proposal/MMT/'

fig = plt.figure()
# for k in range(0, len(ab.clusters)):
for k in range(0, 1):
    sex_cat = ascii.read(ab.sex_dir+f'DECam_merged_SEx_cat_{ab.clusters[k]}_Gal_ext_corrected_{ab.ver}_dqm_edit.txt')
    sex_coords = SkyCoord(sex_cat['ALPHA_J2000'], sex_cat['DELTA_J2000'], unit='deg')

    in_petals = [(x.separation(petals_coords[0]).deg < 0.5) | (x.separation(petals_coords[1]).deg < 0.5) |
                 (x.separation(petals_coords[2]).deg < 0.5) | (x.separation(petals_coords[3]).deg < 0.5)
                 for x in sex_coords]

    in_petals_N = [x.separation(petals_coords[0]).deg < 0.5 for x in sex_coords]
    in_petals_E = [x.separation(petals_coords[1]).deg < 0.5 for x in sex_coords]
    in_petals_W = [x.separation(petals_coords[2]).deg < 0.5 for x in sex_coords]
    in_petals_S = [x.separation(petals_coords[3]).deg < 0.5 for x in sex_coords]

    is_16 = [x[ab.mag_sys+ '_r' ] < 16 for x in sex_cat]
    is_18 = [x[ab.mag_sys + '_r'] < 18 for x in sex_cat]
    is_19 = [x[ab.mag_sys + '_r'] < 19 for x in sex_cat]

    # sep0 = sex_coords.separation(petals_coords[0])
    # sep1 = sex_coords.separation(petals_coords[1])
    # sep2 = sex_coords.separation(petals_coords[2])
    # sep3 = sex_coords.separation(petals_coords[3])
    #
    # in_petals = (sep0.deg < 0.5) | (sep1.deg < 0.5) | sep2.deg < 0.5 or sep3.deg < 0.5

    # sex_cat = sex_cat[in_petals]

    cat = ascii.read("/Users/duhokim/work/abell/spec/WINGS/WINGS_Cava+2009_A754.txt")

    # good_g = (sex_cat['CLASS_STAR'] < 0.2) & (sex_cat[ab.mag_sys+'_g'] < 30) & (sex_cat[ab.mag_sys+'_r'] < xran[1]) & (sex_cat[ab.magerr_sys+'_r'] > 0)

    plt.scatter(sex_cat[ab.mag_sys+'_r'],
                      sex_cat[ab.mag_sys+'_g'] - sex_cat[ab.mag_sys+'_r'],
                      alpha=0.2, s=size)

    red_seq = (sex_cat[ab.mag_sys+'_g'] - sex_cat[ab.mag_sys+'_r'] < -0.02 * sex_cat[ab.mag_sys+'_r'] + 1.26) & \
              (sex_cat[ab.mag_sys+'_g'] - sex_cat[ab.mag_sys+'_r'] > -0.02 * sex_cat[ab.mag_sys+'_r'] + 1.01) & \
             (sex_cat[ab.magerr_sys + '_r'] > 0) & \
              (sex_cat[ab.mag_sys + '_r'] < 19)

    g_r = np.array(sex_cat[ab.mag_sys + '_g'][red_seq] - sex_cat[ab.mag_sys + '_r'][red_seq])
    r_good_g = np.array(sex_cat[ab.mag_sys + '_r'][red_seq])
    merr_good_g = np.array(sex_cat[ab.magerr_sys + '_r'][red_seq])

    coords = []
    for i in range(0, len(cat)):
        coords.append(f"{cat['col8'][i]}:{cat['col9'][i]}:{cat['col10'][i]} {cat['col11'][i]}:{cat['col12'][i]}:{cat['col13'][i]}")
    coords_cat = SkyCoord(coords, unit=(u.hourangle, u.deg))
    coords_sex = SkyCoord(sex_cat['ALPHA_J2000'], sex_cat['DELTA_J2000'], unit='deg')
    idx_cat2sex, d2d, d3d = coords_sex.match_to_catalog_sky(coords_cat)
    matched_sex = (d2d.arcsec < 1.0)

    plt.scatter(sex_cat[ab.mag_sys+'_r'][matched_sex],
                      sex_cat[ab.mag_sys+'_g'][matched_sex] - sex_cat[ab.mag_sys+'_r'][matched_sex],
                      alpha=0.2, s=30, color='red')

    plt.plot(xran, [1, 0.9], '--', alpha=0.8, color='teal')
    plt.plot(xran, [0.75, 0.65], '--', alpha=0.8, color='teal')
    plt.fill_between(xran, [1, 0.9], [0.75, 0.65], alpha=0.1, color='teal')
    plt.axvline(x=16)
    plt.axvline(x=18)
    # plt.axvline(x=19)

    # plt.set_title(ab.clusters[k])
    plt.xlim(xran)
    plt.ylim([0, 1.5])
    plt.xlabel('r', fontsize=20)
    plt.ylabel('g - r', fontsize=20)
    plt.text(13.2, 1.4, 'A754-N:', fontsize=15)
    # plt.text(13.2, 1.4, r'', fontsize=15)
    plt.text(15, 1.4, f'{np.sum([a and b and c for a, b, c in zip(in_petals_N, is_16, red_seq)])} '
                        f'({np.sum([a and b and c and d for a, b, c, d in zip(in_petals_N, is_16, red_seq, matched_sex)])}) ', fontsize=15)
    plt.text(16.5, 1.4, f'{np.sum([a and b and c for a, b, c in zip(in_petals_N, is_18, red_seq)])} '
                        f'({np.sum([a and b and c and d for a, b, c, d in zip(in_petals_N, is_18, red_seq, matched_sex)])}) ', fontsize=15)
    plt.text(17.8, 1.4, f'{np.sum([a and b and c for a, b, c in zip(in_petals_N, is_19, red_seq)])} '
                        f'({np.sum([a and b and c and d for a, b, c, d in zip(in_petals_N, is_19, red_seq, matched_sex)])}) ', fontsize=15)

    plt.text(13.2, 1.3, 'A754-E:', fontsize=15)
    plt.text(15, 1.3, f'{np.sum([a and b and c for a, b, c in zip(in_petals_E, is_16, red_seq)])} '
                        f'({np.sum([a and b and c and d for a, b, c, d in zip(in_petals_E, is_16, red_seq, matched_sex)])}) ',
             fontsize=15)
    plt.text(16.5, 1.3, f'{np.sum([a and b and c for a, b, c in zip(in_petals_E, is_18, red_seq)])} '
                        f'({np.sum([a and b and c and d for a, b, c, d in zip(in_petals_E, is_18, red_seq, matched_sex)])}) ',
             fontsize=15)
    plt.text(17.8, 1.3, f'{np.sum([a and b and c for a, b, c in zip(in_petals_E, is_19, red_seq)])} '
                        f'({np.sum([a and b and c and d for a, b, c, d in zip(in_petals_E, is_19, red_seq, matched_sex)])}) ',
             fontsize=15)

    plt.text(13.2, 1.2, 'A754-W:', fontsize=15)
    plt.text(15, 1.2, f'{np.sum([a and b and c for a, b, c in zip(in_petals_W, is_16, red_seq)])} '
                        f'({np.sum([a and b and c and d for a, b, c, d in zip(in_petals_W, is_16, red_seq, matched_sex)])}) ', fontsize=15)
    plt.text(16.5, 1.2, f'{np.sum([a and b and c for a, b, c in zip(in_petals_W, is_18, red_seq)])} '
                        f'({np.sum([a and b and c and d for a, b, c, d in zip(in_petals_W, is_18, red_seq, matched_sex)])}) ', fontsize=15)
    plt.text(17.8, 1.2, f'{np.sum([a and b and c for a, b, c in zip(in_petals_W, is_19, red_seq)])} '
                        f'({np.sum([a and b and c and d for a, b, c, d in zip(in_petals_W, is_19, red_seq, matched_sex)])}) ', fontsize=15)

    plt.text(13.2, 1.1, 'A754-S:', fontsize=15)
    plt.text(15, 1.1, f'{np.sum([a and b and c for a, b, c in zip(in_petals_S, is_16, red_seq)])} '
                        f'({np.sum([a and b and c and d for a, b, c, d in zip(in_petals_S, is_16, red_seq, matched_sex)])}) ', fontsize=15)
    plt.text(16.5, 1.1, f'{np.sum([a and b and c for a, b, c in zip(in_petals_S, is_18, red_seq)])} '
                        f'({np.sum([a and b and c and d for a, b, c, d in zip(in_petals_S, is_18, red_seq, matched_sex)])}) ', fontsize=15)
    plt.text(17.8, 1.1, f'{np.sum([a and b and c for a, b, c in zip(in_petals_S, is_19, red_seq)])} '
                        f'({np.sum([a and b and c and d for a, b, c, d in zip(in_petals_S, is_19, red_seq, matched_sex)])}) ', fontsize=15)

    plt.tick_params(labelsize=15, axis='both')







fig.savefig(mmt_dir + 'CMD.png')
print(sum(red_seq))
print(sum(matched_sex))
