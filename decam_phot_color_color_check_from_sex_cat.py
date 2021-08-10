from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.coordinates import SkyCoord
import abell_cluster_module as ab

size = 1
alpha = 0.05

tag = 'deblend'

params = ['MEAN_MAG_DIFF', 'DATE', 'EXP_TIME', 'SEEING', 'AIRMASS', 'TOTAL_MATCH_NUM']
num_of_param = len(params)

am_cat = pd.read_csv('/Users/duhokim/work/abell/sex/cat/airmass.csv')

fig, axs = plt.subplots(4, len(ab.clusters), tight_layout=True, figsize=(20, 12))
# fig2, axs2 = plt.subplots(2, 3, tight_layout=True, figsize=(12, 8))

for k in range(0, len(ab.clusters)):
# for k in range(0, 1):
    for i in range(0, 2):
        if i == 0:
            sex_cat_u = ascii.read(ab.short_cat_dir + ab.short_cat_fn[k][0] + '_' + tag + '.cat')
            sex_cat_g = ascii.read(ab.short_cat_dir + ab.short_cat_fn[k][1] + '_' + tag + '.cat')
            sex_cat_r = ascii.read(ab.short_cat_dir + ab.short_cat_fn[k][2] + '_' + tag + '.cat')

            # standardize magnitudes and Milky Way extinction correction
            sex_cat_u[ab.mag_sys] = sex_cat_u[ab.mag_sys] + ab.short_a[k][0] + ab.short_b[k][0] - ab.mw_ext[k][0]
            sex_cat_g[ab.mag_sys] = sex_cat_g[ab.mag_sys] + ab.short_a[k][1] + ab.short_b[k][1] - ab.mw_ext[k][1]
            sex_cat_r[ab.mag_sys] = sex_cat_r[ab.mag_sys] + ab.short_a[k][2] + ab.short_b[k][2] - ab.mw_ext[k][2]
        else:
            sex_cat_u = ascii.read(ab.work_dir + 'sex/cat/stack/' + ab.clusters[k] + '_usi.cat')
            sex_cat_g = ascii.read(ab.work_dir + 'sex/cat/stack/' + ab.clusters[k] + '_gsi.cat')
            sex_cat_r = ascii.read(ab.work_dir + 'sex/cat/stack/' + ab.clusters[k] + '_rsi.cat')

            sex_cat_u[ab.mag_sys] = sex_cat_u[ab.mag_sys] + ab.stack_a[k][0] #- ab.mw_ext[k][0]
            sex_cat_g[ab.mag_sys] = sex_cat_g[ab.mag_sys] + ab.stack_a[k][1] #- ab.mw_ext[k][1]
            sex_cat_r[ab.mag_sys] = sex_cat_r[ab.mag_sys] + ab.stack_a[k][2] #- ab.mw_ext[k][2]

        ### match r to u and g band short exposure ###
        coords_u = SkyCoord(sex_cat_u['ALPHA_J2000'], sex_cat_u['DELTA_J2000'], unit='deg')
        coords_g = SkyCoord(sex_cat_g['ALPHA_J2000'], sex_cat_g['DELTA_J2000'], unit='deg')
        coords_r = SkyCoord(sex_cat_r['ALPHA_J2000'], sex_cat_r['DELTA_J2000'], unit='deg')

        idx_u2r, d2d_r2u, d3d = coords_r.match_to_catalog_sky(coords_u)
        idx_g2r, d2d_r2g, d3d = coords_r.match_to_catalog_sky(coords_g)

        # bright_central_galaxies = (merged_r['CLASS_STAR'] < ab.class_star_lim) & (merged_r[ab.mag_sys] < ab.mag_lim) & \
        #                           (coords_merged_r.separation(ab.coords_cl_cen[k]).value < ab.rad_lim)

        match_r2u = (d2d_r2u.arcsec < ab.max_sep) & (coords_r.separation(ab.coords_cl_cen[k]).value < ab.rad_lim)
        match_r2g = (d2d_r2g.arcsec < ab.max_sep) & (coords_r.separation(ab.coords_cl_cen[k]).value < ab.rad_lim)

        axs[0 + i*2, k].scatter(sex_cat_r[match_r2g][ab.mag_sys],
                          sex_cat_g[idx_g2r[match_r2g]][ab.mag_sys] - sex_cat_r[match_r2g][ab.mag_sys],
                          alpha=alpha, s=size)
        axs[1 + i*2, k].scatter(sex_cat_r[match_r2u][ab.mag_sys],
                          sex_cat_u[idx_u2r[match_r2u]][ab.mag_sys] - sex_cat_r[match_r2u][ab.mag_sys],
                          alpha=alpha, s=size)

        axs[0 + i*2, k].set_title(ab.clusters[k])
        axs[0 + i*2, k].set_xlabel('r')
        axs[0 + i*2, k].set_xlabel('r')
        axs[0 + i*2, k].set_ylabel('g - r')
        axs[1 + i*2, k].set_ylabel('u - r')
        axs[0 + i*2, k].set_xlim([13, 20])
        axs[1 + i*2, k].set_xlim([13, 20])
        axs[0 + i*2, k].set_ylim([-1, 3])
        axs[1 + i*2, k].set_ylim([0, 6])

fig.savefig(ab.plot_dir + 'CMD_from_sex_cat_'+ab.ver+'_'+tag+'_test.png')

