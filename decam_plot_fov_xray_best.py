import numpy as np
import my_module as mm
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from reproject.mosaicking import find_optimal_celestial_wcs, reproject_and_coadd
from reproject import reproject_interp
from astropy.visualization.wcsaxes import SphericalCircle, Quadrangle
from astropy import units as u
import abell_cluster_module as ab
import importlib
importlib.reload(ab)
from astropy.cosmology import Planck18 as Cosmo
import img_scale
# from matplotlib.patches import Rectangle
import scipy.ndimage as nd
from scipy.ndimage.filters import gaussian_filter
from matplotlib.pyplot import contour, show
from astropy.io import ascii
from astropy.coordinates import SkyCoord, Distance
import astropy.wcs
import os
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.patches as patches
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from PIL import Image
from scipy.constants import c, G
import sys
sys.path.append('./MilaDS')
import milaDS
from matplotlib.lines import Line2D
from astropy import constants as const
import glob

def incre_cnt_dsp(feat_type, alloc, vote):
    more_than_two = 2 if vote > 2 else 0
    cnt_feat_sub_dsp[feat_type + more_than_two, int(alloc) + 1] += 1

def phase_region(x, y):
    y = abs(y)
    if y < -4 * x + 2:
        return 5
    elif y < 0.5 and y < -4 * x + 6.0:
        return 4
    elif y < 1.5 and y < -4 * x + 5.0:
        return 3
    elif y < 2.5 and y < -4 * x + 4.0:
        return 2
    elif x > 0.375 and y < -0.571 * x + 2.714:
        return 1
    else:
        return 0

def draw_phase_region(draw_ax, draw_alpha):
    draw_ax.plot([1.375, 1.5], [0.5, 0], ls='--', c='black', alpha=draw_alpha)
    draw_ax.plot([0.375, 1.375], [0.5, 0.5], ls='--', c='black', alpha=draw_alpha)
    draw_ax.plot([0.875, 1.125], [1.5, 0.5], ls='--', c='black', alpha=draw_alpha)
    draw_ax.plot([0.125, 0.875], [1.5, 1.5], ls='--', c='black', alpha=draw_alpha)
    draw_ax.plot([0.375, 0.625], [2.5, 1.5], ls='--', c='black', alpha=draw_alpha)
    draw_ax.plot([0.375, 3.0], [2.5, 1.0], ls='--', c='black', alpha=draw_alpha)
    draw_ax.plot([0, 0.375], [2.5, 2.5], ls='--', c='black', alpha=draw_alpha)
    draw_ax.plot([0, 0.5], [2.0, 0.0], ls='--', c='black', alpha=draw_alpha)
    draw_ax.plot([1.375, 1.5], [-0.5, 0], ls='--', c='black', alpha=draw_alpha)
    draw_ax.plot([0.375, 1.375], [-0.5, -0.5], ls='--', c='black', alpha=draw_alpha)
    draw_ax.plot([0.875, 1.125], [-1.5, -0.5], ls='--', c='black', alpha=draw_alpha)
    draw_ax.plot([0.125, 0.875], [-1.5, -1.5], ls='--', c='black', alpha=draw_alpha)
    draw_ax.plot([0.375, 0.625], [-2.5, -1.5], ls='--', c='black', alpha=draw_alpha)
    draw_ax.plot([0.375, 3.0], [-2.5, -1.0], ls='--', c='black', alpha=draw_alpha)
    draw_ax.plot([0, 0.375], [-2.5, -2.5], ls='--', c='black', alpha=draw_alpha)
    draw_ax.plot([0, 0.5], [-2.0, -0.0], ls='--', c='black', alpha=draw_alpha)
    draw_ax.text(1.5, 2.5, 'out of bound', fontsize=20, c='black', alpha=draw_alpha)
    draw_ax.text(1.8, -0.1, 'A', fontsize=30, c='black', alpha=draw_alpha)
    draw_ax.text(0.2, 1.9, 'B', fontsize=30, c='black', alpha=draw_alpha)
    draw_ax.text(0.5, 0.9, 'C', fontsize=30, c='black', alpha=draw_alpha)
    draw_ax.text(0.8, -0.1, 'D', fontsize=30, c='black', alpha=draw_alpha)
    draw_ax.text(0.1, -0.1, 'E', fontsize=30, c='black', alpha=draw_alpha)
    draw_ax.text(0.2, -1.9, 'B', fontsize=30, c='black', alpha=draw_alpha)
    draw_ax.text(0.5, -0.9, 'C', fontsize=30, c='black', alpha=draw_alpha)

def draw_escape_line(draw_ax, draw_alpha, r200, m200, sig_vr):
    xxx = np.linspace(0, 3, 100)
    yyy = np.zeros(100)
    con = 6
    g_c = (np.log(1 + con) - (con / (1 + con))) ** -1
    for i in range(100):
        s = (np.pi / 2.0) * xxx[i]
        kk = g_c * np.log(1 + con * s) / s
        vesc = np.sqrt((2 * G * m200 * 10 ** 14 * const.M_sun.value * kk) /
                    (r200 * 10 ** 6 * const.pc.value))
        yyy[i] = vesc / np.sqrt(3)  / (sig_vr * 10 ** 3)
    draw_ax.plot(xxx, yyy, ':', c='black', alpha=draw_alpha)
    draw_ax.plot(xxx, -yyy, ':', c='black', alpha=draw_alpha)

is_dsp_first = False

plim_ps = [0.01, 0.05, 0.1, 0.5, 1, 5, 10]

comp_no_main = [[0, 0, 0, 6, 0, 0, 0],     # nth BIC
                [0, 4, 0, 1, 4, 1, 0]]     # component number of main cluster

fov_ra = [ 1.5, 2.3, 2.3, 2.7, 0.0, 0.0, 2.2]
fov_dec = [ 1.5, 2.3, 2.3, 2.3, 0.0, 0.0, 1.4]

scale_x = [180, 219, 224, 170, 0, 0, 176]
scale_y = [194, 229, 235, 215, 0, 0, 197]

colors = ['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'cyan', 'magenta', 'yellow', 'salmon','skyblue',
          'lightgrey', 'darkgrey', 'black']
markers = ['o', 's', 'p', 'P', 'h', '+', 'x', 'D']

fig10, axs10 = plt.subplots(tight_layout=True, figsize=(6, 6))
fig12, axs12 = plt.subplots(tight_layout=True, figsize=(12, 12))
fig122, axs122 = plt.subplots(tight_layout=True, figsize=(6, 6))
fig133, axs133 = plt.subplots(tight_layout=True, figsize=(6, 6))
fig16, axs16 = plt.subplots(tight_layout=True, figsize=(6, 6))

### X-ray peak ###
# coords_x = [SkyCoord(137.3030456, -9.6573548, unit='deg'),
#                  SkyCoord(329.4607702, -7.9094757, unit='deg'),
#                  SkyCoord(358.5533262, -10.41283, unit='deg')]
coords_x = [[SkyCoord(137.3030456, -9.6573548, unit='deg')],
                 [SkyCoord(329.4607702, -7.9094757, unit='deg')],
                 [SkyCoord(358.5533262, -10.41283, unit='deg')],
                 [SkyCoord(201.9768614, -31.5060439, unit='deg'),
                    SkyCoord(202.4251519, -31.6107587, unit='deg'),
                    SkyCoord(202.8692166, -31.7892771, unit='deg'),
                  SkyCoord(203.4053549, -31.6889897, unit='deg')]]

# xray_levels = [[10 ** -6.1, 10 ** -5.9, 10 ** -5.7, 10 ** -5.5],
#                [10 ** -6.1, 10 ** -5.9, 10 ** -5.7, 10 ** -5.5],
#                [10 ** -7.1, 10 ** -6.6, 10 ** -6.1, 10 ** -5.7],
#                [10 ** -6.8, 10 ** -6.3, 10 ** -5.8, 10 ** -5.5],
#                 [10 ** -7.1, 10 ** -6.6, 10 ** -6.1, 10 ** -5.7],
#                 [10 ** -7.1, 10 ** -6.6, 10 ** -6.1, 10 ** -5.7],
#                [10 ** -7.1, 10 ** -6.6, 10 ** -6.1, 10 ** -5.7]]

xray_high_fns = ['Chandra/A754.img', 'XMM/A2399.fits', 'Chandra/A2670.img', 'XMM/A3558.fits', '', '', 'XMM/A3716.fits']
xray_high_sigmas = [5, 0, 5, 5, 5, 5, 30]
xray_low_sigmas = [1.0, 1.2, 1.0, 1.0, 1.0, 1.0, 1.0]
xray_high_levels = [[10 ** -6.2, 10 ** -6.1, 10 ** -6.0, 10 ** -5.9, 10 ** -5.8],
               [12, 20, 30, 40],
               [10 ** -7.5, 10 ** -7.0, 10 ** -6.5, 10 ** -6.0, 10 ** -5.5],
               [10 ** 1.8, 10 ** 2.0, 10 ** 2.15, 10 ** 2.3],
                [10 ** -7.1, 10 ** -6.6, 10 ** -6.1, 10 ** -5.7],
                [10 ** -7.1, 10 ** -6.6, 10 ** -6.1, 10 ** -5.7],
               [17, 34, 50, 100]]
xray_low_levels = [[2, 3, 4, 5],
               [0.5, 0.75, 1],
                   [0.5, 0.75, 1],
               [0.5, 0.75, 1],
                [0.5, 0.75, 1],
                [0.5, 0.75, 1],
               [0.5, 0.75, 1]]
xray_center_box = [ [115, 193, 110, 177],
                    [128, 172, 131, 173],
                    [140, 163, 141, 159],
                    [131, 167, 132, 168],
                    [110, 115, 177, 193],
                    [110, 115, 177, 193],
                    [128, 163, 126, 161]]


### center coordinates for a3558 and a3562
coords_ss = [SkyCoord('13:27:54 -31:30:00', unit=(u.hourangle, u.deg)),
                 SkyCoord('13:33:30 -31:40:00', unit=(u.hourangle, u.deg))]

coords_3716 = [SkyCoord('20:52:00 -52:45:18', unit=(u.hourangle, u.deg)),   # S
               SkyCoord('20:51:57 -52:37:40', unit=(u.hourangle, u.deg))]   # N


### radial velocities and dispersions for a3558 and a3562 in km s^-1
rad_vel_mean_ss = [14463, 14351]
rad_vel_sig_ss = [994, 1188]

rad_vel_mean_3716 = [13390, 14315]
rad_vel_sig_3716 = [965, 651]

cluster_ss = ['A3558', 'A3562']
cluster_3716 = ['A3716S', 'A3716N']

r200_ss = [2.46, 2.94]
m200_ss = [16.2, 19.4]
xx_ss = [0.25, 0.75, 1.5]

r200_3716 = [2.39, 1.61]
m200_3716 = [17.7, 10.6]



num_sub = [] # number of substructure
frac_I = []   # percentage of interacting
frac_PM = []  # percentage of post-merger

frac_I_rad = []
frac_PM_rad = []
weight_rad = []

d_int_tot = []
d_pm_tot = []
d_eith_tot = []

local_peak = []

num_I_in = 0
num_I_out = 0
num_I_mid = 0

# 2D array for number counts of galaxies in each region in a phase-space for each feature type
ncnt_reg = np.array([
    [0, 0, 0, 0, 0, 0],     # no feat
    [0, 0, 0, 0, 0, 0],     # PM
    [0, 0, 0, 0, 0, 0]     # I
])

# same as above but for galaxies selected by at least 3 people not 2 people
ncnt_reg_3 = np.array([
    [0, 0, 0, 0, 0, 0],     # no feat
    [0, 0, 0, 0, 0, 0],     # PM
    [0, 0, 0, 0, 0, 0]     # I
])

# for k in range(0, len(clusters)):
for k in [0, 1, 2, 3, 6]:
# for k in [3]:
# for k in [0, 1, 3]:
#    for plim_p in plim_ps:
    mm.print_time()
    print(f'{ab.clusters[k]}')

    # to print current clusters' number
    if (k == 3 or k == 6):
        ncnt_58 = np.array([
            [0, 0, 0, 0, 0, 0],     # no feat
            [0, 0, 0, 0, 0, 0],     # PM
            [0, 0, 0, 0, 0, 0]     # I
        ])
        ncnt_62 = np.array([
            [0, 0, 0, 0, 0, 0],  # no feat
            [0, 0, 0, 0, 0, 0],  # PM
            [0, 0, 0, 0, 0, 0]  # I
        ])
        ncnt_else = np.array([0, 0, 0]) # no feat, PM, I

        # for galaxies by 3 inspectors
        ncnt_58_3 = np.array([
            [0, 0, 0, 0, 0, 0],     # no feat
            [0, 0, 0, 0, 0, 0],     # PM
            [0, 0, 0, 0, 0, 0]     # I
        ])
        ncnt_62_3 = np.array([
            [0, 0, 0, 0, 0, 0],  # no feat
            [0, 0, 0, 0, 0, 0],  # PM
            [0, 0, 0, 0, 0, 0]  # I
        ])
        ncnt_else_3 = np.array([0, 0, 0]) # no feat, PM, I
    else:
        ncnt_this = np.array([
            [0, 0, 0, 0, 0, 0],  # no feat
            [0, 0, 0, 0, 0, 0],  # PM
            [0, 0, 0, 0, 0, 0]  # I
        ])
        ncnt_this_3 = np.array([
            [0, 0, 0, 0, 0, 0],  # no feat
            [0, 0, 0, 0, 0, 0],  # PM
            [0, 0, 0, 0, 0, 0]  # I
        ])

    vis_cat = ascii.read(
        ab.work_dir + f'sex/cat/DECam_merged_SEx_cat_{ab.clusters[k]}_Gal_ext_corrected_20rmag_psf.txt')
    spec_cat = ascii.read(ab.work_dir + f'spec/{ab.clusters[k]}_spec_ra_dec_z_zran_0.05_rran_1.5.txt')
    res_cat = ascii.read(ab.work_dir + f'vis/combined/{ab.clusters[k]}.txt')

    this_num_sub = 0

    fig = plt.figure(figsize=(10, 10))
    fig_ds = plt.figure(figsize=(10, 10))

    if k == 3:
        fig1_58 = plt.figure(figsize=(10, 10))
        fig1_62 = plt.figure(figsize=(10, 10))
        fig1_58_ds = plt.figure(figsize=(10, 10))
        fig1_62_ds = plt.figure(figsize=(10, 10))
        ax1_58 = fig1_58.add_subplot()
        ax1_62 = fig1_62.add_subplot()
        ax1_58_ds = fig1_58_ds.add_subplot()
        ax1_62_ds = fig1_62_ds.add_subplot()
        ax1_in_58 = ax1_58.inset_axes([0.7, 0.05, 0.25, 0.25])
        ax1_in_62 = ax1_62.inset_axes([0.7, 0.05, 0.25, 0.25])
        ax1_in_58_ds = ax1_58_ds.inset_axes([0.7, 0.05, 0.25, 0.25])
        ax1_in_62_ds = ax1_62_ds.inset_axes([0.7, 0.05, 0.25, 0.25])
    elif k == 6:
        fig1_58 = plt.figure(figsize=(10, 10))
        fig1_62 = plt.figure(figsize=(10, 10))
        fig1_ds = plt.figure(figsize=(10, 10))
        ax1_58 = fig1_58.add_subplot()
        ax1_62 = fig1_62.add_subplot()
        ax1_58_ds = fig1_58_ds.add_subplot()
        ax1_62_ds = fig1_62_ds.add_subplot()
        ax1_in_58 = ax1_58.inset_axes([0.7, 0.05, 0.25, 0.25])
        ax1_in_62 = ax1_62.inset_axes([0.7, 0.05, 0.25, 0.25])
    else:
        fig1 = plt.figure(figsize=(10, 10))
        fig1_ds = plt.figure(figsize=(10, 10))
        ax1 = fig1.add_subplot()
        ax1_ds = fig1_ds.add_subplot()
        ax1_in = ax1.inset_axes([0.7, 0.05, 0.25, 0.25])

    file_name = glob.glob(ab.work_dir + f'spec/{ab.clusters[k]}_{comp_no_main[0][k] + 1}_*_class_uncer_zran_0.05_rran_1.5.csv')
    r_cat = ascii.read(file_name[0], format='csv')
    r_cat = r_cat['classification']
    this_num_sub = max(r_cat)

    num_sub.append(this_num_sub)
    num_I = len(np.where(res_cat['col2'] == 2)[0])
    frac_I.append(num_I / len(res_cat))
    num_PM = len(np.where(res_cat['col2'] == 1)[0])
    frac_PM.append(num_PM / len(res_cat))

    coords_vis = SkyCoord(vis_cat['ALPHA_J2000'], vis_cat['DELTA_J2000'], unit='deg')
    coords_spec = SkyCoord(spec_cat['col2'], spec_cat['col3'], unit='deg')
    idx_vis_spec, d2d_spec, d3d = coords_spec.match_to_catalog_sky(coords_vis)
    matched_spec = (d2d_spec.arcsec < ab.max_sep)

    ra, dec = ab.coords_cl_cen[k].ra.deg, ab.coords_cl_cen[k].dec.deg

    z = ab.redshifts[k]  # redshift
    kpc_arcmin = Cosmo.kpc_proper_per_arcmin(z)  # xxx kpc / arcmin
    arcmin_Mpc = 1e3 / kpc_arcmin  # xxx arcmin / Mpc
    r200_deg = ab.r200[k] * 1e3 / kpc_arcmin / 6e1
    pixel_Mpc = arcmin_Mpc * 60.0 / 45.0    # xxx pixel / Mpc

    ax_list = []
    ax_ds_list = []
    ### draw low-resolution x-ray image with wider FoV
    with fits.open(ab.work_dir + 'fits/rosat/' + ab.clusters[k] + '.fits') as hdu_x_rosat:
        wcs_x_rosat, shape_x = find_optimal_celestial_wcs(hdu_x_rosat)  # has only CompImageHDU files
        ax = fig.add_subplot(projection=wcs_x_rosat)
        ax_ds = fig_ds.add_subplot(projection=wcs_x_rosat)
        data_rosat = hdu_x_rosat[0].data
        data_rosat[xray_center_box[k][0]:xray_center_box[k][1], xray_center_box[k][2]:xray_center_box[k][3]] = np.nan
        data_rosat = gaussian_filter(data_rosat, xray_low_sigmas[k])
        ax.contour(data_rosat, levels=[0.5, 0.75, 1], projection=wcs_x_rosat, alpha=0.1, linewidths=5)
        ax_ds.contour(data_rosat, levels=[0.5, 0.75, 1], projection=wcs_x_rosat, alpha=0.1, linewidths=5)

        x_min, y_min = wcs_x_rosat.all_world2pix(ra - fov_ra[k] / 2, dec - fov_dec[k] / 2, 0)
        x_max, y_max = wcs_x_rosat.all_world2pix(ra + fov_ra[k] / 2, dec + fov_dec[k] / 2, 0)
        ax.plot([scale_x[k], scale_x[k] + pixel_Mpc.value], [scale_y[k], scale_y[k]], linewidth=3, color='purple')
        ax.text(0.75, 0.95, f'1 Mpc at z={z:5.3f}', fontsize=15, transform=ax.transAxes)
        ax_ds.plot([scale_x[k], scale_x[k] + pixel_Mpc.value], [scale_y[k], scale_y[k]], linewidth=3, color='purple')
        ax_ds.text(0.75, 0.95, f'1 Mpc at z={z:5.3f}', fontsize=15, transform=ax_ds.transAxes)

    ### draw high-resolution x-ray image at the center first
    with fits.open(ab.work_dir + 'fits/' + xray_high_fns[k]) as hdu_x:
        wcs_x, shape_x = find_optimal_celestial_wcs(hdu_x)  # has only CompImageHDU files
        data = hdu_x[0].data
        if xray_high_sigmas[k]:
            data = gaussian_filter(data, xray_high_sigmas[k])
        ax.contour(data, levels=xray_high_levels[k], transform=ax.get_transform(wcs_x), alpha=0.4, linewidths=3, cmap='seismic')
        ax_ds.contour(data, levels=xray_high_levels[k], transform=ax_ds.get_transform(wcs_x), alpha=0.4, linewidths=3, cmap='seismic')

    #### DS+ test
    # save substructure membership for all plim_p values
    if k == 3:
        dsp_alloc = np.zeros(len(spec_cat))
        this_num_sub_dsp = np.zeros(2, dtype=int)
        if is_dsp_first:
            for i, j, l in zip([2, 5], ['a3558', 'a3562'], [0, 1]):
                xcoor = (coords_spec[r_cat == i].ra - ab.coords_cl_cen[k].ra).arcmin * kpc_arcmin.value
                ycoor = (coords_spec[r_cat == i].dec - ab.coords_cl_cen[k].dec).arcmin * kpc_arcmin.value
                vlos = c * spec_cat['col4'][r_cat == i]
                data_DSp, data_grs_alloc, summary_DSp_grs = milaDS.DSp_groups(Xcoor=xcoor,
                                                                              Ycoor=ycoor,
                                                                              Vlos=vlos,
                                                                              Zclus=z,
                                                                              cluster_name=j,
                                                                              nsims=1000,
                                                                              Plim_P=0.1)
                np.savez(f'{j}_dsp_result_0.1.npz', data_DSp=data_DSp, data_grs_alloc=data_grs_alloc,
                         summary_DSp_grs=summary_DSp_grs)
                dsp_alloc[r_cat == i] = data_grs_alloc[:, 8]
                this_num_sub_dsp[l] = int(max(dsp_alloc[r_cat == i]) + 2)
            print(f'{this_num_sub_dsp}')
        else:
            for i, j, l in zip([2, 5], ['a3558', 'a3562'], [0, 1]):
                dsp_arrays = np.load(f'{j}_dsp_result_0.1.npz')
                data_DSp = dsp_arrays['data_DSp']
                data_grs_alloc = dsp_arrays['data_grs_alloc']
                summary_DSp_grs = dsp_arrays['summary_DSp_grs']
                dsp_alloc[r_cat == i] = data_grs_alloc[:, 8]
                this_num_sub_dsp[l] = int(max(dsp_alloc[r_cat == i]) + 2)
    else:
        dsp_alloc = np.zeros(len(spec_cat))
        this_num_sub_dsp = 0
        if is_dsp_first:
            xcoor = (coords_spec.ra - ab.coords_cl_cen[k].ra).arcmin * kpc_arcmin.value
            ycoor = (coords_spec.dec - ab.coords_cl_cen[k].dec).arcmin * kpc_arcmin.value
            vlos = c * spec_cat['col4']
            data_DSp, data_grs_alloc, summary_DSp_grs = milaDS.DSp_groups(Xcoor=xcoor,
                                                                              Ycoor=ycoor,
                                                                              Vlos=vlos,
                                                                              Zclus=z,
                                                                              cluster_name=ab.clusters[k],
                                                                              nsims=1000,
                                                                              Plim_P=0.1)
            np.savez(f'{ab.clusters[k]}_dsp_result_0.1.npz', data_DSp=data_DSp, data_grs_alloc=data_grs_alloc,
                     summary_DSp_grs = summary_DSp_grs)
            dsp_alloc = data_grs_alloc[:, 8]
            this_num_sub_dsp = int(max(dsp_alloc) + 2)
        else:
            dsp_arrays = np.load(f'{ab.clusters[k]}_dsp_result_0.1.npz')
            data_DSp = dsp_arrays['data_DSp']
            data_grs_alloc = dsp_arrays['data_grs_alloc']
            summary_DSp_grs = dsp_arrays['summary_DSp_grs']
            dsp_alloc = data_grs_alloc[:, 8]
            this_num_sub_dsp = int(max(dsp_alloc) + 2)

    if k == 3:
        r200_deg_58 = r200_ss[0] * 1e3 / kpc_arcmin / 6e1
        r200_deg_62 = r200_ss[1] * 1e3 / kpc_arcmin / 6e1
        r_in_r200_58 = coords_spec.separation(coords_ss[0]).value / r200_deg_58.value
        r_in_r200_62 = coords_spec.separation(coords_ss[1]).value / r200_deg_62.value
    elif k == 6:
        r200_deg_58 = r200_3716[0] * 1e3 / kpc_arcmin / 6e1
        r200_deg_62 = r200_3716[1] * 1e3 / kpc_arcmin / 6e1
        r_in_r200_58 = coords_spec.separation(coords_3716[0]).value / r200_deg_58.value     # S
        r_in_r200_62 = coords_spec.separation(coords_3716[1]).value / r200_deg_62.value     # N
        r_in_r200 = coords_spec.separation(ab.coords_cl_cen[k]).value / r200_deg.value
    else:
        r_in_r200 = coords_spec.separation(ab.coords_cl_cen[k]).value / r200_deg.value

    ### for record subgroups by mclust
    cnt_feat_sub2 = np.zeros((3, 9))      # number of no feat. in each subs
    cnt_feat_sub3 = np.zeros((3, 9))

    ### for record subgroups by DS+
    if k == 3:
        cnt_feat_sub_dsp = np.zeros((5, sum(this_num_sub_dsp) - 1))      # number of no feat. PM, M PM3, M3 in sub or main
    else:
        cnt_feat_sub_dsp = np.zeros((5, this_num_sub_dsp))  # number of no feat. PM, M PM3, M3 in sub or main

    for i in range(len(r_cat)):
        this_mk = '.'
        this_size = 120
        this_alpha = 0.4
        linewidth = 3
        if k == 3:
            delv58 = (spec_cat['col4'][i] * c / 1e3 - rad_vel_mean_ss[0]) / rad_vel_sig_ss[0]
            delv62 = (spec_cat['col4'][i] * c / 1e3 - rad_vel_mean_ss[1]) / rad_vel_sig_ss[1]
        elif k == 6:
            delv58 = (spec_cat['col4'][i] * c / 1e3 - rad_vel_mean_3716[0]) / rad_vel_sig_3716[0]
            delv62 = (spec_cat['col4'][i] * c / 1e3 - rad_vel_mean_3716[1]) / rad_vel_sig_3716[1]
        delv = (spec_cat['col4'][i] * c / 1e3 - ab.vr[k]) / ab.sig_vr[k]
        coord = SkyCoord(ra=spec_cat['col2'][i], dec=spec_cat['col3'][i], unit='deg')
        sep = ab.coords_cl_cen[k].separation(coord)
        if matched_spec[i]:
            if vis_cat[idx_vis_spec[i]]['NUMBER'] in res_cat['col1']:
                ind, = np.where(res_cat['col1'] == vis_cat[idx_vis_spec[i]]['NUMBER'])
                if res_cat['col2'][ind] == 0:  # classified as no feature
                    this_mk = 'o'
                    this_alpha = 0.4
                    cnt_feat_sub2[0, r_cat[i] - 1] += 1
                    cnt_feat_sub3[0, r_cat[i] - 1] += 1
                    if k == 3:      # subgroup 2 and 5 for a3558
                        if (r_cat[i] == 2):     # Green component in 7th BIC
                            ncnt_reg[0, phase_region(r_in_r200_58[i], delv58)] += 1
                            ncnt_reg_3[0, phase_region(r_in_r200_58[i], delv58)] += 1
                            ncnt_58[0, phase_region(r_in_r200_58[i], delv58)] += 1
                            ncnt_58_3[0, phase_region(r_in_r200_58[i], delv58)] += 1
                            incre_cnt_dsp(0, dsp_alloc[i], 2)
                        elif (r_cat[i] == 5):       # Brown component in 6th BIC
                            ncnt_reg[0, phase_region(r_in_r200_62[i], delv62)] += 1
                            ncnt_reg_3[0, phase_region(r_in_r200_62[i], delv62)] += 1
                            ncnt_62[0, phase_region(r_in_r200_62[i], delv62)] += 1
                            ncnt_62_3[0, phase_region(r_in_r200_62[i], delv62)] += 1
                            if dsp_alloc[i] < 0:
                                incre_cnt_dsp(0, -1, 2)
                            else:
                                incre_cnt_dsp(0, dsp_alloc[i] + this_num_sub_dsp[0] - 1, 2)
                        else:
                            ncnt_else[0] += 1
                            ncnt_else_3[0] += 1
                    elif k == 6:      # subgroup 1 for A3716S
                        if (r_cat[i] == 1):     # Orange component in 1st BIC   S
                            ncnt_reg[0, phase_region(r_in_r200_58[i], delv58)] += 1
                            ncnt_reg_3[0, phase_region(r_in_r200_58[i], delv58)] += 1
                            ncnt_58[0, phase_region(r_in_r200_58[i], delv58)] += 1
                            ncnt_58_3[0, phase_region(r_in_r200_58[i], delv58)] += 1
                        elif (r_cat[i] == 2):       # subgroup 2 for A3716N Green
                            ncnt_reg[0, phase_region(r_in_r200_62[i], delv62)] += 1
                            ncnt_reg_3[0, phase_region(r_in_r200_62[i], delv62)] += 1
                            ncnt_62[0, phase_region(r_in_r200_62[i], delv62)] += 1
                            ncnt_62_3[0, phase_region(r_in_r200_62[i], delv62)] += 1
                        else:
                            ncnt_else[0] += 1
                            ncnt_else_3[0] += 1
                        incre_cnt_dsp(0, dsp_alloc[i], 2)
                    else:
                        ncnt_reg[0, phase_region(r_in_r200[i], delv)] += 1
                        ncnt_reg_3[0, phase_region(r_in_r200[i], delv)] += 1
                        ncnt_this[0, phase_region(r_in_r200[i], delv)] += 1
                        ncnt_this_3[0, phase_region(r_in_r200[i], delv)] += 1
                        incre_cnt_dsp(0, dsp_alloc[i], 2)
                elif res_cat['col2'][ind] == 1:  # classified as PM
                    this_mk = 'P'
                    this_size = 300
                    this_alpha = 0.6
                    cnt_feat_sub2[1, r_cat[i] - 1] += 1
                    if res_cat['col3'][ind] > 2:
                        cnt_feat_sub3[1, r_cat[i] - 1] += 1
                    else:
                        cnt_feat_sub3[0, r_cat[i] - 1] += 1
                    if k == 3:
                        if (r_cat[i] == 2):
                            ncnt_reg[1, phase_region(r_in_r200_58[i], delv58)] += 1
                            ncnt_58[1, phase_region(r_in_r200_58[i], delv58)] += 1
                            incre_cnt_dsp(1, dsp_alloc[i], res_cat['col3'][ind])
                            if res_cat['col3'][ind] > 2:
                                ncnt_reg_3[1, phase_region(r_in_r200_58[i], delv58)] += 1
                                ncnt_58_3[1, phase_region(r_in_r200_58[i], delv58)] += 1
                            else:
                                ncnt_reg_3[0, phase_region(r_in_r200_58[i], delv58)] += 1
                                ncnt_58_3[0, phase_region(r_in_r200_58[i], delv58)] += 1
                        elif (r_cat[i] == 5):       # subgroup 1 and 9 for a3562
                            ncnt_reg[1, phase_region(r_in_r200_62[i], delv62)] += 1
                            ncnt_62[1, phase_region(r_in_r200_62[i], delv62)] += 1
                            if dsp_alloc[i] < 0:
                                incre_cnt_dsp(1, dsp_alloc[i], res_cat['col3'][ind])
                            else:
                                incre_cnt_dsp(1, dsp_alloc[i] + this_num_sub_dsp[0] - 1, res_cat['col3'][ind])
                            if res_cat['col3'][ind] > 2:
                                ncnt_reg_3[1, phase_region(r_in_r200_62[i], delv62)] += 1
                                ncnt_62_3[1, phase_region(r_in_r200_62[i], delv62)] += 1
                            else:
                                ncnt_reg_3[0, phase_region(r_in_r200_62[i], delv62)] += 1
                                ncnt_62_3[0, phase_region(r_in_r200_62[i], delv62)] += 1
                        else:
                            ncnt_else[1] += 1
                    elif k == 6:
                        if (r_cat[i] == 1):
                            ncnt_reg[1, phase_region(r_in_r200_58[i], delv58)] += 1
                            ncnt_58[1, phase_region(r_in_r200_58[i], delv58)] += 1
                            if res_cat['col3'][ind] > 2:
                                ncnt_reg_3[1, phase_region(r_in_r200_58[i], delv58)] += 1
                                ncnt_58_3[1, phase_region(r_in_r200_58[i], delv58)] += 1
                            else:
                                ncnt_reg_3[0, phase_region(r_in_r200_58[i], delv58)] += 1
                                ncnt_58_3[0, phase_region(r_in_r200_58[i], delv58)] += 1
                        elif (r_cat[i] == 2):       # subgroup 2 for A3716N
                            ncnt_reg[1, phase_region(r_in_r200_62[i], delv62)] += 1
                            ncnt_62[1, phase_region(r_in_r200_62[i], delv62)] += 1
                            if res_cat['col3'][ind] > 2:
                                ncnt_reg_3[1, phase_region(r_in_r200_62[i], delv62)] += 1
                                ncnt_62_3[1, phase_region(r_in_r200_62[i], delv62)] += 1
                            else:
                                ncnt_reg_3[0, phase_region(r_in_r200_62[i], delv62)] += 1
                                ncnt_62_3[0, phase_region(r_in_r200_62[i], delv62)] += 1
                        else:
                            ncnt_else[1] += 1
                            if res_cat['col3'][ind] > 2:
                                ncnt_else_3[1] += 1
                            else:
                                ncnt_else_3[0] += 1
                        incre_cnt_dsp(1, dsp_alloc[i], res_cat['col3'][ind])
                    else:
                        ncnt_reg[1, phase_region(r_in_r200[i], delv)] += 1
                        ncnt_this[1, phase_region(r_in_r200[i], delv)] += 1
                        if res_cat['col3'][ind] > 2:
                            ncnt_reg_3[1, phase_region(r_in_r200[i], delv)] += 1
                            ncnt_this_3[1, phase_region(r_in_r200[i], delv)] += 1
                        else:
                            ncnt_reg_3[0, phase_region(r_in_r200[i], delv)] += 1
                            ncnt_this_3[0, phase_region(r_in_r200[i], delv)] += 1
                        incre_cnt_dsp(1, dsp_alloc[i], res_cat['col3'][ind])
                elif res_cat['col2'][ind] == 2:  # classified as I
                    this_mk = '*'
                    this_size = 300
                    this_alpha = 0.9
                    cnt_feat_sub2[2, r_cat[i]-1] += 1
                    if res_cat['col3'][ind] > 2:
                        cnt_feat_sub3[2, r_cat[i] - 1] += 1
                    else:
                        cnt_feat_sub3[0, r_cat[i] - 1] += 1
                    if k == 3:
                        if (r_cat[i] == 2):
                            ncnt_reg[2, phase_region(r_in_r200_58[i], delv58)] += 1
                            ncnt_58[2, phase_region(r_in_r200_58[i], delv58)] += 1
                            incre_cnt_dsp(2, dsp_alloc[i], res_cat['col3'][ind])
                            if res_cat['col3'][ind] > 2:
                                ncnt_reg_3[2, phase_region(r_in_r200_58[i], delv58)] += 1
                                ncnt_58_3[2, phase_region(r_in_r200_58[i], delv58)] += 1
                            else:
                                ncnt_reg_3[0, phase_region(r_in_r200_58[i], delv58)] += 1
                                ncnt_58_3[0, phase_region(r_in_r200_58[i], delv58)] += 1
                        elif (r_cat[i] == 5):       # subgroup 1 and 9 for a3562
                            ncnt_reg[2, phase_region(r_in_r200_62[i], delv62)] += 1
                            ncnt_62[2, phase_region(r_in_r200_62[i], delv62)] += 1
                            if dsp_alloc[i] < 0:
                                incre_cnt_dsp(2, dsp_alloc[i], res_cat['col3'][ind])
                            else:
                                incre_cnt_dsp(2, dsp_alloc[i] + this_num_sub_dsp[0] - 1, res_cat['col3'][ind])
                            if res_cat['col3'][ind] > 2:
                                ncnt_reg_3[2, phase_region(r_in_r200_62[i], delv62)] += 1
                                ncnt_62_3[2, phase_region(r_in_r200_62[i], delv62)] += 1
                            else:
                                ncnt_reg_3[0, phase_region(r_in_r200_62[i], delv62)] += 1
                                ncnt_62_3[0, phase_region(r_in_r200_62[i], delv62)] += 1
                        else:
                            ncnt_else[2] += 1
                            if res_cat['col3'][ind] > 2:
                                ncnt_else_3[2] += 1
                            else:
                                ncnt_else_3[0] += 1
                    elif k == 6:
                        if (r_cat[i] == 1):
                            ncnt_reg[2, phase_region(r_in_r200_58[i], delv58)] += 1
                            ncnt_58[2, phase_region(r_in_r200_58[i], delv58)] += 1
                            if res_cat['col3'][ind] > 2:
                                ncnt_reg_3[2, phase_region(r_in_r200_58[i], delv58)] += 1
                                ncnt_58_3[2, phase_region(r_in_r200_58[i], delv58)] += 1
                            else:
                                ncnt_reg_3[0, phase_region(r_in_r200_58[i], delv58)] += 1
                                ncnt_58_3[0, phase_region(r_in_r200_58[i], delv58)] += 1
                        elif (r_cat[i] == 2):       # subgroup 2 for a3716N
                            ncnt_reg[2, phase_region(r_in_r200_62[i], delv62)] += 1
                            ncnt_62[2, phase_region(r_in_r200_62[i], delv62)] += 1
                            if res_cat['col3'][ind] > 2:
                                ncnt_reg_3[2, phase_region(r_in_r200_62[i], delv62)] += 1
                                ncnt_62_3[2, phase_region(r_in_r200_62[i], delv62)] += 1
                            else:
                                ncnt_reg_3[0, phase_region(r_in_r200_62[i], delv62)] += 1
                                ncnt_62_3[0, phase_region(r_in_r200_62[i], delv62)] += 1
                        else:
                            ncnt_else[2] += 1
                        incre_cnt_dsp(2, dsp_alloc[i], res_cat['col3'][ind])
                    else:
                        ncnt_reg[2, phase_region(r_in_r200[i], delv)] += 1
                        ncnt_this[2, phase_region(r_in_r200[i], delv)] += 1
                        if res_cat['col3'][ind] > 2:
                            ncnt_reg_3[2, phase_region(r_in_r200[i], delv)] += 1
                            ncnt_this_3[2, phase_region(r_in_r200[i], delv)] += 1
                        else:
                            ncnt_reg_3[0, phase_region(r_in_r200[i], delv)] += 1
                            ncnt_this_3[0, phase_region(r_in_r200[i], delv)] += 1
                        incre_cnt_dsp(2, dsp_alloc[i], res_cat['col3'][ind])

        ax.scatter(spec_cat['col2'][i], spec_cat['col3'][i], transform=ax.get_transform('fk5'),
                               marker=this_mk, s=this_size, alpha=this_alpha,
                               linewidth=linewidth, facecolor='none', edgecolor=colors[r_cat[i]])
        if k == 3:
            if ((r_cat[i] - 1) == comp_no_main[1][k]):            # a3558
                ax1_58.scatter(r_in_r200_58[i], delv58, marker=this_mk, s=this_size, alpha=this_alpha,
                            linewidth=linewidth, facecolor='none', edgecolor=colors[r_cat[i]])
                ax1_58_ds.scatter(r_in_r200_58[i], delv58, marker=this_mk, s=this_size,
                               edgecolor=colors[int(dsp_alloc[i]) % len(colors)],
                               alpha=this_alpha, linewidth=linewidth, facecolor='none')
                ax_ds.scatter(spec_cat['col2'][i], spec_cat['col3'][i], transform=ax_ds.get_transform('fk5'),
                              marker=this_mk, s=this_size, alpha=this_alpha, linewidth=linewidth,
                              facecolor='none', edgecolor=colors[int(dsp_alloc[i]) % len(colors)])
            elif ((r_cat[i] - 1) == comp_no_main[1][k+1]):       #  a3562
                ax1_62.scatter(r_in_r200_62[i], delv62, marker=this_mk, s=this_size, alpha=this_alpha,
                               linewidth=linewidth, facecolor='none', edgecolor=colors[r_cat[i]])
                ax1_62_ds.scatter(r_in_r200_62[i], delv62, marker=this_mk, s=this_size,
                               edgecolor=colors[int(dsp_alloc[i]) % len(colors)],
                               alpha=this_alpha, linewidth=linewidth, facecolor='none')
                ax_ds.scatter(spec_cat['col2'][i], spec_cat['col3'][i], transform=ax_ds.get_transform('fk5'),
                              marker=this_mk, s=this_size, alpha=this_alpha, linewidth=linewidth,
                              facecolor='none', edgecolor=colors[int(dsp_alloc[i]) % len(colors)])
        elif k == 6:
            if ((r_cat[i] - 1) == comp_no_main[1][k-1]):
                ax1_58.scatter(r_in_r200_58[i], delv58, marker=this_mk, s=this_size, alpha=this_alpha,
                            linewidth=linewidth, facecolor='none', edgecolor=colors[r_cat[i]])
            elif ((r_cat[i] - 1) == comp_no_main[1][k]):
                ax1_62.scatter(r_in_r200_62[i], delv62, marker=this_mk, s=this_size, alpha=this_alpha,
                               linewidth=linewidth, facecolor='none', edgecolor=colors[r_cat[i]])
            ax1_ds.scatter(r_in_r200[i], delv, marker=this_mk, s=this_size,
                              edgecolor=colors[int(dsp_alloc[i]) % len(colors)],
                              alpha=this_alpha, linewidth=linewidth, facecolor='none')
            ax_ds.scatter(spec_cat['col2'][i], spec_cat['col3'][i], transform=ax_ds.get_transform('fk5'),
                          marker=this_mk, s=this_size, alpha=this_alpha, linewidth=linewidth,
                          facecolor='none', edgecolor=colors[int(dsp_alloc[i]) % len(colors)])
        else:
            ax1.scatter(sep / r200_deg, delv, marker=this_mk, s=this_size, alpha=this_alpha,
                            linewidth=linewidth, facecolor='none', edgecolor=colors[r_cat[i]])
            ax1_ds.scatter(r_in_r200[i], delv, marker=this_mk, s=this_size,
                           edgecolor=colors[int(dsp_alloc[i]) % len(colors)],
                           alpha=this_alpha, linewidth=linewidth, facecolor='none')
            ax_ds.scatter(spec_cat['col2'][i], spec_cat['col3'][i], transform=ax_ds.get_transform('fk5'),
                          marker=this_mk, s=this_size, alpha=this_alpha, linewidth=linewidth,
                          facecolor='none', edgecolor=colors[int(dsp_alloc[i]) % len(colors)])



    ax.set_xlabel('R.A.', fontsize=25)
    ax.set_ylabel('Decl.', fontsize=25, labelpad=-1)
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.tick_params(direction='in')
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.invert_xaxis()

    ax_ds.set_xlim(x_min, x_max)
    ax_ds.set_ylim(y_min, y_max)
    ax_ds.set_xlabel('R.A.', fontsize=25)
    ax_ds.set_ylabel('Decl.', fontsize=25, labelpad=-1)
    ax_ds.tick_params(axis='both', which='major', labelsize=20)
    ax_ds.tick_params(direction='in')
    ax_ds.text(0.03, 0.95, f'{ab.clusters[k]}', transform=ax_ds.transAxes, fontsize=25)
    ax_ds.invert_xaxis()
    if k == 3:
        ax1_58_ds.set_ylim([-4, 4])
        ax1_58_ds.set_xlabel(r'R/R$_{200}$', fontsize=25)
        ax1_58_ds.set_ylabel(r'$\Delta$v$_{r}$ / $\sigma_{v_{r}}$', fontsize=25, labelpad=-1)
        ax1_58_ds.tick_params(axis='both', which='major', labelsize=20)
        ax1_58_ds.text(0.03, 0.95, 'a3558', transform=ax1_58_ds.transAxes, fontsize=25)
        ax1_62_ds.set_ylim([-4, 4])
        ax1_62_ds.set_xlabel(r'R/R$_{200}$', fontsize=25)
        ax1_62_ds.set_ylabel(r'$\Delta$v$_{r}$ / $\sigma_{v_{r}}$', fontsize=25, labelpad=-1)
        ax1_62_ds.tick_params(axis='both', which='major', labelsize=20)
        ax1_62_ds.text(0.03, 0.95, 'a3562', transform=ax1_62_ds.transAxes, fontsize=25)
    else:
        ax1_ds.set_ylim([-4, 4])
        ax1_ds.set_xlabel(r'R/R$_{200}$', fontsize=25)
        ax1_ds.set_ylabel(r'$\Delta$v$_{r}$ / $\sigma_{v_{r}}$', fontsize=25, labelpad=-1)
        ax1_ds.tick_params(axis='both', which='major', labelsize=20)
        ax1_ds.text(0.03, 0.95, f'{ab.clusters[k]}', transform=ax1_ds.transAxes, fontsize=25)

    if k == 3:
        c2_58 = SphericalCircle((coords_ss[0].ra, coords_ss[0].dec), r200_ss[0] * arcmin_Mpc.value * u.arcmin,
                             edgecolor='green', facecolor='none', linestyle='--', transform=ax.get_transform('fk5'),
                             linewidth=3, alpha=0.2)
        c2_62 = SphericalCircle((coords_ss[1].ra, coords_ss[1].dec), r200_ss[1] * arcmin_Mpc.value * u.arcmin,
                             edgecolor='brown', facecolor='none', linestyle='--', transform=ax.get_transform('fk5'),
                             linewidth=3,alpha=0.2)
        ax.add_patch(c2_58)
        ax.add_patch(c2_62)
        c2_58_ds = SphericalCircle((coords_ss[0].ra, coords_ss[0].dec), r200_ss[0] * arcmin_Mpc.value * u.arcmin,
                             edgecolor='green', facecolor='none', linestyle='--',
                                   transform=ax_ds.get_transform('fk5'), linewidth=3, alpha=0.2)
        c2_62_ds = SphericalCircle((coords_ss[1].ra, coords_ss[1].dec), r200_ss[1] * arcmin_Mpc.value * u.arcmin,
                             edgecolor='brown', facecolor='none', linestyle='--',
                                   transform=ax_ds.get_transform('fk5'), linewidth=3, alpha=0.2)
        ax_ds.add_patch(c2_58_ds)
        ax_ds.add_patch(c2_62_ds)
    elif k == 6:
        c2_58 = SphericalCircle((coords_3716[0].ra, coords_3716[0].dec), r200_3716[0] * arcmin_Mpc.value * u.arcmin,
                             edgecolor='orange', facecolor='none', linestyle='--',
                             transform=ax.get_transform('fk5'), linewidth=3, alpha=0.2)
        c2_62 = SphericalCircle((coords_3716[1].ra, coords_3716[1].dec), r200_3716[1] * arcmin_Mpc.value * u.arcmin,
                             edgecolor='green', facecolor='none', linestyle='--',
                             transform=ax.get_transform('fk5'), linewidth=3, alpha=0.2)
        ax.add_patch(c2_58)
        ax.add_patch(c2_62)
        c2_58_ds = SphericalCircle((coords_3716[0].ra, coords_3716[0].dec), r200_3716[0] * arcmin_Mpc.value * u.arcmin,
                             edgecolor='orange', facecolor='none', linestyle='--',
                             transform=ax_ds.get_transform('fk5'), linewidth=3, alpha=0.2)
        c2_62_ds = SphericalCircle((coords_3716[1].ra, coords_3716[1].dec), r200_3716[1] * arcmin_Mpc.value * u.arcmin,
                             edgecolor='green', facecolor='none', linestyle='--',
                             transform=ax_ds.get_transform('fk5'), linewidth=3, alpha=0.2)
        ax_ds.add_patch(c2_58_ds)
        ax_ds.add_patch(c2_62_ds)
    else:
        c2_ds = SphericalCircle((ra, dec) * u.deg,
                                ab.r200[k] * arcmin_Mpc.value * u.arcmin,
                                edgecolor='black',
                                facecolor='none',
                                linestyle='--',
                                transform=ax_ds.get_transform('fk5'),
                                linewidth=3,
                                alpha=0.2)
        ax_ds.add_patch(c2_ds)
        c2 = SphericalCircle((ra, dec) * u.deg,
                             ab.r200[k] * arcmin_Mpc.value * u.arcmin,
                             edgecolor='black',
                             facecolor='none',
                             linestyle='--',
                             transform=ax.get_transform('fk5'),
                             linewidth=3,
                             alpha=0.2)
        ax.add_patch(c2)
        ax.text(0.03, 0.95, f'{ab.clusters[k]}', transform=ax.transAxes, fontsize=25)

    if k == 3 or k == 6:
        frac_PM_58 = np.zeros(6)
        frac_I_58 = np.zeros(6)
        frac_PM_62 = np.zeros(6)
        frac_I_62 = np.zeros(6)
        for i in range(6):
            print(f'i58:{i}, no:{ncnt_58[0, i]}, pm:{ncnt_58[1, i]}, I:{ncnt_58[2, i]}')
            print(f'i62:{i}, no:{ncnt_62[0, i]}, pm:{ncnt_62[1, i]}, I:{ncnt_62[2, i]}')
            print(f'i58_3:{i}, no:{ncnt_58_3[0, i]}, pm:{ncnt_58_3[1, i]}, I:{ncnt_58_3[2, i]}')
            print(f'i62_3:{i}, no:{ncnt_62_3[0, i]}, pm:{ncnt_62_3[1, i]}, I:{ncnt_62_3[2, i]}')
            frac_PM_58[i] = ncnt_58[1, i] / np.sum(ncnt_58[0:3, i]) * 100
            frac_I_58[i] = ncnt_58[2, i] / np.sum(ncnt_58[0:3, i]) * 100
            frac_PM_62[i] = ncnt_62[1, i] / np.sum(ncnt_62[0:3, i]) * 100
            frac_I_62[i] = ncnt_62[2, i] / np.sum(ncnt_62[0:3, i]) * 100
        print(f'else: no:{ncnt_else[0]}, pm:{ncnt_else[1]}, I:{ncnt_else[2]}')

        ax1_in_58.plot(range(6), frac_PM_58, ls=':', c='orange', label='PM', lw=2)
        ax1_in_58.plot(range(6), frac_I_58, c='blue', label='M', lw=2)
        ax1_in_58.legend()
        ax1_in_58.set_xticks(range(6))
        ax1_in_58.set_xticklabels(['ob', 'A', 'B', 'C', 'D', 'E'], fontsize=10)
        ax1_in_58.set_ylabel('%', fontsize=10)
        ax1_in_58.tick_params(direction='in')

        ax1_in_62.plot(range(6), frac_PM_62, ls=':', c='orange', label='PM', lw=2)
        ax1_in_62.plot(range(6), frac_I_62, c='blue', label='M', lw=2)
        ax1_in_62.legend()
        ax1_in_62.set_xticks(range(6))
        ax1_in_62.set_xticklabels(['ob', 'A', 'B', 'C', 'D', 'E'], fontsize=10)
        ax1_in_62.set_ylabel('%', fontsize=10)
        ax1_in_62.tick_params(direction='in')

        ax1_58.set_xlabel(r'R/R$_{200}$', fontsize=25)
        ax1_58.set_ylabel(r'$\Delta$v$_{r}$ / $\sigma_{v_{r}}$', fontsize=25, labelpad=-1)
        ax1_58.tick_params(axis='both', which='major', labelsize=20)
        ax1_62.set_xlabel(r'R/R$_{200}$', fontsize=25)
        ax1_62.set_ylabel(r'$\Delta$v$_{r}$ / $\sigma_{v_{r}}$', fontsize=25, labelpad=-1)
        ax1_62.tick_params(axis='both', which='major', labelsize=20)
    else:
        frac_PM_this = np.zeros(6)
        frac_I_this = np.zeros(6)
        for i in range(6):
            print(f'i:{i}, no:{ncnt_this[0, i]}, pm:{ncnt_this[1, i]}, I:{ncnt_this[2, i]}')
            print(f'i_3:{i}, no:{ncnt_this_3[0, i]}, pm:{ncnt_this_3[1, i]}, I:{ncnt_this_3[2, i]}')
            frac_PM_this[i] = ncnt_this[1, i] / np.sum(ncnt_this[0:3, i]) * 100
            frac_I_this[i] = ncnt_this[2, i] / np.sum(ncnt_this[0:3, i]) * 100
        ax1_in.plot(range(6), frac_PM_this, ls=':', c='orange', label='PM', lw=2)
        ax1_in.plot(range(6), frac_I_this, c='blue', label='M', lw=2)
        ax1_in.legend()
        ax1_in.set_xticks(range(6))
        ax1_in.set_xticklabels(['ob', 'A', 'B', 'C', 'D', 'E'], fontsize=10)
        ax1_in.set_ylabel('%', fontsize=10)
        ax1_in.tick_params(direction='in')
        ax1.set_ylim([-4, 4])
        ax1.text(0.03, 0.95, f'{ab.clusters[k]}', transform=ax1.transAxes, fontsize=25)
        ax1.set_xlabel(r'R/R$_{200}$', fontsize=25)
        ax1.set_ylabel(r'$\Delta$v$_{r}$ / $\sigma_{v_{r}}$', fontsize=25, labelpad=-1)
        ax1.tick_params(axis='both', which='major', labelsize=20)

    ins_axes_I = inset_axes(ax,
                            width='100%',
                            height='100%',
                            bbox_to_anchor=(0.05, 0.02, 0.3, 0.1),
                            bbox_transform=ax.transAxes)
    ins_axes_I.set_title(r'% of    (Merging)', fontsize=18)
    ins_axes_I.tick_params(axis='x', which='major', labelsize=0)
    ins_axes_I.set_xticks([])
    ins_axes_PM = inset_axes(ax,
                             width='100%',
                             height='100%',
                             bbox_to_anchor=(0.7, 0.02, 0.3, 0.1),
                             bbox_transform=ax.transAxes)
    ins_axes_PM.set_title(r'% of    (Post-merging)', fontsize=18)
    ins_axes_PM.tick_params(axis='x', which='major', labelsize=0)
    ins_axes_PM.set_xticks([])
    ax.scatter(0.15, 0.134, marker='*', s=220, transform=ax.transAxes, facecolor='none', edgecolor='black', linewidth=3)
    ax.scatter(0.76, 0.134, marker='P', s=200, transform=ax.transAxes, facecolor='none', edgecolor='black', linewidth=3)

    ins_axes_I_ds = inset_axes(ax_ds,
                            width='100%',
                            height='100%',
                            bbox_to_anchor=(0.05, 0.02, 0.3, 0.1),
                            bbox_transform=ax_ds.transAxes)
    ins_axes_I_ds.set_title(r'% of    (Merging)', fontsize=18)
    ins_axes_I_ds.tick_params(axis='x', which='major', labelsize=0)
    ins_axes_I_ds.set_xticks([])
    ins_axes_PM_ds = inset_axes(ax_ds,
                             width='100%',
                             height='100%',
                             bbox_to_anchor=(0.7, 0.02, 0.3, 0.1),
                             bbox_transform=ax_ds.transAxes)
    ins_axes_PM_ds.set_title(r'% of    (Post-merging)', fontsize=18)
    ins_axes_PM_ds.tick_params(axis='x', which='major', labelsize=0)
    ins_axes_PM_ds.set_xticks([])
    ax_ds.scatter(0.15, 0.134, marker='*', s=220, transform=ax_ds.transAxes, facecolor='none', edgecolor='black', linewidth=3)
    ax_ds.scatter(0.76, 0.134, marker='P', s=200, transform=ax_ds.transAxes, facecolor='none', edgecolor='black', linewidth=3)

    fsub_fmain = np.zeros((4, 2))  # 0: M-type by two, 1: PM-type by two, 2: M-type by three, 3: PM-type by three

    n_t_main, n_m_main, n_pm_main, n_t_main3, n_m_main3, n_pm_main3 = 0, 0, 0, 0, 0, 0
    n_t_sub, n_m_sub, n_pm_sub, n_t_sub3, n_m_sub3, n_pm_sub3 = 0, 0, 0, 0, 0, 0
    sub_frac_I = np.zeros(this_num_sub)
    sub_frac_PM = np.zeros(this_num_sub)
    for i in range(this_num_sub):
        sub_frac_I[i] = cnt_feat_sub2[2, i] / sum(cnt_feat_sub2[:, i]) * 100
        sub_frac_PM[i] = cnt_feat_sub2[1, i] / sum(cnt_feat_sub2[:, i]) * 100
        if k == 3:
            if i == comp_no_main[1][k] or i == comp_no_main[1][k+1]:
                n_t_main += sum(cnt_feat_sub2[:, i])
                n_m_main += cnt_feat_sub2[2, i]
                n_pm_main += cnt_feat_sub2[1, i]
                n_t_main3 += sum(cnt_feat_sub3[:, i])
                n_m_main3 += cnt_feat_sub3[2, i]
                n_pm_main3 += cnt_feat_sub3[1, i]
            else:
                n_t_sub += sum(cnt_feat_sub2[:, i])
                n_m_sub += cnt_feat_sub2[2, i]
                n_pm_sub += cnt_feat_sub2[1, i]
                n_t_sub3 += sum(cnt_feat_sub3[:, i])
                n_m_sub3 += cnt_feat_sub3[2, i]
                n_pm_sub3 += cnt_feat_sub3[1, i]
        else:
            if i == comp_no_main[1][k]:
                n_t_main += sum(cnt_feat_sub2[:, i])
                n_m_main += cnt_feat_sub2[2, i]
                n_pm_main += cnt_feat_sub2[1, i]
                n_t_main3 += sum(cnt_feat_sub3[:, i])
                n_m_main3 += cnt_feat_sub3[2, i]
                n_pm_main3 += cnt_feat_sub3[1, i]
            else:
                n_t_sub += sum(cnt_feat_sub2[:, i])
                n_m_sub += cnt_feat_sub2[2, i]
                n_pm_sub += cnt_feat_sub2[1, i]
                n_t_sub3 += sum(cnt_feat_sub3[:, i])
                n_m_sub3 += cnt_feat_sub3[2, i]
                n_pm_sub3 += cnt_feat_sub3[1, i]
    fsub_fmain[0, 0] = (n_m_sub / n_t_sub - n_m_main / n_t_main) * 100
    fsub_fmain[1, 0] = (n_pm_sub / n_t_sub - n_pm_main / n_t_main) * 100
    fsub_fmain[0, 1] = (n_m_sub / n_t_sub ** 2 + n_m_main / n_t_main ** 2) ** 0.5 * 100
    fsub_fmain[1, 1] = (n_pm_sub / n_t_sub ** 2 + n_pm_main / n_t_main ** 2) ** 0.5 * 100
    fsub_fmain[2, 0] = (n_m_sub3 / n_t_sub3 - n_m_main3 / n_t_main3) * 100
    fsub_fmain[3, 0] = (n_pm_sub3 / n_t_sub3 - n_pm_main3 / n_t_main3) * 100

    barlistI = ins_axes_I.bar(range(this_num_sub), sub_frac_I, width=0.5)
    barlistPM = ins_axes_PM.bar(range(this_num_sub), sub_frac_PM, width=0.5)
    ins_axes_I.set_ylim(0, 35)
    ins_axes_PM.set_ylim(0, 35)
    for i in range(this_num_sub):
        barlistI[i].set_color(colors[i + 1])
        barlistPM[i].set_color(colors[i + 1])

    fsub_fmain_dsp = np.zeros((4, 2))  # 0: M-type by two, 1: PM-type by two, 2: M-type by three, 3: PM-type by three

    n_t_main, n_m_main, n_pm_main, n_m_main3, n_pm_main3 = 0, 0, 0, 0, 0
    n_t_sub, n_m_sub, n_pm_sub, n_m_sub3, n_pm_sub3 = 0, 0, 0, 0, 0
    sub_frac_I = []
    sub_frac_PM = []

    for i in range(len(cnt_feat_sub_dsp[0])):
        n_t_this = sum(cnt_feat_sub_dsp[:, i])
        if n_t_this > 0:
            sub_frac_I.append((cnt_feat_sub_dsp[2, i] + cnt_feat_sub_dsp[4, i]) / n_t_this * 100)
            sub_frac_PM.append((cnt_feat_sub_dsp[1, i] + cnt_feat_sub_dsp[3, i]) / n_t_this * 100)
            if i == 0:
                n_t_main += n_t_this
                n_m_main += cnt_feat_sub_dsp[2, i] + cnt_feat_sub_dsp[4, i]
                n_pm_main += cnt_feat_sub_dsp[1, i] + cnt_feat_sub_dsp[3, i]
                n_m_main3 += cnt_feat_sub_dsp[4, i]
                n_pm_main3 += cnt_feat_sub_dsp[3, i]
            else:
                n_t_sub += n_t_this
                n_m_sub += cnt_feat_sub_dsp[2, i] + cnt_feat_sub_dsp[4, i]
                n_pm_sub += cnt_feat_sub_dsp[1, i] + cnt_feat_sub_dsp[3, i]
                n_m_sub3 += cnt_feat_sub_dsp[2, i]
                n_pm_sub3 += cnt_feat_sub_dsp[1, i]

    if n_t_sub > 0:
        fsub_fmain_dsp[0, 0] = (n_m_sub / n_t_sub - n_m_main / n_t_main) * 100
        fsub_fmain_dsp[1, 0] = (n_pm_sub / n_t_sub - n_pm_main / n_t_main) * 100
        fsub_fmain_dsp[0, 1] = (n_m_sub / n_t_sub ** 2 + n_m_main / n_t_main ** 2) ** 0.5 * 100
        fsub_fmain_dsp[1, 1] = (n_pm_sub / n_t_sub ** 2 + n_pm_main / n_t_main **2) ** 0.5 * 100
        fsub_fmain_dsp[2, 0] = (n_m_sub3 / n_t_sub - n_m_main3 / n_t_main) * 100
        fsub_fmain_dsp[3, 0] = (n_pm_sub3 / n_t_sub - n_pm_main3 / n_t_main) * 100

        barlistI_ds = ins_axes_I_ds.bar(range(len(sub_frac_I)), sub_frac_I, width=0.5)
        barlistPM_ds = ins_axes_PM_ds.bar(range(len(sub_frac_PM)), sub_frac_PM, width=0.5)
        ins_axes_I_ds.set_ylim(0, 35)
        ins_axes_PM_ds.set_ylim(0, 35)
        for i in range(len(sub_frac_I)):
            barlistI_ds[i].set_color(colors[i % len(colors) - 1])  # black for background is at the last of the color array
            barlistPM_ds[i].set_color(colors[i % len(colors) - 1])

        print(f'{this_num_sub_dsp}')
        ### plot mclust and DS+ in one plot
        offset = 0.25 if fsub_fmain[0, 0] == fsub_fmain[2, 0] and fsub_fmain_dsp[0, 0] == fsub_fmain_dsp[2, 0] else 0
        axs16.errorbar(fsub_fmain[0, 0] - offset, fsub_fmain_dsp[0, 0], xerr=fsub_fmain[0, 1], yerr=fsub_fmain_dsp[0, 1],
                           c='blue', fmt=markers[k], label=ab.clusters[k]+' M', markersize=10)
        axs16.errorbar(fsub_fmain[1, 0], fsub_fmain_dsp[1, 0], xerr=fsub_fmain[1, 1], yerr=fsub_fmain_dsp[1, 1],
                           c='orange', fmt=markers[k], label=ab.clusters[k]+' PM', markersize=10)
        axs16.scatter(fsub_fmain[2, 0] + offset, fsub_fmain_dsp[2, 0], edgecolors='blue', marker=markers[k], s=100, facecolors='none')
        axs16.scatter(fsub_fmain[3, 0], fsub_fmain_dsp[3, 0], edgecolors='orange', marker=markers[k], s=100, facecolors='none')
        print(fsub_fmain)
        print(fsub_fmain_dsp)

    if k == 3:
        ax1_58.text(0.03, 0.95, 'a3558', transform=ax.transAxes, fontsize=25)
        ax1_62.text(0.03, 0.95, 'a3562', transform=ax.transAxes, fontsize=25)
    elif k == 6:
        ax1_58.text(0.03, 0.95, 'a3716S', transform=ax.transAxes, fontsize=25)
        ax1_62.text(0.03, 0.95, 'a3716N', transform=ax.transAxes, fontsize=25)

    if k == 3:
        for i in range(len(coords_ss)):
            ax.scatter(coords_ss[i].ra.value, coords_ss[i].dec.value, marker='x', s=150, c='black',
                       transform=ax.get_transform('fk5'))
            ax.text(coords_ss[i].ra.value, coords_ss[i].dec.value, cluster_ss[i], c='black',
                    transform=ax.get_transform('fk5'))
            ax_ds.scatter(coords_ss[i].ra.value, coords_ss[i].dec.value, marker='x', s=150, c='black',
                       transform=ax_ds.get_transform('fk5'))
            ax_ds.text(coords_ss[i].ra.value, coords_ss[i].dec.value, cluster_ss[i], c='black',
                    transform=ax_ds.get_transform('fk5'))

    fig_ds.savefig(ab.work_dir + 'plot/' + ab.clusters[k] + f'_rosat_substruc_ds.png')

    # draw PPS region lines
    if k == 3:
        draw_phase_region(ax1_58, 0.3)
        draw_escape_line(ax1_58, 0.3, r200_ss[0], m200_ss[0], rad_vel_sig_ss[0])
        draw_phase_region(ax1_62, 0.3)
        draw_escape_line(ax1_62, 0.3, r200_ss[1], m200_ss[1], rad_vel_sig_ss[1])
        fig1_58.savefig(ab.work_dir + 'plot/' + ab.clusters[k] + f'_phase_58.png')
        fig1_62.savefig(ab.work_dir + 'plot/' + ab.clusters[k] + f'_phase_62.png')
        fig1_58_ds.savefig(ab.work_dir + 'plot/' + ab.clusters[k] + f'_phase_ds_58.png')
        fig1_62_ds.savefig(ab.work_dir + 'plot/' + ab.clusters[k] + f'_phase_ds_62.png')
        fig.savefig(ab.work_dir + 'plot/' + ab.clusters[k] + f'_rosat_substruc.png')
    elif k == 6:
        draw_phase_region(ax1_58, 0.3)
        draw_escape_line(ax1_58, 0.3, r200_3716[0], m200_3716[0], rad_vel_sig_3716[0])
        draw_phase_region(ax1_62, 0.3)
        draw_escape_line(ax1_62, 0.3, r200_3716[1], m200_3716[1], rad_vel_sig_3716[1])
        fig1_58.savefig(ab.work_dir + 'plot/' + ab.clusters[k] + f'_phase_S.png')
        fig1_62.savefig(ab.work_dir + 'plot/' + ab.clusters[k] + f'_phase_N.png')
        fig1_ds.savefig(ab.work_dir + 'plot/' + ab.clusters[k] + f'_phase_ds.png')
        fig.savefig(ab.work_dir + 'plot/' + ab.clusters[k] + f'_rosat_substruc.png')
    else:
        draw_phase_region(ax1, 0.3)
        draw_escape_line(ax1, 0.3, ab.r200[k], ab.m200[k], ab.sig_vr[k])
        fig1.savefig(ab.work_dir + 'plot/' + ab.clusters[k] + f'_phase.png')
        fig.savefig(ab.work_dir + 'plot/' + ab.clusters[k] + f'_rosat_substruc.png')
        fig1_ds.savefig(ab.work_dir + 'plot/' + ab.clusters[k] + f'_phase_ds.png')



pm_frac_reg = np.zeros(6)
i_frac_reg = np.zeros(6)
pmerr_frac_reg = np.zeros(6)
ierr_frac_reg = np.zeros(6)
pm_frac_reg_3 = np.zeros(6)
i_frac_reg_3 = np.zeros(6)
pmerr_frac_reg_3 = np.zeros(6)
ierr_frac_reg_3 = np.zeros(6)
for i in range(6):
    total_num = ncnt_reg[0, i] + ncnt_reg[1, i] + ncnt_reg[2, i]
    total_num_3 = ncnt_reg_3[0, i] + ncnt_reg_3[1, i] + ncnt_reg_3[2, i]
    print(f'{total_num} {total_num_3}')
    pm_frac_reg[i] = ncnt_reg[1, i] / total_num * 100
    pm_frac_reg_3[i] = ncnt_reg_3[1, i] / total_num * 100
    i_frac_reg[i] = ncnt_reg[2, i] / total_num * 100
    i_frac_reg_3[i] = ncnt_reg_3[2, i] / total_num * 100
    pmerr_frac_reg[i] = np.sqrt(ncnt_reg[1, i]) / total_num * 100
    pmerr_frac_reg_3[i] = np.sqrt(ncnt_reg_3[1, i]) / total_num * 100
    ierr_frac_reg[i] = np.sqrt(ncnt_reg[2, i]) / total_num * 100
    ierr_frac_reg_3[i] = np.sqrt(ncnt_reg_3[2, i]) / total_num * 100
axs133.plot(range(6), pm_frac_reg, c='orange', label='Post-merging (2)', lw=5)
axs133.plot(range(6), i_frac_reg, c='blue', label='Merging (2)', lw=5)
axs133.plot(range(6), pm_frac_reg_3, c='orange', ls=':', label='Post-merging (3)', lw=5)
axs133.plot(range(6), i_frac_reg_3, c='blue', ls=':', label='Merging (3)', lw=5)
axs133.fill_between(range(6),
                           [a - b for a, b in zip(pm_frac_reg, pmerr_frac_reg)],
                           [a + b for a, b in zip(pm_frac_reg, pmerr_frac_reg)],
                           alpha=.15, facecolor='orange')
axs133.fill_between(range(6),
                           [a - b for a, b in zip(i_frac_reg, ierr_frac_reg)],
                           [a + b for a, b in zip(i_frac_reg, ierr_frac_reg)],
                           alpha=.15, facecolor='blue')

axs133.set_xticks(range(6))
axs133.set_xticklabels(['ob', 'A', 'B', 'C', 'D', 'E'],
                      fontsize=20)
axs133.set_ylabel('%', fontsize=20)
axs133.tick_params(axis='both', which='major', labelsize=20)
axs133.tick_params(axis='both', which='both', right=True, labelright=False, top=True, labeltop=False)
axs133.tick_params(direction='in')
axs133.legend(fontsize=15)
fig133.savefig(ab.plot_dir + f'frac_vs_reg.png')


axs122.set_xlim(0, 3.0)
axs122.set_ylim(0, 3.0)

axs122.plot([1.375, 1.5], [0.5, 0], ls='--', c='black')
axs122.plot([0.375, 1.375], [0.5, 0.5], ls='--', c='black')
axs122.plot([0.875, 1.125], [1.5, 0.5], ls='--', c='black')
axs122.plot([0.125, 0.875], [1.5, 1.5], ls='--', c='black')
axs122.plot([0.375, 0.625], [2.5, 1.5], ls='--', c='black')
axs122.plot([0.375, 3.0], [2.5, 1.0], ls='--', c='black')
axs122.plot([0, 0.375], [2.5, 2.5], ls='--', c='black')
axs122.plot([0, 0.5], [2.0, 0.0], ls='--', c='black')

axs122.tick_params(axis='both', which='major', labelsize=20)
axs122.tick_params(axis='both', which='both', right=True, labelright=False, top=True, labeltop=False, direction='in')
axs122.tick_params(direction='in')
axs122.set_xlabel(r'r$_{proj}$/R$_{200}$', fontsize=20)
axs122.set_ylabel(r'V$_{LOS}$ / $\sigma_{LOS}$', fontsize=20, labelpad=10)
axs122.text(1.5, 2.5, 'out of bound', fontsize=20)
axs122.text(1.8, 2.35, f'{ncnt_reg[0, 0]}/{ncnt_reg[1, 0]}/{ncnt_reg[2, 0]}', fontsize=12)
axs122.text(1.8, 0.8, 'A', fontsize=30)
axs122.text(1.7, 0.65, f'{ncnt_reg[0, 1]}/{ncnt_reg[1, 1]}/{ncnt_reg[2, 1]}', fontsize=12)
axs122.text(0.2, 1.9, 'B', fontsize=30)
axs122.text(0.12, 1.75, f'{ncnt_reg[0, 2]}/{ncnt_reg[1, 2]}/{ncnt_reg[2, 2]}', fontsize=12)
axs122.text(0.5, 0.9, 'C', fontsize=30)
axs122.text(0.4, 0.75, f'{ncnt_reg[0, 3]}/{ncnt_reg[1, 3]}/{ncnt_reg[2, 3]}', fontsize=12)
axs122.text(0.8, 0.2, 'D', fontsize=30)
axs122.text(0.7, 0.05, f'{ncnt_reg[0, 4]}/{ncnt_reg[1, 4]}/{ncnt_reg[2, 4]}', fontsize=12)
axs122.text(0.1, 0.3, 'E', fontsize=30)
axs122.text(0.0, 0.1, f'{ncnt_reg[0, 5]}/{ncnt_reg[1, 5]}/{ncnt_reg[2, 5]}', fontsize=12)
axs122.text(0.1, 2.8, 'N(no feature) / N(Post-merger) / N(Interacting)', fontsize=12)


fig122.savefig(ab.plot_dir + f'frac_I_vs_r_infall.png')

## for feature types selected by 2 out of 4 people
axs16.set_title(r'%$_{substructures}$ - %$_{main}$', fontsize=20)
axs16.text(0.65, 0.95, 'more in sub', fontsize=20, transform=axs16.transAxes)
# axs16.text(21.5, 16.5, 'sub', fontsize=20)
axs16.text(0.01, 0.01, 'more in main', fontsize=20, transform=axs16.transAxes)
# axs16.text(-10, -6.5, 'in main', fontsize=20)
#                       fontsize=10, rotation=60)
axs16.tick_params(axis='both', which='major', labelsize=15)
axs16.tick_params(direction='in')
axs16.axhline(0, linestyle='--', color='grey')
axs16.axvline(0, linestyle='--', color='grey')
axs16.set_ylabel('DS+', fontsize=20)
axs16.set_xlabel('mclust', fontsize=20)
axs16.legend(fontsize=10, loc='upper left')
axs16.set_xlim(-20, 20)
axs16.set_ylim(-20, 20)
fig16.savefig(ab.plot_dir + f'frac_vs_sub_all_total_frac.png')


plt.close('all')