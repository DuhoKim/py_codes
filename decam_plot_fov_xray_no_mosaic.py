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


is_dsp = True
is_dsp_first = False

plim_p = 0.1

fov_ra = [ 1.5, 2.3, 2.3, 2.7, 0.0, 0.0, 2.2]
fov_dec = [ 1.5, 2.3, 2.3, 2.3, 0.0, 0.0, 1.4]

scale_x = [180, 219, 224, 170, 0, 0, 176]
scale_y = [194, 229, 235, 215, 0, 0, 197]

colors = ['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'cyan', 'magenta', 'yellow', 'salmon','skyblue',
          'lightgrey', 'darkgrey', 'black']
markers = ['o', 's', 'p', 'P', 'h', '+', 'x', 'D']

fig2 = plt.figure(figsize=(5, 5))
ax2 = fig2.add_subplot()

fig10, axs10 = plt.subplots(tight_layout=True, figsize=(6, 6))
fig12, axs12 = plt.subplots(tight_layout=True, figsize=(12, 12))
fig122, axs122 = plt.subplots(tight_layout=True, figsize=(6, 6))
fig133, axs133 = plt.subplots(tight_layout=True, figsize=(6, 6))
fig16, axs16 = plt.subplots(tight_layout=True, figsize=(6, 6))
fig166, axs166 = plt.subplots(tight_layout=True, figsize=(6, 6))

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

coords_3716 = [SkyCoord('20:52:00 -52:45:18', unit=(u.hourangle, u.deg)),
               SkyCoord('20:51:57 -52:37:40', unit=(u.hourangle, u.deg))]


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

images = []
is_I = []
img_locs = []
img_sizes = []

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

gs = ['2', '4', '5', '6', '', '', '2']

# for k in range(0, len(clusters)):
for k in [0, 1, 2, 3, 6]:
# for k in [0]:
# for k in [0, 1, 3]:
    mm.print_time()
    print(f'{ab.clusters[k]}')

    # to print current clusters' number
    if k == 3 or k == 6:
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

    fig = plt.figure(figsize=(10, 10))
    fig_ds = plt.figure(figsize=(10, 10))
    fig1_ds = plt.figure(figsize=(10, 10))
    ax1_ds = fig1_ds.add_subplot()
    if k == 3 or k == 6:
        fig1_58 = plt.figure(figsize=(10, 10))
        fig1_62 = plt.figure(figsize=(10, 10))
        ax1_58 = fig1_58.add_subplot()
        ax1_62 = fig1_62.add_subplot()
        ax1_in_58 = ax1_58.inset_axes([0.7, 0.05, 0.25, 0.25])
        ax1_in_62 = ax1_62.inset_axes([0.7, 0.05, 0.25, 0.25])
    else:
        fig1 = plt.figure(figsize=(10, 10))
        ax1 = fig1.add_subplot()
        ax1_in = ax1.inset_axes([0.7, 0.05, 0.25, 0.25])

    vis_cat = ascii.read(
        ab.work_dir + f'sex/cat/DECam_merged_SEx_cat_{ab.clusters[k]}_Gal_ext_corrected_20rmag_psf.txt')
    spec_cat = ascii.read(ab.work_dir + f'spec/{ab.clusters[k]}_spec_ra_dec_z_zran_0.05_rran_1.5.txt')
    r_cat = ascii.read(ab.work_dir + f'spec/{ab.clusters[k]}_G{gs[k]}_class_uncer_zran_0.05_rran_1.5.csv', format='csv')
    res_cat = ascii.read(ab.work_dir + f'vis/combined/{ab.clusters[k]}.txt')

    this_num_sub = max(r_cat['classification'])

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
        ax_ds.text(0.75, 0.95, f'1 Mpc at z={z:5.3f}', fontsize=15, transform=ax.transAxes)

    ### draw high-resolution x-ray image at the center first
    with fits.open(ab.work_dir + 'fits/' + xray_high_fns[k]) as hdu_x:
        wcs_x, shape_x = find_optimal_celestial_wcs(hdu_x)  # has only CompImageHDU files
        data = hdu_x[0].data
        if xray_high_sigmas[k]:
            data = gaussian_filter(data, xray_high_sigmas[k])
        ax.contour(data, levels=xray_high_levels[k], transform=ax.get_transform(wcs_x), alpha=0.4, linewidths=3, cmap='seismic')
        ax_ds.contour(data, levels=xray_high_levels[k], transform=ax_ds.get_transform(wcs_x), alpha=0.4, linewidths=3, cmap='seismic')


    if is_dsp:
        #### DS+ test
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
                                                                          Plim_P=plim_p)
            np.savez(f'{ab.clusters[k]}_dsp_result_{plim_p}.npz', data_DSp=data_DSp, data_grs_alloc=data_grs_alloc,
                     summary_DSp_grs = summary_DSp_grs)
        else:
            dsp_arrays = np.load(f'{ab.clusters[k]}_dsp_result_{plim_p}.npz')
            data_DSp = dsp_arrays['data_DSp']
            data_grs_alloc = dsp_arrays['data_grs_alloc']
            summary_DSp_grs = dsp_arrays['summary_DSp_grs']

    if k == 3:
        r200_deg_58 = r200_ss[0] * 1e3 / kpc_arcmin / 6e1
        r200_deg_62 = r200_ss[1] * 1e3 / kpc_arcmin / 6e1
        r_in_r200_58 = coords_spec.separation(coords_ss[0]).value / r200_deg_58.value
        r_in_r200_62 = coords_spec.separation(coords_ss[1]).value / r200_deg_62.value
    elif k == 6:
        r200_deg_58 = r200_3716[0] * 1e3 / kpc_arcmin / 6e1
        r200_deg_62 = r200_3716[1] * 1e3 / kpc_arcmin / 6e1
        r_in_r200_58 = coords_spec.separation(coords_3716[0]).value / r200_deg_58.value
        r_in_r200_62 = coords_spec.separation(coords_3716[1]).value / r200_deg_62.value
    else:
        r_in_r200 = coords_spec.separation(ab.coords_cl_cen[k]).value / r200_deg.value

    ### for record subgroups by mclust
    num_0 = np.zeros(this_num_sub)      # number of no feat. in each subs
    num_1 = np.zeros(this_num_sub)      # number of PM feat. in each subs
    num_2 = np.zeros(this_num_sub)      # number of I feat. in each subs

    # for feature type selected by 3 people
    num_0_3 = np.zeros(this_num_sub)      # number of no feat. in each subs
    num_1_3 = np.zeros(this_num_sub)      # number of PM feat. in each subs
    num_2_3 = np.zeros(this_num_sub)      # number of I feat. in each subs

    ### for record subgroups by DS+
    this_num_sub_dsp = int(max(data_grs_alloc[:,8])) + 2
    num_0_dsp = np.zeros(this_num_sub_dsp)      # number of no feat. in each subs
    num_1_dsp = np.zeros(this_num_sub_dsp)      # number of PM feat. in each subs
    num_2_dsp = np.zeros(this_num_sub_dsp)      # number of I feat. in each subs

    # for feature type selected by 3 people
    num_0_dsp_3 = np.zeros(this_num_sub_dsp)      # number of no feat. in each subs
    num_1_dsp_3 = np.zeros(this_num_sub_dsp)      # number of PM feat. in each subs
    num_2_dsp_3 = np.zeros(this_num_sub_dsp)      # number of I feat. in each subs

    for i in range(len(r_cat)):
        this_mk = '.'
        this_size = 120
        this_alpha = 0.4
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
                    num_0[r_cat[i][0]-1] += 1
                    num_0_3[r_cat[i][0] - 1] += 1
                    num_0_dsp[int(data_grs_alloc[i, 8]) + 1] += 1
                    num_0_dsp_3[int(data_grs_alloc[i, 8]) + 1] += 1
                    if k == 3:      # subgroup 2 and 5 for a3558
                        if (r_cat[i][0] == 2):
                            ncnt_reg[0, phase_region(r_in_r200_58[i], delv58)] += 1
                            ncnt_reg_3[0, phase_region(r_in_r200_58[i], delv58)] += 1
                            ncnt_58[0, phase_region(r_in_r200_58[i], delv58)] += 1
                            ncnt_58_3[0, phase_region(r_in_r200_58[i], delv58)] += 1
                        elif (r_cat[i][0] == 5):       # subgroup 1 and 9 for a3562
                            ncnt_reg[0, phase_region(r_in_r200_62[i], delv62)] += 1
                            ncnt_reg_3[0, phase_region(r_in_r200_62[i], delv62)] += 1
                            ncnt_62[0, phase_region(r_in_r200_62[i], delv62)] += 1
                            ncnt_62_3[0, phase_region(r_in_r200_62[i], delv62)] += 1
                        else:
                            ncnt_else[0] += 1
                            ncnt_else_3[0] += 1
                    elif k == 6:      # subgroup 1 for A3716S
                        if (r_cat[i][0] == 1):
                            ncnt_reg[0, phase_region(r_in_r200_58[i], delv58)] += 1
                            ncnt_reg_3[0, phase_region(r_in_r200_58[i], delv58)] += 1
                            ncnt_58[0, phase_region(r_in_r200_58[i], delv58)] += 1
                            ncnt_58_3[0, phase_region(r_in_r200_58[i], delv58)] += 1
                        elif (r_cat[i][0] == 2):       # subgroup 2 for A3716N
                            ncnt_reg[0, phase_region(r_in_r200_62[i], delv62)] += 1
                            ncnt_reg_3[0, phase_region(r_in_r200_62[i], delv62)] += 1
                            ncnt_62[0, phase_region(r_in_r200_62[i], delv62)] += 1
                            ncnt_62_3[0, phase_region(r_in_r200_62[i], delv62)] += 1
                        else:
                            ncnt_else[0] += 1
                            ncnt_else_3[0] += 1
                    else:
                        ncnt_reg[0, phase_region(r_in_r200[i], delv)] += 1
                        ncnt_reg_3[0, phase_region(r_in_r200[i], delv)] += 1
                        ncnt_this[0, phase_region(r_in_r200[i], delv)] += 1
                        ncnt_this_3[0, phase_region(r_in_r200[i], delv)] += 1
                elif res_cat['col2'][ind] == 1:  # classified as PM
                    this_mk = 'P'
                    this_size = 300
                    this_alpha = 0.6
                    num_1[r_cat[i][0] - 1] += 1
                    num_1_dsp[int(data_grs_alloc[i, 8]) + 1] += 1
                    if res_cat['col3'][ind] > 2:
                        num_1_3[r_cat[i][0] - 1] += 1
                        num_1_dsp_3[int(data_grs_alloc[i, 8]) + 1] += 1
                    else:
                        num_0_3[r_cat[i][0] - 1] += 1
                        num_0_dsp_3[int(data_grs_alloc[i, 8]) + 1] += 1
                    img = Image.open(ab.work_dir + f'pics/{ab.clusters[k]}_for_all/'
                                                     f'{res_cat["col1"][ind][0]}_crop_PM.png')
                    if k == 3:
                        if (r_cat[i][0] == 2):
                            images.append(img)
                            is_I.append(False)  # 0 for PM
                            img_locs.append((r_in_r200_58[i], delv58))
                            ncnt_reg[1, phase_region(r_in_r200_58[i], delv58)] += 1
                            ncnt_58[1, phase_region(r_in_r200_58[i], delv58)] += 1
                            if res_cat['col3'][ind] > 2:
                                ncnt_reg_3[1, phase_region(r_in_r200_58[i], delv58)] += 1
                                ncnt_58_3[1, phase_region(r_in_r200_58[i], delv58)] += 1
                            else:
                                ncnt_reg_3[0, phase_region(r_in_r200_58[i], delv58)] += 1
                                ncnt_58_3[0, phase_region(r_in_r200_58[i], delv58)] += 1
                            img_sizes.append(img.size)
                            print(f"PM, sub:{r_cat[i][0]}, R/R200:{r_in_r200_58[i]}, id:{res_cat['col1'][ind][0]}")
                        elif (r_cat[i][0] == 5):       # subgroup 1 and 9 for a3562
                            images.append(img)
                            is_I.append(False)  # 0 for PM
                            img_locs.append((r_in_r200_62[i], delv62))
                            ncnt_reg[1, phase_region(r_in_r200_62[i], delv62)] += 1
                            ncnt_62[1, phase_region(r_in_r200_62[i], delv62)] += 1
                            if res_cat['col3'][ind] > 2:
                                ncnt_reg_3[1, phase_region(r_in_r200_62[i], delv62)] += 1
                                ncnt_62_3[1, phase_region(r_in_r200_62[i], delv62)] += 1
                            else:
                                ncnt_reg_3[0, phase_region(r_in_r200_62[i], delv62)] += 1
                                ncnt_62_3[0, phase_region(r_in_r200_62[i], delv62)] += 1
                            img_sizes.append(img.size)
                            print(f"PM, sub:{r_cat[i][0]}, R/R200:{r_in_r200_62[i]}, id:{res_cat['col1'][ind][0]}")
                        else:
                            ncnt_else[1] += 1
                    elif k == 6:
                        if (r_cat[i][0] == 1):
                            images.append(img)
                            is_I.append(False)  # 0 for PM
                            img_locs.append((r_in_r200_58[i], delv58))
                            ncnt_reg[1, phase_region(r_in_r200_58[i], delv58)] += 1
                            ncnt_58[1, phase_region(r_in_r200_58[i], delv58)] += 1
                            if res_cat['col3'][ind] > 2:
                                ncnt_reg_3[1, phase_region(r_in_r200_58[i], delv58)] += 1
                                ncnt_58_3[1, phase_region(r_in_r200_58[i], delv58)] += 1
                            else:
                                ncnt_reg_3[0, phase_region(r_in_r200_58[i], delv58)] += 1
                                ncnt_58_3[0, phase_region(r_in_r200_58[i], delv58)] += 1
                            img_sizes.append(img.size)
                            print(f"PM, sub:{r_cat[i][0]}, R/R200:{r_in_r200_58[i]}, id:{res_cat['col1'][ind][0]}")
                        elif (r_cat[i][0] == 2):       # subgroup 2 for A3716N
                            images.append(img)
                            is_I.append(False)  # 0 for PM
                            img_locs.append((r_in_r200_62[i], delv62))
                            ncnt_reg[1, phase_region(r_in_r200_62[i], delv62)] += 1
                            ncnt_62[1, phase_region(r_in_r200_62[i], delv62)] += 1
                            if res_cat['col3'][ind] > 2:
                                ncnt_reg_3[1, phase_region(r_in_r200_62[i], delv62)] += 1
                                ncnt_62_3[1, phase_region(r_in_r200_62[i], delv62)] += 1
                            else:
                                ncnt_reg_3[0, phase_region(r_in_r200_62[i], delv62)] += 1
                                ncnt_62_3[0, phase_region(r_in_r200_62[i], delv62)] += 1
                            img_sizes.append(img.size)
                            print(f"PM, sub:{r_cat[i][0]}, R/R200:{r_in_r200_62[i]}, id:{res_cat['col1'][ind][0]}")
                        else:
                            ncnt_else[1] += 1
                            if res_cat['col3'][ind] > 2:
                                ncnt_else_3[1] += 1
                            else:
                                ncnt_else_3[0] += 1
                    else:
                        images.append(img)
                        is_I.append(False)  # 0 for PM
                        img_locs.append((r_in_r200[i], delv))
                        ncnt_reg[1, phase_region(r_in_r200[i], delv)] += 1
                        ncnt_this[1, phase_region(r_in_r200[i], delv)] += 1
                        if res_cat['col3'][ind] > 2:
                            ncnt_reg_3[1, phase_region(r_in_r200[i], delv)] += 1
                            ncnt_this_3[1, phase_region(r_in_r200[i], delv)] += 1
                        else:
                            ncnt_reg_3[0, phase_region(r_in_r200[i], delv)] += 1
                            ncnt_this_3[0, phase_region(r_in_r200[i], delv)] += 1
                        img_sizes.append(img.size)
                        print(f"PM, sub:{r_cat[i][0]}, R/R200:{r_in_r200[i]}, id:{res_cat['col1'][ind][0]}")
                elif res_cat['col2'][ind] == 2:  # classified as I
                    this_mk = '*'
                    this_size = 300
                    this_alpha = 0.9
                    num_2[r_cat[i][0]-1] += 1
                    num_2_dsp[int(data_grs_alloc[i, 8]) + 1] += 1
                    if res_cat['col3'][ind] > 2:
                        num_2_3[r_cat[i][0] - 1] += 1
                        num_2_dsp_3[int(data_grs_alloc[i, 8]) + 1] += 1
                    else:
                        num_0_3[r_cat[i][0] - 1] += 1
                        num_0_dsp_3[int(data_grs_alloc[i, 8]) + 1] += 1
                    img = Image.open(ab.work_dir + f'pics/{ab.clusters[k]}_for_all/'
                                                   f'{res_cat["col1"][ind][0]}_crop.png')
                    if k == 3:
                        if (r_cat[i][0] == 2):
                            images.append(img)
                            is_I.append(True)
                            img_locs.append((r_in_r200_58[i], delv58))
                            ncnt_reg[2, phase_region(r_in_r200_58[i], delv58)] += 1
                            ncnt_58[2, phase_region(r_in_r200_58[i], delv58)] += 1
                            if res_cat['col3'][ind] > 2:
                                ncnt_reg_3[2, phase_region(r_in_r200_58[i], delv58)] += 1
                                ncnt_58_3[2, phase_region(r_in_r200_58[i], delv58)] += 1
                            else:
                                ncnt_reg_3[0, phase_region(r_in_r200_58[i], delv58)] += 1
                                ncnt_58_3[0, phase_region(r_in_r200_58[i], delv58)] += 1
                            img_sizes.append(img.size)
                            print(f"sub:{r_cat[i][0]}, R/R200:{r_in_r200_58[i]}, id:{res_cat['col1'][ind][0]}")
                        elif (r_cat[i][0] == 5):       # subgroup 1 and 9 for a3562
                            images.append(img)
                            is_I.append(True)
                            img_locs.append((r_in_r200_62[i], delv62))
                            ncnt_reg[2, phase_region(r_in_r200_62[i], delv62)] += 1
                            ncnt_62[2, phase_region(r_in_r200_62[i], delv62)] += 1
                            if res_cat['col3'][ind] > 2:
                                ncnt_reg_3[2, phase_region(r_in_r200_62[i], delv62)] += 1
                                ncnt_62_3[2, phase_region(r_in_r200_62[i], delv62)] += 1
                            else:
                                ncnt_reg_3[0, phase_region(r_in_r200_62[i], delv62)] += 1
                                ncnt_62_3[0, phase_region(r_in_r200_62[i], delv62)] += 1
                            img_sizes.append(img.size)
                            print(f"sub:{r_cat[i][0]}, R/R200:{r_in_r200_62[i]}, id:{res_cat['col1'][ind][0]}")
                        else:
                            ncnt_else[2] += 1
                            if res_cat['col3'][ind] > 2:
                                ncnt_else_3[2] += 1
                            else:
                                ncnt_else_3[0] += 1
                    elif k == 6:
                        if (r_cat[i][0] == 1):
                            images.append(img)
                            is_I.append(True)
                            img_locs.append((r_in_r200_58[i], delv58))
                            ncnt_reg[2, phase_region(r_in_r200_58[i], delv58)] += 1
                            ncnt_58[2, phase_region(r_in_r200_58[i], delv58)] += 1
                            if res_cat['col3'][ind] > 2:
                                ncnt_reg_3[2, phase_region(r_in_r200_58[i], delv58)] += 1
                                ncnt_58_3[2, phase_region(r_in_r200_58[i], delv58)] += 1
                            else:
                                ncnt_reg_3[0, phase_region(r_in_r200_58[i], delv58)] += 1
                                ncnt_58_3[0, phase_region(r_in_r200_58[i], delv58)] += 1
                            img_sizes.append(img.size)
                            print(f"sub:{r_cat[i][0]}, R/R200:{r_in_r200_58[i]}, id:{res_cat['col1'][ind][0]}")
                        elif (r_cat[i][0] == 2):       # subgroup 2 for a3716N
                            images.append(img)
                            is_I.append(True)
                            img_locs.append((r_in_r200_62[i], delv62))
                            ncnt_reg[2, phase_region(r_in_r200_62[i], delv62)] += 1
                            ncnt_62[2, phase_region(r_in_r200_62[i], delv62)] += 1
                            if res_cat['col3'][ind] > 2:
                                ncnt_reg_3[2, phase_region(r_in_r200_62[i], delv62)] += 1
                                ncnt_62_3[2, phase_region(r_in_r200_62[i], delv62)] += 1
                            else:
                                ncnt_reg_3[0, phase_region(r_in_r200_62[i], delv62)] += 1
                                ncnt_62_3[0, phase_region(r_in_r200_62[i], delv62)] += 1
                            img_sizes.append(img.size)
                            print(f"sub:{r_cat[i][0]}, R/R200:{r_in_r200_62[i]}, id:{res_cat['col1'][ind][0]}")
                        else:
                            ncnt_else[2] += 1
                    else:
                        images.append(img)
                        is_I.append(True)
                        img_locs.append((r_in_r200[i], delv))
                        ncnt_reg[2, phase_region(r_in_r200[i], delv)] += 1
                        ncnt_this[2, phase_region(r_in_r200[i], delv)] += 1
                        if res_cat['col3'][ind] > 2:
                            ncnt_reg_3[2, phase_region(r_in_r200[i], delv)] += 1
                            ncnt_this_3[2, phase_region(r_in_r200[i], delv)] += 1
                        else:
                            ncnt_reg_3[0, phase_region(r_in_r200[i], delv)] += 1
                            ncnt_this_3[0, phase_region(r_in_r200[i], delv)] += 1
                        img_sizes.append(img.size)
                        print(f'sub:{r_cat[i][0]}, R/R200:{r_in_r200[i]}, id:{res_cat["col1"][ind][0]}')

                if k == 3:
                    if (r_cat[i][0] == 2):
                        ax1_58.scatter(r_in_r200_58[i], delv58,
                                    marker=this_mk,
                                    s=this_size,
                                    alpha=this_alpha,
                                    linewidth=linewidth,
                                    facecolor='none',
                                    edgecolor=colors[r_cat[i][0]])
                    elif (r_cat[i][0] == 5):       # subgroup 1 and 9 for a3562
                        ax1_62.scatter(r_in_r200_62[i], delv62,
                                       marker=this_mk,
                                       s=this_size,
                                       alpha=this_alpha,
                                       linewidth=linewidth,
                                       facecolor='none',
                                       edgecolor=colors[r_cat[i][0]])
                elif k == 6:
                    if (r_cat[i][0] == 1):
                        ax1_58.scatter(r_in_r200_58[i], delv58,
                                    marker=this_mk,
                                    s=this_size,
                                    alpha=this_alpha,
                                    linewidth=linewidth,
                                    facecolor='none',
                                    edgecolor=colors[r_cat[i][0]])
                    elif (r_cat[i][0] == 2):
                        ax1_62.scatter(r_in_r200_62[i], delv62,
                                       marker=this_mk,
                                       s=this_size,
                                       alpha=this_alpha,
                                       linewidth=linewidth,
                                       facecolor='none',
                                       edgecolor=colors[r_cat[i][0]])

        if k == 3:
            linewidth = 2
            this_size -= 30
        elif k == 6:
            linewidth = 2
            this_size -= 30
        else:
            linewidth = 3
            ax1.scatter(sep / r200_deg, delv,
                        marker=this_mk,
                        s=this_size,
                        alpha=this_alpha,
                        linewidth=linewidth,
                        facecolor='none',
                        edgecolor=colors[r_cat[i][0]])
        if is_dsp:
            ax1_ds.scatter(sep / r200_deg, delv,
                           marker=this_mk,
                           s=this_size,
                           edgecolor=colors[int(data_grs_alloc[i, 8]) % len(colors)],
                           alpha=this_alpha,
                           linewidth=linewidth,
                           facecolor='none'
                           )
        ax.scatter(spec_cat['col2'][i], spec_cat['col3'][i], transform=ax.get_transform('fk5'),
                   marker=this_mk,
                   s=this_size,
                   alpha=this_alpha,
                   linewidth=linewidth,
                   facecolor='none',
                   edgecolor=colors[r_cat[i][0]])

        ax_ds.scatter(spec_cat['col2'][i], spec_cat['col3'][i], transform=ax.get_transform('fk5'),
                       marker=this_mk,
                       s=this_size,
                       # c=colors[r_cat[i][0]],
                       alpha=this_alpha,
                       linewidth=linewidth,
                       facecolor='none',
                       edgecolor=colors[int(data_grs_alloc[i, 8]) % len(colors)])


    ax.set_xlabel('R.A.', fontsize=25)
    ax.set_ylabel('Decl.', fontsize=25, labelpad=-1)
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.tick_params(direction='in')

    if is_dsp:
        ax_ds.set_xlabel('R.A.', fontsize=25)
        ax_ds.set_ylabel('Decl.', fontsize=25, labelpad=-1)
        ax_ds.tick_params(axis='both', which='major', labelsize=20)
        ax_ds.tick_params(direction='in')

        ax1_ds.set_xlabel(r'R/R$_{200}$', fontsize=25)
        ax1_ds.set_ylabel(r'$\Delta$v$_{r}$ / $\sigma_{v_{r}}$', fontsize=25, labelpad=-1)
        ax1_ds.tick_params(axis='both', which='major', labelsize=20)
        ax_ds.set_xlim(x_min, x_max)
        ax_ds.set_ylim(y_min, y_max)


    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)

    if k == 3:
        c2_58 = SphericalCircle((coords_ss[0].ra, coords_ss[0].dec),
                             r200_ss[0] * arcmin_Mpc.value * u.arcmin,
                             edgecolor='green',
                             facecolor='none',
                             linestyle='--',
                             transform=ax.get_transform('fk5'),
                             linewidth=3,
                             alpha=0.2)

        c2_62 = SphericalCircle((coords_ss[1].ra, coords_ss[1].dec) * u.deg,
                             r200_ss[1] * arcmin_Mpc.value * u.arcmin,
                             edgecolor='brown',
                             facecolor='none',
                             linestyle='--',
                             transform=ax.get_transform('fk5'),
                             linewidth=3,
                             alpha=0.2)
        if is_dsp:
            c2_58_ds = SphericalCircle((coords_ss[0].ra, coords_ss[0].dec),
                                 r200_ss[0] * arcmin_Mpc.value * u.arcmin,
                                 edgecolor='green',
                                 facecolor='none',
                                 linestyle='--',
                                 transform=ax.get_transform('fk5'),
                                 linewidth=3,
                                 alpha=0.2)

            c2_62_ds = SphericalCircle((coords_ss[1].ra, coords_ss[1].dec) * u.deg,
                                 r200_ss[1] * arcmin_Mpc.value * u.arcmin,
                                 edgecolor='brown',
                                 facecolor='none',
                                 linestyle='--',
                                 transform=ax.get_transform('fk5'),
                                 linewidth=3,
                                 alpha=0.2)
            ax_ds.add_patch(c2_58_ds)
            ax_ds.add_patch(c2_62_ds)
        ax.add_patch(c2_58)
        ax.add_patch(c2_62)
    elif k == 6:
        c2_58 = SphericalCircle((coords_3716[0].ra, coords_3716[0].dec),
                             r200_3716[0] * arcmin_Mpc.value * u.arcmin,
                             edgecolor='orange',
                             facecolor='none',
                             linestyle='--',
                             transform=ax.get_transform('fk5'),
                             linewidth=3,
                             alpha=0.2)

        c2_62 = SphericalCircle((coords_3716[1].ra, coords_3716[1].dec) * u.deg,
                             r200_3716[1] * arcmin_Mpc.value * u.arcmin,
                             edgecolor='green',
                             facecolor='none',
                             linestyle='--',
                             transform=ax.get_transform('fk5'),
                             linewidth=3,
                             alpha=0.2)
        if is_dsp:
            c2_58_ds = SphericalCircle((coords_3716[0].ra, coords_3716[0].dec),
                                 r200_3716[0] * arcmin_Mpc.value * u.arcmin,
                                 edgecolor='orange',
                                 facecolor='none',
                                 linestyle='--',
                                 transform=ax.get_transform('fk5'),
                                 linewidth=3,
                                 alpha=0.2)

            c2_62_ds = SphericalCircle((coords_3716[1].ra, coords_3716[1].dec) * u.deg,
                                 r200_3716[1] * arcmin_Mpc.value * u.arcmin,
                                 edgecolor='green',
                                 facecolor='none',
                                 linestyle='--',
                                 transform=ax.get_transform('fk5'),
                                 linewidth=3,
                                 alpha=0.2)
            ax_ds.add_patch(c2_58_ds)
            ax_ds.add_patch(c2_62_ds)
        ax.add_patch(c2_58)
        ax.add_patch(c2_62)
    else:
        c2 = SphericalCircle((ra, dec) * u.deg,
                             ab.r200[k] * arcmin_Mpc.value * u.arcmin,
                             edgecolor='black',
                             facecolor='none',
                             linestyle='--',
                             transform=ax.get_transform('fk5'),
                             linewidth=3,
                             alpha=0.2)
        if is_dsp:
            c2_ds = SphericalCircle((ra, dec) * u.deg,
                                 ab.r200[k] * arcmin_Mpc.value * u.arcmin,
                                 edgecolor='black',
                                 facecolor='none',
                                 linestyle='--',
                                 transform=ax.get_transform('fk5'),
                                 linewidth=3,
                                 alpha=0.2)
            ax_ds.add_patch(c2_ds)
            print('after')
        ax.add_patch(c2)
        print('before')



    if k == 3:
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
    elif k == 6:
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

    ins_axes_I_ds = inset_axes(ax_ds,
                            width='100%',
                            height='100%',
                            bbox_to_anchor=(0.05, 0.02, 0.3, 0.1),
                            bbox_transform=ax.transAxes)
    ins_axes_I_ds.set_title(r'% of    (Merging)', fontsize=18)
    ins_axes_I_ds.tick_params(axis='x', which='major', labelsize=0)
    ins_axes_I_ds.set_xticks([])
    ins_axes_PM_ds = inset_axes(ax_ds,
                             width='100%',
                             height='100%',
                             bbox_to_anchor=(0.7, 0.02, 0.3, 0.1),
                             bbox_transform=ax.transAxes)
    ins_axes_PM_ds.set_title(r'% of    (Post-merging)', fontsize=18)
    ins_axes_PM_ds.tick_params(axis='x', which='major', labelsize=0)
    ins_axes_PM_ds.set_xticks([])

    ax.scatter(0.15, 0.134, marker='*', s=220, transform=ax.transAxes, facecolor='none', edgecolor='black', linewidth=3)
    ax.scatter(0.76, 0.134, marker='P', s=200, transform=ax.transAxes, facecolor='none', edgecolor='black', linewidth=3)
    ax_ds.scatter(0.15, 0.134, marker='*', s=220, transform=ax_ds.transAxes, facecolor='none', edgecolor='black', linewidth=3)
    ax_ds.scatter(0.76, 0.134, marker='P', s=200, transform=ax_ds.transAxes, facecolor='none', edgecolor='black', linewidth=3)

    sub_frac_I = np.zeros(this_num_sub)
    sub_frac_PM = np.zeros(this_num_sub)
    sub_frac_Ierr = np.zeros(this_num_sub)
    sub_frac_PMerr = np.zeros(this_num_sub)

    main_total_no = 0
    main_PM_no = 0
    main_I_no = 0

    sub_total_no = 0
    sub_PM_no = 0
    sub_I_no = 0

    main_total_no_3 = 0
    main_PM_no_3 = 0
    main_I_no_3 = 0

    sub_total_no_3 = 0
    sub_PM_no_3 = 0
    sub_I_no_3 = 0

    for i in range(this_num_sub):
       sub_frac_I[i] = num_2[i] / (num_0[i] + num_1[i] + num_2[i])
       sub_frac_Ierr[i] = np.sqrt(num_2[i]) / (num_0[i] + num_1[i] + num_2[i])
       sub_frac_PM[i] = num_1[i] / (num_0[i] + num_1[i] + num_2[i])
       sub_frac_PMerr[i] = np.sqrt(num_1[i]) / (num_0[i] + num_1[i] + num_2[i])
       if k == 3:
           if i == 1 or i == 4:
               main_total_no += (num_0[i] + num_1[i] + num_2[i])
               main_PM_no += num_1[i]
               main_I_no += num_2[i]

               main_total_no_3 += (num_0_3[i] + num_1_3[i] + num_2_3[i])
               main_PM_no_3 += num_1_3[i]
               main_I_no_3 += num_2_3[i]
           else:
               sub_total_no += (num_0[i] + num_1[i] + num_2[i])
               sub_PM_no += num_1[i]
               sub_I_no += num_2[i]

               sub_total_no_3 += (num_0_3[i] + num_1_3[i] + num_2_3[i])
               sub_PM_no_3 += num_1_3[i]
               sub_I_no_3 += num_2_3[i]
       elif k == 0:
           if i == 1:
               main_total_no += (num_0[i] + num_1[i] + num_2[i])
               main_PM_no += num_1[i]
               main_I_no += num_2[i]

               main_total_no_3 += (num_0_3[i] + num_1_3[i] + num_2_3[i])
               main_PM_no_3 += num_1_3[i]
               main_I_no_3 += num_2_3[i]
           else:
               sub_total_no += (num_0[i] + num_1[i] + num_2[i])
               sub_PM_no += num_1[i]
               sub_I_no += num_2[i]

               sub_total_no_3 += (num_0_3[i] + num_1_3[i] + num_2_3[i])
               sub_PM_no_3 += num_1_3[i]
               sub_I_no_3 += num_2_3[i]
       else:
           if i == 0:
               main_total_no += (num_0[i] + num_1[i] + num_2[i])
               main_PM_no += num_1[i]
               main_I_no += num_2[i]

               main_total_no_3 += (num_0_3[i] + num_1_3[i] + num_2_3[i])
               main_PM_no_3 += num_1_3[i]
               main_I_no_3 += num_2_3[i]
           else:
               sub_total_no += (num_0[i] + num_1[i] + num_2[i])
               sub_PM_no += num_1[i]
               sub_I_no += num_2[i]

               sub_total_no_3 += (num_0_3[i] + num_1_3[i] + num_2_3[i])
               sub_PM_no_3 += num_1_3[i]
               sub_I_no_3 += num_2_3[i]

    barlistI = ins_axes_I.bar(range(this_num_sub),
                                    sub_frac_I * 1e2,
                                    width=0.5)
    barlistPM = ins_axes_PM.bar(range(this_num_sub),
                                      sub_frac_PM * 1e2,
                                      width=0.5)
    ins_axes_I.errorbar(range(this_num_sub),
                        sub_frac_Ierr * 1e2,
                        yerr=sub_frac_Ierr * 1e2,
                        fmt='|')
    ins_axes_PM.errorbar(range(this_num_sub),
                         sub_frac_PM * 1e2,
                         yerr=sub_frac_PMerr * 1e2,
                         fmt='|')


    ### fraction comparison between main vs. substructure (mclust)
    sub_frac_I_dsp = np.zeros(this_num_sub_dsp)
    sub_frac_PM_dsp = np.zeros(this_num_sub_dsp)
    sub_frac_Ierr_dsp = np.zeros(this_num_sub_dsp)
    sub_frac_PMerr_dsp = np.zeros(this_num_sub_dsp)

    main_total_no_dsp = 0
    main_PM_no_dsp = 0
    main_I_no_dsp = 0

    sub_total_no_dsp = 0
    sub_PM_no_dsp = 0
    sub_I_no_dsp = 0

    main_total_no_dsp_3 = 0
    main_PM_no_dsp_3 = 0
    main_I_no_dsp_3 = 0

    sub_total_no_dsp_3 = 0
    sub_PM_no_dsp_3 = 0
    sub_I_no_dsp_3 = 0

    ### calculate the fraction in substructure by DS+
    for i in range(this_num_sub_dsp):
        sub_frac_I_dsp[i] = num_2_dsp[i] / (num_0_dsp[i] + num_1_dsp[i] + num_2_dsp[i])
        sub_frac_Ierr_dsp[i] = np.sqrt(num_2_dsp[i]) / (num_0_dsp[i] + num_1_dsp[i] + num_2_dsp[i])
        sub_frac_PM_dsp[i] = num_1_dsp[i] / (num_0_dsp[i] + num_1_dsp[i] + num_2_dsp[i])
        sub_frac_PMerr_dsp[i] = np.sqrt(num_1_dsp[i]) / (num_0_dsp[i] + num_1_dsp[i] + num_2_dsp[i])
        if i == 0:
            main_total_no_dsp += (num_0_dsp[i] + num_1_dsp[i] + num_2_dsp[i])
            main_PM_no_dsp += num_1_dsp[i]
            main_I_no_dsp += num_2_dsp[i]

            main_total_no_dsp_3 += (num_0_dsp_3[i] + num_1_dsp_3[i] + num_2_dsp_3[i])
            main_PM_no_dsp_3 += num_1_dsp_3[i]
            main_I_no_dsp_3 += num_2_dsp_3[i]
        else:
            sub_total_no_dsp += (num_0_dsp[i] + num_1_dsp[i] + num_2_dsp[i])
            sub_PM_no_dsp += num_1_dsp[i]
            sub_I_no_dsp += num_2_dsp[i]

            sub_total_no_dsp_3 += (num_0_dsp_3[i] + num_1_dsp_3[i] + num_2_dsp_3[i])
            sub_PM_no_dsp_3 += num_1_dsp_3[i]
            sub_I_no_dsp_3 += num_2_dsp_3[i]

    barlistI_ds = ins_axes_I_ds.bar(range(this_num_sub_dsp),
                                    sub_frac_I_dsp * 1e2,
                                    width=0.5)
    barlistPM_ds = ins_axes_PM_ds.bar(range(this_num_sub_dsp),
                                      sub_frac_PM_dsp * 1e2,
                                      width=0.5)
    # ins_axes_I_ds.errorbar(range(this_num_sub_dsp),
    #                     sub_frac_Ierr_dsp * 1e2,
    #                     yerr=sub_frac_Ierr_dsp * 1e2,
    #                     fmt='|')
    # ins_axes_PM_ds.errorbar(range(this_num_sub_dsp),
    #                      sub_frac_PM_dsp * 1e2,
    #                      yerr=sub_frac_PMerr_dsp * 1e2,
    #                      fmt='|')

    ### plot mclust and DS+ in one plot
    if k == 3:
        sub_frac_I = [np.mean([sub_frac_I[i] for i in [1, 4]]) - sub_frac_I[ii] for ii in [0, 2, 3, 5]]
        sub_frac_PM = [np.mean([sub_frac_PM[i] for i in [1, 4]]) - sub_frac_PM[ii] for ii in [0, 2, 3, 5]]
    else:
        for i in range(this_num_sub-1):
            sub_frac_I[i+1] = (sub_frac_I[i+1] - sub_frac_I[0]) * 1e2
            sub_frac_PM[i + 1] = (sub_frac_PM[i + 1] - sub_frac_PM[0]) * 1e2

    for i in range(this_num_sub_dsp-1):
        sub_frac_I_dsp[i + 1] = (sub_frac_I_dsp[i + 1] - sub_frac_I_dsp[0]) * 1e2
        sub_frac_PM_dsp[i + 1] = (sub_frac_PM_dsp[i + 1] - sub_frac_PM_dsp[0]) * 1e2

    xerr = ((sub_I_no ** 0.5 / sub_total_no * 100) ** 2 + (main_I_no ** 0.5 / main_total_no * 100) ** 2) ** 0.5
    yerr = ((sub_I_no_dsp ** 0.5 / sub_total_no_dsp * 100) ** 2 + (main_I_no_dsp ** 0.5 / main_total_no_dsp * 100) ** 2) ** 0.5
    axs16.errorbar(sub_I_no / sub_total_no * 100 - main_I_no / main_total_no * 100,
                   sub_I_no_dsp / sub_total_no_dsp * 100 - main_I_no_dsp / main_total_no_dsp * 100,
                   xerr=xerr, yerr=yerr, c='blue', fmt=markers[k], label=ab.clusters[k]+' M', markersize=10)

    xerr = ((sub_PM_no ** 0.5 / sub_total_no * 100) ** 2 + (main_PM_no ** 0.5 / main_total_no * 100) ** 2) ** 0.5
    yerr = ((sub_PM_no_dsp ** 0.5 / sub_total_no_dsp * 100) ** 2 + (main_PM_no_dsp ** 0.5 / main_total_no_dsp * 100) ** 2) ** 0.5
    axs16.errorbar(sub_PM_no / sub_total_no * 100 - main_PM_no / main_total_no * 100,
                   sub_I_no_dsp / sub_total_no_dsp * 100 - main_I_no_dsp / main_total_no_dsp * 100,
                   xerr=xerr, yerr=yerr, c='orange', fmt=markers[k], label=ab.clusters[k]+' PM', markersize=10)

    xerr = ((sub_I_no_3 ** 0.5 / sub_total_no_3 * 100) ** 2 + (main_I_no_3 ** 0.5 / main_total_no_3 * 100) ** 2) ** 0.5
    yerr = ((sub_I_no_dsp_3 ** 0.5 / sub_total_no_dsp_3 * 100) ** 2 + (main_I_no_dsp_3 ** 0.5 / main_total_no_dsp_3 * 100) ** 2) ** 0.5
    axs166.errorbar(sub_I_no_3 / sub_total_no_3 * 100 - main_I_no_3 / main_total_no_3 * 100,
                   sub_I_no_dsp_3 / sub_total_no_dsp_3 * 100 - main_I_no_dsp_3 / main_total_no_dsp_3 * 100,
                   xerr=xerr, yerr=yerr, c='blue', fmt=markers[k], label=ab.clusters[k]+' M', markersize=10)

    xerr = ((sub_PM_no_3 ** 0.5 / sub_total_no_3 * 100) ** 2 + (main_PM_no_3 ** 0.5 / main_total_no_3 * 100) ** 2) ** 0.5
    yerr = ((sub_PM_no_dsp_3 ** 0.5 / sub_total_no_dsp_3 * 100) ** 2 + (main_PM_no_dsp_3 ** 0.5 / main_total_no_dsp_3 * 100) ** 2) ** 0.5
    axs166.errorbar(sub_PM_no_3 / sub_total_no_3 * 100 - main_PM_no_3 / main_total_no_3 * 100,
                   sub_I_no_dsp_3 / sub_total_no_dsp_3 * 100 - main_I_no_dsp_3 / main_total_no_dsp_3 * 100,
                   xerr=xerr, yerr=yerr, c='orange', fmt=markers[k], label=ab.clusters[k]+' PM', markersize=10)

    axs16.scatter(sub_I_no_3 / sub_total_no_3 * 100 - main_I_no_3 / main_total_no_3 * 100,
                   sub_I_no_dsp_3 / sub_total_no_dsp_3 * 100 - main_I_no_dsp_3 / main_total_no_dsp_3 * 100,
                   edgecolors='blue', marker=markers[k], s=100, facecolors='none')
    axs16.scatter(sub_PM_no_3 / sub_total_no_3 * 100 - main_PM_no_3 / main_total_no_3 * 100,
                   sub_I_no_dsp_3 / sub_total_no_dsp_3 * 100 - main_I_no_dsp_3 / main_total_no_dsp_3 * 100,
                   edgecolors='orange', marker=markers[k], s=100, facecolors='none')



    # axs16.errorbar(np.nanmean(sub_frac_I[1:]), np.nanmean(sub_frac_I_dsp[1:]),
    #                yerr=np.nanstd(sub_frac_I[1:])/np.sqrt(this_num_sub-1),
    #                xerr=np.nanstd(sub_frac_I_dsp[1:])/np.sqrt(this_num_sub_dsp-1),
    #                c='blue', fmt=markers[k], label=ab.clusters[k]+' I', markersize=10)
    # axs16.errorbar(np.nanmean(sub_frac_PM[1:]), np.nanmean(sub_frac_PM_dsp[1:]),
    #                    yerr=np.nanstd(sub_frac_PM[1:])/np.sqrt(this_num_sub-1),
    #                    xerr=np.nanstd(sub_frac_PM_dsp[1:])/np.sqrt(this_num_sub_dsp-1),
    #                    c='orange', fmt=markers[k], label=ab.clusters[k]+' PM', markersize=10)



    ### errorbar
    ins_axes_I.set_ylim(0, 35)
    ins_axes_PM.set_ylim(0, 35)
    for i in range(this_num_sub):
        barlistI[i].set_color(colors[i+1])
        barlistPM[i].set_color(colors[i+1])

    ins_axes_I_ds.set_ylim(0, 35)
    ins_axes_PM_ds.set_ylim(0, 35)
    for i in range(this_num_sub_dsp):
        barlistI_ds[i].set_color(colors[i % len(colors) - 1])   # black for background is at the last of the color array
        barlistPM_ds[i].set_color(colors[i % len(colors) - 1])

    ax.text(0.03, 0.95, f'{ab.clusters[k]}', transform=ax.transAxes, fontsize=25)
    if k == 3:
        ax1_58.text(0.03, 0.95, 'a3558', transform=ax.transAxes, fontsize=25)
        ax1_62.text(0.03, 0.95, 'a3562', transform=ax.transAxes, fontsize=25)
    elif k == 6:
        ax1_58.text(0.03, 0.95, 'a3716S', transform=ax.transAxes, fontsize=25)
        ax1_62.text(0.03, 0.95, 'a3716N', transform=ax.transAxes, fontsize=25)
    else:
        ax1.text(0.03, 0.95, f'{ab.clusters[k]}', transform=ax.transAxes, fontsize=25)
    if is_dsp:
        ax_ds.text(0.03, 0.95, f'{ab.clusters[k]}', transform=ax.transAxes, fontsize=25)
        ax1_ds.text(0.03, 0.95, f'{ab.clusters[k]}', transform=ax.transAxes, fontsize=25)

    if k == 3:
        for i in range(len(coords_ss)):
            ax.scatter(coords_ss[i].ra.value, coords_ss[i].dec.value, marker='x', s=150, c='black',
                       transform=ax.get_transform('fk5'))
            ax.text(coords_ss[i].ra.value, coords_ss[i].dec.value, cluster_ss[i], c='black',
                    transform=ax.get_transform('fk5'))
            if is_dsp:
                ax_ds.scatter(coords_ss[i].ra.value, coords_ss[i].dec.value, marker='x', s=150, c='black',
                           transform=ax.get_transform('fk5'))
                ax_ds.text(coords_ss[i].ra.value, coords_ss[i].dec.value, cluster_ss[i], c='black',
                        transform=ax.get_transform('fk5'))

    if k != 3 and k != 6:
        ax1.set_ylim([-4, 4])

    ax1_ds.set_ylim([-4, 4])
    ax.invert_xaxis()
    ax_ds.invert_xaxis()
    fig.savefig(ab.work_dir + 'plot/' + ab.clusters[k] + f'_rosat_substruc.png')
    fig_ds.savefig(ab.work_dir + 'plot/' + ab.clusters[k] + f'_rosat_substruc_ds_{plim_p}.png')
    fig1_ds.savefig(ab.work_dir + 'plot/' + ab.clusters[k] + f'_phase_ds_{plim_p}.png')

    # draw PPS region lines
    if k == 3:
        draw_phase_region(ax1_58, 0.3)
        draw_escape_line(ax1_58, 0.3, r200_ss[0], m200_ss[0], rad_vel_sig_ss[0])
        draw_phase_region(ax1_62, 0.3)
        draw_escape_line(ax1_62, 0.3, r200_ss[1], m200_ss[1], rad_vel_sig_ss[1])
        fig1_58.savefig(ab.work_dir + 'plot/' + ab.clusters[k] + '_phase_58.png')
        fig1_62.savefig(ab.work_dir + 'plot/' + ab.clusters[k] + '_phase_62.png')
    elif k == 6:
        draw_phase_region(ax1_58, 0.3)
        draw_escape_line(ax1_58, 0.3, r200_3716[0], m200_3716[0], rad_vel_sig_3716[0])
        draw_phase_region(ax1_62, 0.3)
        draw_escape_line(ax1_62, 0.3, r200_3716[1], m200_3716[1], rad_vel_sig_3716[1])
        fig1_58.savefig(ab.work_dir + 'plot/' + ab.clusters[k] + '_phase_S.png')
        fig1_62.savefig(ab.work_dir + 'plot/' + ab.clusters[k] + '_phase_N.png')
    else:
        draw_phase_region(ax1, 0.3)
        draw_escape_line(ax1, 0.3, ab.r200[k], ab.m200[k], ab.sig_vr[k])
        fig1.savefig(ab.work_dir + 'plot/' + ab.clusters[k] + '_phase.png')


barWidth = 0.25
num_sub[2] += 0.5
ax2.bar(num_sub, frac_I, label='I', width=barWidth)
ax2.bar([x + 0.25 for x in num_sub], frac_PM, label='PM', width=barWidth)
ax2.set_xlabel('Number of substructure')
ax2.set_ylabel('Fraction')
ax2.legend()
fig2.savefig(ab.work_dir + f'plot/frac_vs_subnum_{plim_p}.png')


axs12.set_xlim(-0.5, 3.0)
axs12.set_ylim(-4, 4)

img_sizes2 = []

zoom = 0.2

for i in range(len(images)):
    imagebox = OffsetImage(images[i], zoom=zoom)
    abox = AnnotationBbox(imagebox, xy=img_locs[i], xycoords='data', frameon=False, pad=0)
    extent = abox.get_window_extent(axs12.figure.canvas.get_renderer())
    extent_data = extent.transformed(axs12.transData.inverted())
    x_size = abs(extent_data.x1 - extent_data.x0)
    y_size = abs(extent_data.y1 - extent_data.y0)
    img_sizes2.append((x_size, y_size))

ddx = 0.6
ddy = 0.9
for i in range(len(images)):
    imagebox = OffsetImage(images[i], zoom=zoom)
    # imagebox.image.axes = axs12
    x, y = img_locs[i]
    width_this, height_this = img_sizes2[i]
    is_broken = False
    dxx, dyy = (0, 0)
    for dx, dy in [(0, 0), (0, ddy/2), (0, -ddy/2), (0, ddy), (0, -ddy), (0, ddy*1.5), (0, -ddy*1.5), (ddx/2, 0), (-ddx/2, 0),
                   (ddx/2, ddy/2), (ddx/2, -ddy/2), (-ddx/2, ddy/2), (-ddx/2, -ddy/2),
                   (ddx/2, ddy), (ddx/2, -ddy), (-ddx/2, ddy), (-ddx/2, -ddy),
                   (ddx / 2, ddy*1.5), (ddx / 2, -ddy*1.5), (-ddx / 2, ddy*1.5), (-ddx / 2, -ddy*1.5),
                   (ddx, ddy / 2), (ddx, -ddy / 2), (ddx, ddy), (ddx, -ddy),
                   (ddx, ddy * 1.5), (ddx, -ddy * 1.5), (ddx*1.5, 0), (ddx * 1.5, ddy * 1.5),
                   (-ddx, ddy * 1.5), (-ddx, -ddy * 1.5), (-ddx, 0), (-ddx, ddy * 1.5),
                   (ddx, ddy * 2.0), (ddx, -ddy * 2.0), (ddx * 2.0, 0), (ddx * 2.0, ddy * 2.0),
                   (ddx, ddy * 3.0), (ddx, -ddy * 3.0), (ddx * 3.0, 0), (ddx * 3.0, ddy * 3.0)
                   ]:
        xs, ys = zip(*[loc for j, loc in enumerate(img_locs) if j != i])
        widths, heights = zip(*[size for j, size in enumerate(img_sizes2) if j != i])
        if all((abs((x + dx) - xi) > (width/2 + width_this/2)) or (abs((y + dy) - yi) > (height/2 + height_this/2))
               for xi, yi, width, height in zip(xs, ys, widths, heights)):
            x, y, dxx, dyy = x, y, dx, dy
            # print(f'i:{i}, x:{x}, dx:{dx}, y:{y}, dy:{dy}')
            is_broken = True
            break
    if is_broken:
        print(f'i:{i}, x:{x}, width:{width_this}, y:{y}, height:{height_this}')
    else:
        print(f'not broken!! i:{i}, x:{x}, width:{width_this}, y:{y}, height:{height_this}')
    img_locs[i] = (x + dxx, y + dyy)
    abox = AnnotationBbox(imagebox, (x, y), xybox=(x + dxx, y + dyy), boxcoords="data", xycoords='data',
                          bboxprops=dict(edgecolor='blue' if is_I[i] else 'orange', linewidth=5),
                          arrowprops=dict(arrowstyle="->", alpha=0.2), pad=0)
    axs12.add_artist(abox)


axs12.axvline(x=0.5, linestyle='--', color='grey')
axs12.axvline(x=1.0, linestyle='--', color='grey')
axs12.set_xticks([0.0, 0.5, 1.0, 2.0])
# axs12.set_xticklabels(['in', 'R200', 'out'], fontsize=20)
axs12.tick_params(axis='both', which='major', labelsize=30)
axs12.tick_params(direction='in')
# axs12.set_xscale('log')

axs12.set_xlabel(r'R/R$_{200}$', fontsize=30)
axs12.set_ylabel(r'$\Delta$v$_{r}$ / $\sigma_{v_{r}}$', fontsize=30, labelpad=30)

blue_box = patches.Rectangle((0, 0), 1, 1, facecolor='blue')
orange_box = patches.Rectangle((0, 0), 1, 1, facecolor='orange')

# Create a custom legend with the two boxes
axs12.legend([blue_box, orange_box], ['Interacting', 'Post merging'], loc='lower right', fontsize=30)
fig12.savefig(ab.plot_dir + f'frac_I_vs_r_{plim_p}.png')

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


# Create a custom legend with the two boxes
# axs122.legend()
# axs122.yaxis.tick_right()
# axs122.yaxis.set_ticks_position('both')
# axs122.xaxis.tick_top()
# axs122.xaxis.set_ticks_position('both')
# axs122.axis.set_ticks_position('both')
fig122.savefig(ab.plot_dir + f'frac_I_vs_r_{plim_p}_infall.png')

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
fig16.savefig(ab.plot_dir + f'frac_vs_sub_all_{plim_p}_total_frac.png')

## for feature types selected by 3 out of 4 people
axs166.set_title(r'%$_{substructures}$ - %$_{main}$', fontsize=20)
axs166.text(0.65, 0.95, 'more in sub', fontsize=20, transform=axs16.transAxes)
# axs16.text(21.5, 16.5, 'sub', fontsize=20)
axs166.text(0.01, 0.01, 'more in main', fontsize=20, transform=axs16.transAxes)
# axs16.text(-10, -6.5, 'in main', fontsize=20)
#                       fontsize=10, rotation=60)
axs166.tick_params(axis='both', which='major', labelsize=15)
axs166.tick_params(direction='in')
axs166.axhline(0, linestyle='--', color='grey')
axs166.axvline(0, linestyle='--', color='grey')
axs166.set_ylabel('DS+', fontsize=20)
axs166.set_xlabel('mclust', fontsize=20)
axs166.legend(fontsize=10, loc='upper left')
fig166.savefig(ab.plot_dir + f'frac_vs_sub_all_{plim_p}_total_frac_3agree.png')

plt.close('all')