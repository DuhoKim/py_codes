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

fov_ra = [ 1.5, 3.0, 3.0, 3.0, 0.0, 0.0, 2.0]
fov_dec = [ 1.5, 3.0, 3.0, 2.5, 0.0, 0.0, 1.2]

colors = ['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink']

fig2 = plt.figure(figsize=(5, 5))
ax2 = fig2.add_subplot()

fig10, axs10 = plt.subplots(tight_layout=True, figsize=(6, 6))
fig11, axs11 = plt.subplots(tight_layout=True, figsize=(6, 6))
fig12 = plt.figure(tight_layout=True, figsize=(12, 12))
axs12 = plt.subplot2grid((6, 6), (1, 1), rowspan=4, colspan=4)
fig13, axs13 = plt.subplots(tight_layout=True, figsize=(6, 6))

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

coords_ss = [SkyCoord('13:27:54 -31:30:00', unit=(u.hourangle, u.deg)),
                 SkyCoord('13:33:30 -31:40:00', unit=(u.hourangle, u.deg))]

cluster_ss = ['A3558', 'A3562']

r200_ss = [2.09, 2.43]
xx_ss = [0.25, 0.75, 1.5]

num_sub = [] # number of substructure
frac_I = []   # percentage of interacting
frac_PM = []  # percentage of post-merger


d_int_tot = []
d_pm_tot = []
d_eith_tot = []


local_peak = []

num_I_in = 0
num_I_out = 0
num_I_mid = 0


# for k in range(0, len(clusters)):
for k in [0, 1, 2, 3, 6]:
#for k in [0, 1, 2, 3]:
# for k in [3]:
    mm.print_time()
    print(f'{ab.clusters[k]}')

    fig = plt.figure(figsize=(10, 10))
    fig1 = plt.figure(figsize=(10, 10))
    ax1 = fig1.add_subplot()

    vis_cat = ascii.read(
        ab.work_dir + f'sex/cat/DECam_merged_SEx_cat_{ab.clusters[k]}_Gal_ext_corrected_20rmag_psf.txt')
    spec_cat = ascii.read(ab.work_dir + f'spec/{ab.clusters[k]}_spec_ra_dec_z.txt')
    r_cat = ascii.read(ab.work_dir + f'spec/{ab.clusters[k]}_class_uncer.csv', format='csv')
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

    # if k == 3:
    #     dist_vis2 = coords_spec.separation(coords_x[k][0]).value  # in degrees
    #     dist_vis3 = coords_spec.separation(coords_x[k][1]).value  # in degrees
    #     dist_vis4 = coords_spec.separation(coords_x[k][2]).value  # in degrees
    #     dist_vis5 = coords_spec.separation(coords_x[k][3]).value  # in degrees
    # else:
    #     dist_vis = coords_spec.separation(coords_x[k][0]).value  # in degrees

    ra, dec = ab.coords_cl_cen[k].ra.deg, ab.coords_cl_cen[k].dec.deg

    with fits.open(ab.work_dir + 'fits/rosat/' + ab.clusters[k] + '.fits') as hdu_x:
        wcs_x, shape_x = find_optimal_celestial_wcs(hdu_x)  # has only CompImageHDU files
        ax = fig.add_subplot(projection=wcs_x)
        data = hdu_x[0].data
        sigma = 1.2  # this depends on how noisy your data is, play with it!
        data = gaussian_filter(data, sigma)
        ax.contour(data, projection=wcs_x, alpha=0.3, linewidths=3)

        x_min, y_min = wcs_x.all_world2pix(ra - fov_ra[k] / 2, dec - fov_dec[k] / 2, 0)
        x_max, y_max = wcs_x.all_world2pix(ra + fov_ra[k] / 2, dec + fov_dec[k] / 2, 0)

    z = ab.redshifts[k]  # redshift
    kpc_arcmin = Cosmo.kpc_proper_per_arcmin(z)  # xxx kpc / arcmin
    arcmin_Mpc = 1e3 / kpc_arcmin  # xxx arcmin / 5 Mpc
    r200_deg = ab.r200[k] * 1e3 / kpc_arcmin / 6e1

    if k == 3:
        r200_deg_58 = r200_ss[0] * 1e3 / kpc_arcmin / 6e1
        r200_deg_62 = r200_ss[1] * 1e3 / kpc_arcmin / 6e1
        r_in_r200_58 = coords_spec.separation(coords_ss[0]).value / r200_deg_58.value
        r_in_r200_62 = coords_spec.separation(coords_ss[1]).value / r200_deg_62.value
    else:
        r_in_r200 = coords_spec.separation(ab.coords_cl_cen[k]).value / r200_deg.value

    num_0 = np.zeros(this_num_sub)
    num_1 = np.zeros(this_num_sub)
    num_2 = np.zeros(this_num_sub)

    # number count from clustocentric r200
    # n_three_cen = [0 for x in range(6)]  # number counts for three for in and out of I,PM,other from clustocentric dist
    if k == 3:
        n_three_cen58 = [0 for x in range(9)]  # number counts for three for in and out of I,PM,other from clustocentric dist
        n_three_cen62 = [0 for x in range(9)]  # number counts for three for in and out of I,PM,other from clustocentric dist
    else:
        n_three_cen = [0 for x in range(9)]  # number counts for three for in and out of I,PM,other from clustocentric dist

    for i in range(len(r_cat)):
        this_mk = '.'
        this_size = 20
        this_alpha = 0.3
        if matched_spec[i]:
            if vis_cat[idx_vis_spec[i]]['NUMBER'] in res_cat['col1']:
                ind, = np.where(res_cat['col1'] == vis_cat[idx_vis_spec[i]]['NUMBER'])
                if res_cat['col2'][ind] == 0:  # classified as no feature
                    this_mk = 'o'
                    this_alpha = 0.4
                    num_0[r_cat[i][0]-1] += 1
                    if k == 3:
                        if r_cat[i][0] == 2:
                            if r_in_r200_58[i] < 0.5:
                                n_three_cen58[0] += 1
                            elif r_in_r200_58[i] < 1.0:
                                n_three_cen58[1] += 1
                            else:
                                n_three_cen58[2] += 1
                        elif r_cat[i][0] == 5:
                            if r_in_r200_62[i] < 0.5:
                                n_three_cen62[0] += 1
                            elif r_in_r200_62[i] < 1.0:
                                n_three_cen62[1] += 1
                            else:
                                n_three_cen62[2] += 1
                    else:
                        if r_in_r200[i] < 0.5:
                            n_three_cen[0] += 1
                        elif r_in_r200[i] < 1.0:
                            n_three_cen[1] += 1
                        else:
                            n_three_cen[2] += 1
                elif res_cat['col2'][ind] == 1:  # classified as PM
                    this_mk = 'P'
                    this_size = 80
                    this_alpha = 0.6
                    num_1[r_cat[i][0]-1] += 1
                    if k == 3:
                        if r_cat[i][0] == 2:
                            if r_in_r200_58[i] < 0.5:
                                n_three_cen58[3] += 1
                            elif r_in_r200_58[i] < 1.0:
                                n_three_cen58[4] += 1
                            else:
                                n_three_cen58[5] += 1
                        elif r_cat[i][0] == 5:
                            if r_in_r200_62[i] < 0.5:
                                n_three_cen62[3] += 1
                            elif r_in_r200_62[i] < 1.0:
                                n_three_cen62[4] += 1
                            else:
                                n_three_cen62[5] += 1
                    else:
                        if r_in_r200[i] < 0.5:
                            n_three_cen[3] += 1
                        elif r_in_r200[i] < 1.0:
                            n_three_cen[4] += 1
                        else:
                            n_three_cen[5] += 1
                elif res_cat['col2'][ind] == 2:  # classified as I
                    this_mk = '*'
                    this_size = 120
                    this_alpha = 0.6
                    num_2[r_cat[i][0]-1] += 1
                    img = mpimg.imread(ab.work_dir + f'pics/{ab.clusters[k]}_for_all/'
                                                     f'{res_cat["col1"][ind][0]}_crop.png')
                    #aspect_ratio = img.shape[0] / img.shape[1]
                    if k == 3:
                        if r_cat[i][0] == 2:
                            if r_in_r200_58[i] < 0.5:
                                n_three_cen58[6] += 1
                            elif r_in_r200_58[i] < 1.0:
                                n_three_cen58[7] += 1
                                png_ax = plt.subplot2grid((6, 6), (0, 1 + num_I_mid), rowspan=1, colspan=1, fig=fig12)
                                png_ax.imshow(img)
                                png_ax.axis('off')
                                num_I_mid += 1
                                # line = plt.Line2D([0.15 * num_I_mid + 0.1, 0.5], [0.9, 0.6], transform=fig12.transFigure,
                                #                   color='r', linestyle=':')
                                rect = patches.Rectangle((0, 0), img.shape[1], img.shape[0], linewidth=10,
                                                         edgecolor='r', facecolor='none')
                                png_ax.add_patch(rect)
                            else:
                                n_three_cen58[8] += 1
                            print(f"sub:{r_cat[i][0]}, R/R200:{r_in_r200_58[i]}, id:{res_cat['col1'][ind][0]}")
                        elif r_cat[i][0] == 5:
                            if r_in_r200_62[i] < 0.5:
                                n_three_cen62[6] += 1
                                png_ax = plt.subplot2grid((6, 6), (5, 1), rowspan=1, colspan=1, fig=fig12)
                                # line = plt.Line2D([0.2, 0.2], [0.1, 0.2], transform=fig12.transFigure,
                                #                   color='violet', linestyle=':')
                                rect = patches.Rectangle((0, 0), img.shape[1], img.shape[0], linewidth=10,
                                                         edgecolor='violet', facecolor='none')
                                png_ax.add_patch(rect)
                                png_ax.imshow(img)
                                png_ax.axis('off')
                                num_I_in += 1
                            elif r_in_r200_62[i] < 1.0:
                                n_three_cen62[7] += 1
                            else:
                                n_three_cen62[8] += 1
                            print(f"sub:{r_cat[i][0]}, R/R200:{r_in_r200_62[i]}, id:{res_cat['col1'][ind][0]}")
                    else:
                        if r_in_r200[i] < 0.5:
                            n_three_cen[6] += 1
                            if k == 0:
                                num_grid = num_I_in + 3
                            else:
                                num_grid = num_I_in - 3
                            png_ax = plt.subplot2grid((6, 6), (num_grid, 0), rowspan=1, colspan=1, fig=fig12)
                            # line = plt.Line2D([0.1, 0.2], [0.9 - 0.2 * num_I_in, 0.5], transform=fig12.transFigure,
                            #                   color=colors[k], linestyle=':')
                            png_ax.imshow(img)
                            rect = patches.Rectangle((0, 0), img.shape[1], img.shape[0], linewidth=10,
                                                     edgecolor=colors[k], facecolor='none')
                            png_ax.add_patch(rect)
                            png_ax.axis('off')
                            num_I_in += 1
                        elif r_in_r200[i] < 1.0:
                            n_three_cen[7] += 1
                        else:
                            n_three_cen[8] += 1
                            if k == 1:
                                num_grid = num_I_out + 5
                            else:
                                num_grid = num_I_out + 1
                            png_ax = plt.subplot2grid((6, 6), (num_grid, 5), rowspan=1, colspan=1, fig=fig12)
                            # line = plt.Line2D([0.9, 0.8], [0.9 - 0.2 * num_I_out - 0.4, 0.1], transform=fig12.transFigure,
                            #                   color=colors[k], linestyle=':')
                            rect = patches.Rectangle((0, 0), img.shape[1], img.shape[0], linewidth=10,
                                                     edgecolor=colors[k], facecolor='none')
                            png_ax.add_patch(rect)
                            png_ax.imshow(img)
                            png_ax.axis('off')
                            num_I_out += 1
                        print(f'sub:{r_cat[i][0]}, R/R200:{r_in_r200[i]}, id:{res_cat["col1"][ind][0]}')
                    # fig12.lines.extend([line])

        ax.scatter(spec_cat['col2'][i], spec_cat['col3'][i], transform=ax.get_transform('fk5'),
                   marker=this_mk,
                   s=this_size,
                   c=colors[r_cat[i][0]],
                   alpha=this_alpha,
                   linewidth=1)

        coord = SkyCoord(ra=spec_cat['col2'][i], dec=spec_cat['col3'][i], unit='deg')
        sep = ab.coords_cl_cen[k].separation(coord)
        # dist = Distance(sep, unit='deg')

        ax1.scatter(sep/r200_deg, spec_cat['col4'][i],
                    marker=this_mk,
                    s=this_size,
                    c=colors[r_cat[i][0]],
                    alpha=this_alpha,
                    linewidth=1)

    ax.set_xlabel('R.A.', fontsize=25)
    ax.set_ylabel('Decl.', fontsize=25, labelpad=-1)
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.tick_params(direction='in')

    ax1.set_xlabel('R/R200', fontsize=25)
    ax1.set_ylabel('z', fontsize=25, labelpad=-1)
    ax1.tick_params(axis='both', which='major', labelsize=20)
    # ax1.tick_params(direction='in')

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)

    c2 = SphericalCircle((ra, dec) * u.deg,
                        ab.r200[k] * arcmin_Mpc.value * u.arcmin,
                        edgecolor='black',
                        facecolor='none',
                        linestyle='--',
                        transform=ax.get_transform('fk5'),
                        linewidth=3,
                        alpha=0.2)
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

        ax.add_patch(c2_58)
        ax.add_patch(c2_62)

    ax.add_patch(c2)

    if k == 3:
        xx58 = []
        xx62 = []
        pm58 = []
        pme58 = []
        pm62 = []
        pme62 = []
        i58 = []
        ie58 = []
        i62 = []
        ie62 = []
        for i in range(3):
            y58 = n_three_cen58[i] + n_three_cen58[i + 3] + n_three_cen58[i + 6]
            print(f'i58:{i}, no:{n_three_cen58[i]}, pm:{n_three_cen58[i + 3]}, I:{n_three_cen58[i + 6]}')
            y62 = n_three_cen62[i] + n_three_cen62[i + 3] + n_three_cen62[i + 6]
            print(f'i62:{i}, no:{n_three_cen62[i]}, pm:{n_three_cen62[i + 3]}, I:{n_three_cen62[i + 6]}')
            if y58 > 1:
                xx58.append(xx_ss[i])
                x = n_three_cen58[i + 3]
                y = y58
                pm58.append(x / y * 100)
                # calculate error of the percentage using Jacobian matrices
                #dz = (100 / y) * np.sqrt(x) + (-100 * x / y ** 2) * (np.sqrt(y))
                dz = (100 / y) * np.sqrt(x) if x > 0 else (100 / y)
                pme58.append(dz)
                x = n_three_cen58[i + 6]
                i58.append(x / y * 100)
                dz = (100 / y) * np.sqrt(x) if x > 0 else (100 / y)
                ie58.append(dz)
            if y62 > 1:
                xx62.append(xx_ss[i])
                x = n_three_cen62[i + 3]
                y = y62
                pm62.append(x / y * 100)
                # calculate error of the percentage using Jacobian matrices
                dz = (100 / y) * np.sqrt(x) if x > 0 else (100 / y)
                pme62.append(dz)
                x = n_three_cen62[i + 6]
                i62.append(x / y * 100)
                dz = (100 / y) * np.sqrt(x) + (-100 * x / y ** 2) * (np.sqrt(y))
                ie62.append(dz)
        axs12.plot(xx58, i58, linestyle='solid', c=colors[k], label='a3558')
        axs12.fill_between(xx58,
                           [a - b for a, b in zip(i58, ie58)],
                           [a + b for a, b in zip(i58, ie58)],
                           color=colors[k], alpha=.15)
        axs13.plot(xx58, pm58, linestyle='dashed', c=colors[k], label='a3558')
        axs13.fill_between(xx58,
                           [a - b for a, b in zip(pm58, pme58)],
                           [a + b for a, b in zip(pm58, pme58)],
                           color=colors[k], alpha=.15)
        axs12.plot(xx62, i62, linestyle='solid', c=colors[k+1], label='a3562')
        axs12.fill_between(xx62,
                           [a - b for a, b in zip(i62, ie62)],
                           [a + b for a, b in zip(i62, ie62)],
                           color=colors[k+1], alpha=.15)
        axs13.plot(xx62, pm62, linestyle='dashed', c=colors[k+1], label='a3562')
        axs13.fill_between(xx62,
                           [a - b for a, b in zip(pm62, pme62)],
                           [a + b for a, b in zip(pm62, pme62)],
                           color=colors[k+1], alpha=.15)
    else:
        xx = []
        pm = []
        pme = []
        ii = []
        iie = []
        for i in range(3):
            y = n_three_cen[i] + n_three_cen[i + 3] + n_three_cen[i + 6]
            print(f'i:{i}, no:{n_three_cen[i]}, pm:{n_three_cen[i + 3]}, I:{n_three_cen[i + 6]}')
            if y > 2:
                xx.append(xx_ss[i])
                x = n_three_cen[i + 3]
                pm.append(x / y * 100)
                dz = (100 / y) * np.sqrt(x) if x > 0 else (100 / y)
                pme.append(dz)
                x = n_three_cen[i + 6]
                ii.append(x / y * 100)
                dz = (100 / y) * np.sqrt(x) if x > 0 else (100 / y)
                iie.append(dz)

        axs12.plot(xx, ii, linestyle='solid', c=colors[k], label=ab.clusters[k])
        axs12.fill_between(xx,
                           [a - b for a, b in zip(ii, iie)],
                           [a + b for a, b in zip(ii, iie)],
                           color=colors[k], alpha=.15)
        axs13.plot(xx, pm, linestyle='dashed', c=colors[k], label=ab.clusters[k])
        axs13.fill_between(xx,
                           [a - b for a, b in zip(pm, pme)],
                           [a + b for a, b in zip(pm, pme)],
                           color=colors[k], alpha=.15)



        # axs12.text(0.2, f_i_in, f'{n_three_cen[4]}/{num_tot_in}', fontsize=15, c=colors[k])
        # axs12.text(1.6, f_i_out, f'{n_three_cen[5]}/{num_tot_out}', fontsize=15, c=colors[k])
        # axs12.text(0.2, f_pm_in, f'{n_three_cen[2]}/{num_tot_in}', fontsize=15, c=colors[k])
        # axs12.text(1.6, f_pm_out, f'{n_three_cen[3]}/{num_tot_out}', fontsize=15, c=colors[k])

    ins_axes_I = inset_axes(ax,
                            width='100%',
                            height='100%',
                            bbox_to_anchor=(0.05, 0.02, 0.3, 0.1),
                            bbox_transform=ax.transAxes)
    ins_axes_I.set_title(r'% of      (Interacting)')
    ins_axes_I.tick_params(axis='x', which='major', labelsize=0)
    ins_axes_I.set_xticks([])
    ins_axes_PM = inset_axes(ax,
                             width='100%',
                             height='100%',
                             bbox_to_anchor=(0.7, 0.02, 0.3, 0.1),
                             bbox_transform=ax.transAxes)
    ins_axes_PM.set_title(r'% of      (Post-merging)')
    ins_axes_PM.tick_params(axis='x', which='major', labelsize=0)
    ins_axes_PM.set_xticks([])

    ax.scatter(0.15, 0.13, marker='*', s=120, c='gray', transform=ax.transAxes)
    ax.scatter(0.785, 0.13, marker='P', s=80, c='gray', transform=ax.transAxes)

    sub_frac_I = np.zeros(this_num_sub)
    sub_frac_PM = np.zeros(this_num_sub)

    for i in range(this_num_sub):
       sub_frac_I[i] = num_2[i] / (num_0[i] + num_1[i] + num_2[i])
       sub_frac_PM[i] = num_1[i] / (num_0[i] + num_1[i] + num_2[i])

    barlistI = ins_axes_I.bar(range(this_num_sub),
                             sub_frac_I * 1e2,
                             width=0.5)
    barlistPM = ins_axes_PM.bar(range(this_num_sub),
                             sub_frac_PM * 1e2,
                             width=0.5)
    # ins_axes_I.set_ylim(0, 35)
    # ins_axes_PM.set_ylim(0, 35)
    for i in range(this_num_sub):
        barlistI[i].set_color(colors[i+1])
        barlistPM[i].set_color(colors[i+1])

    ax.text(0.03, 0.95, f'{ab.clusters[k]}', transform=ax.transAxes, fontsize=25)
    ax1.text(0.03, 0.95, f'{ab.clusters[k]}', transform=ax.transAxes, fontsize=25)

    if k == 3:
        for i in range(len(coords_ss)):
            ax.scatter(coords_ss[i].ra.value, coords_ss[i].dec.value, marker='x', s=150, c='black',
                       transform=ax.get_transform('fk5'))
            ax.text(coords_ss[i].ra.value, coords_ss[i].dec.value, cluster_ss[i], c='black',
                    transform=ax.get_transform('fk5'))

    # plt.plot([1000, 1000 + pixel_Mpc.value], [1000, 1000], linewidth=3, color='purple')
    # plt.text(1000, 1500, f'1Mpc at z={z:5.3f}', fontsize=25)

    fig.savefig(ab.work_dir + 'plot/' + ab.clusters[k] + '_rosat_substruc.png')
    fig1.savefig(ab.work_dir + 'plot/' + ab.clusters[k] + '_phase.png')

barWidth = 0.25
num_sub[2] += 0.5
ax2.bar(num_sub, frac_I, label='I', width=barWidth)
ax2.bar([x + 0.25 for x in num_sub], frac_PM, label='PM', width=barWidth)
ax2.set_xlabel('Number of substructure')
ax2.set_ylabel('Fraction')
ax2.legend()
fig2.savefig(ab.work_dir + 'plot/frac_vs_subnum.png')

axs11.set_ylabel(r'd$_{norm}$', fontsize=20)
# plt.setp(axs11, xticks=[1, 2, 3], xticklabels=['I', 'PM', 'I+PM'])
axs11.set_xticklabels(['I', 'PM', 'I+PM'], fontsize=20)
axs11.tick_params(axis='y', labelsize=15)
axs11.axhline(1, linestyle='--', color='grey')
axs11.boxplot([d_int_tot, d_pm_tot, d_eith_tot])
axs11.legend(ncol=3, loc='upper right', fontsize=8)

fig11.savefig(ab.plot_dir + f'clcendist_comp_xray.png')

# axs12.axvline(x=1, linestyle='--', color='grey')
axs12.set_xticks([0.25, 0.75, 1.5])
# axs12.set_xticklabels(['in', 'R200', 'out'], fontsize=20)
axs12.tick_params(axis='both', which='major', labelsize=15)
axs12.tick_params(direction='in')
axs12.set_xlim(0.1, 1.6)
axs12.set_ylim(0, 20)
axs12.set_xlabel(r'R/R$_{200}$', fontsize=20)
axs12.set_ylabel('% of Interacting gals', fontsize=20)
axs12.legend(fontsize=15, loc='upper left')
fig12.savefig(ab.plot_dir + f'frac_I_vs_r.png')

axs13.set_xticks([0.25, 0.75, 1.5])
# axs12.set_xticklabels(['in', 'R200', 'out'], fontsize=20)
axs13.tick_params(axis='both', which='major', labelsize=15)
axs13.tick_params(direction='in')
axs13.set_xlim(0.1, 1.6)
axs13.set_ylim(0, 30)
axs13.set_xlabel(r'R/R$_{200}$', fontsize=20)
axs13.set_ylabel('% of gals w/ Post-merging feat.', fontsize=20)
axs13.legend(fontsize=15, loc='upper left')
fig13.savefig(ab.plot_dir + f'frac_PM_vs_r.png')

plt.close('all')