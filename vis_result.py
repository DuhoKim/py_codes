import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
# import matplotlib
# matplotlib.use('TkAgg')     # https://stackoverflow.com/questions/28757348/how-to-clear-memory-completely-of-all-matplotlib-plots
from reproject.mosaicking import find_optimal_celestial_wcs, reproject_and_coadd
from reproject import reproject_interp
from astropy.visualization.wcsaxes import SphericalCircle
from astropy import units as u
import abell_cluster_module as ab
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.visualization import make_lupton_rgb
from astropy.io import ascii
from astropy import wcs
import my_module as mm
import img_scale
import importlib
importlib.reload(img_scale)
importlib.reload(mm)
import pylab as py
import sys
from astropy.convolution import Gaussian2DKernel
from scipy.signal import convolve as scipy_convolve
from astropy.convolution import convolve
from astropy.cosmology import Planck18 as Cosmo
from scipy import constants as const
import os
from astropy.table import vstack

clusters = ['A754', 'A2399', 'A2670', 'A3558', 'A3716']
redshifts = [0.05450, 0.05790, 0.07619, 0.04800, 0.04620]
sigmas = [920, 730, 750, 1000, 850]
r200 = [2.2, 1.7, 1.8, 2.4, 2.0]
m200 = [10.2, 6.5, 5.3, 16.9, 10.2]
sb = [28.6, 28.7, 28.8, 28.7, 28.3]
pec_vel = [32, 1518, 352, 1, 903]

coords_cl = [SkyCoord('09:09:08.4 -09:39:58', unit=(u.hourangle, u.deg)),
                 SkyCoord('21:57:25.8 -07:47:41', unit=(u.hourangle, u.deg)),
                 SkyCoord('23:54:13.7 -10:25:08', unit=(u.hourangle, u.deg)),
                 SkyCoord('13:27:57.5 -31:30:09', unit=(u.hourangle, u.deg)),
                 SkyCoord('20:51:16 -52:41.7', unit=(u.hourangle, u.deg))]

save_dir = ab.work_dir + f'plot/'

mag_gap = np.zeros(len(clusters)) # brighteset - penaltimatum

tot_cnt = np.zeros((len(clusters), 5))  # count for total sample, pm, interacting, pair, faint companion
spec_cnt = np.zeros((len(clusters), 5))  # count for spec member
rs_cnt = np.zeros((len(clusters), 5))  # count for red sequence

fig = plt.figure(figsize=(5, 5))

rate1 = []
rate2 = []
rate3 = []
rate4 = []

spec_rate1 = []
spec_rate2 = []
spec_rate3 = []
spec_rate4 = []

rs_rate1 = []
rs_rate2 = []
rs_rate3 = []
rs_rate4 = []

star_class_star = []
galaxy_class_star = []

star_spread_model = []
galaxy_spread_model = []

is_spread = True

# for k in range(0, len(clusters)):
for k in range(3, 4):
    with open(ab.work_dir + f'spec/{clusters[k]}_spec_ind.npy', 'rb') as spec_ind, \
        open(ab.work_dir+f'spec/{clusters[k]}_rs_each_ind.npy', 'rb') as rs_ind,  \
        open(ab.work_dir+f'pics/{clusters[k]}_gala_back/{clusters[k]}_rs_or_spec_duho.vis', 'r') as vis:

        if is_spread:
            psf_sex_dir = ab.work_dir + f'sex/run/{clusters[k]}_rsi_'
            sex_cat1 = ascii.read(psf_sex_dir + '1_spread.cat')
            sex_cat2 = ascii.read(psf_sex_dir + '2_spread.cat')
            sex_cat3 = ascii.read(psf_sex_dir + '3_spread.cat')
            sex_cat4 = ascii.read(psf_sex_dir + '4_spread.cat')
            sex_cat5 = ascii.read(psf_sex_dir + '5_spread.cat')
            sex_cat6 = ascii.read(psf_sex_dir + '6_spread.cat')
            sex_cat7 = ascii.read(psf_sex_dir + '7_spread.cat')
            sex_cat8 = ascii.read(psf_sex_dir + '8_spread.cat')
            sex_cat9 = ascii.read(psf_sex_dir + '9_spread.cat')

            spr_sex_cat = vstack([sex_cat1, sex_cat2, sex_cat3, sex_cat4, sex_cat5, sex_cat6, sex_cat7, sex_cat8, sex_cat9])

        sex_cat = ascii.read(ab.sex_dir + f'DECam_merged_SEx_cat_{clusters[k]}_Gal_ext_corrected_{ab.ver}.txt')
        ind_spec = np.load(spec_ind)
        ind_rs = np.load(rs_ind)
        ind_total = ind_spec | ind_rs

        if is_spread:
            coords_spr = SkyCoord(spr_sex_cat['ALPHA_J2000'], spr_sex_cat['DELTA_J2000'], unit='deg')
        coords_sex = SkyCoord(sex_cat['ALPHA_J2000'], sex_cat['DELTA_J2000'], unit='deg')

        if is_spread:
            idx_spr2sex, d2d_sex2spr, d3d = coords_sex.match_to_catalog_sky(coords_spr)
            match_sex2spr = (d2d_sex2spr.arcsec < ab.max_sep)

        r_mag_spec = sex_cat[ind_spec]
        r_mag_spec.sort('MAG_AUTO_r')

        mag_gap[k] = r_mag_spec['MAG_AUTO_r'][1] - r_mag_spec['MAG_AUTO_r'][0]

        print(f"{r_mag_spec['MAG_AUTO_r'][0]} {r_mag_spec['MAG_AUTO_r'][1]}")

        tot_rad = []
        tot_int = []
        spec_rad = []
        spec_int = []
        rs_rad = []
        rs_int = []

        for line in vis:
            item = line.split(' ')
            if len(item) > 2:
                if int(item[1]) != 4:
                    # galaxy_class_star.append(sex_cat['CLASS_STAR'][int(item[0])-1])
                    if is_spread:
                        galaxy_spread_model.append(spr_sex_cat[idx_spr2sex][int(item[0])-1]['SPREAD_MODEL'])
                    galaxy_class_star.append(spr_sex_cat[idx_spr2sex][int(item[0]) - 1]['CLASS_STAR'])
                    tot_cnt[k][0] += 1
                    item[0] = int(item[0])
                    item[2] = int(item[2])
                    coord = SkyCoord(sex_cat['ALPHA_J2000'][item[0]-1], sex_cat['DELTA_J2000'][item[0]-1], unit='deg')
                    dist = mm.clcendist(redshifts[k], redshifts[k], coords_cl[k], coord)
                    tot_rad.append(dist/r200[k])
                    if item[2] & 2 ** 7:
                        tot_cnt[k][2] += 1
                        tot_int.append(1)
                    else:
                        tot_int.append(0)
                    if (item[2] & 2 ** 6) or (item[2] & 2 ** 5):
                        tot_cnt[k][1] += 1
                    if item[2] & 2 ** 8:
                        tot_cnt[k][3] += 1
                    if item[2] & 2 ** 9:
                        tot_cnt[k][4] += 1
                    if item[0] in sex_cat['NUMBER'][ind_spec]:
                        spec_cnt[k][0] += 1
                        spec_rad.append(dist/r200[k])
                        if item[2] & 2 ** 7:
                            spec_cnt[k][2] += 1
                            spec_int.append(1)
                        else:
                            spec_int.append(0)
                        if (item[2] & 2 ** 6) or (item[2] & 2 ** 5):
                            spec_cnt[k][1] += 1
                        if item[2] & 2 ** 8:
                            spec_cnt[k][3] += 1
                        if item[2] & 2 ** 9:
                            spec_cnt[k][4] += 1
                    if item[0] in sex_cat['NUMBER'][ind_rs]:
                        rs_cnt[k][0] += 1
                        rs_rad.append(dist/r200[k])
                        if item[2] & 2 ** 7:
                            rs_cnt[k][2] += 1
                            rs_int.append(1)
                        else:
                            rs_int.append(0)
                        if (item[2] & 2 ** 6) or (item[2] & 2 ** 5):
                            rs_cnt[k][1] += 1
                        if item[2] & 2 ** 8:
                            rs_cnt[k][3] += 1
                        if item[2] & 2 ** 9:
                            rs_cnt[k][4] += 1
                else:
                    # star_class_star.append(sex_cat['CLASS_STAR'][int(item[0]) - 1])
                    if is_spread:
                        star_spread_model.append(spr_sex_cat[idx_spr2sex][int(item[0]) - 1]['SPREAD_MODEL'])
                    star_class_star.append(spr_sex_cat[idx_spr2sex][int(item[0]) - 1]['CLASS_STAR'])

        print(f'{clusters[k]} tot_cnt:{tot_cnt[k][0]}, tot_pm_cnt:{tot_cnt[k][1]}, tot_int_cnt:{tot_cnt[k][2]} , tot_pair_cnt:{tot_cnt[k][3]} , tot_fc_cnt:{tot_cnt[k][4]}')
        print(f' rs_cnt:{rs_cnt[k][0]}, rs_pm_cnt:{rs_cnt[k][1]}, rs_int_cnt:{rs_cnt[k][2]}, rs_pair_cnt:{rs_cnt[k][3]}, rs_fc_cnt:{rs_cnt[k][4]}')
        print(f' spec_cnt:{spec_cnt[k][0]}, spec_pm_cnt:{spec_cnt[k][1]}, spec_int_cnt:{spec_cnt[k][2]}, spec_pair_cnt:{spec_cnt[k][3]}, spec_fc_cnt:{spec_cnt[k][4]}')

        tot_rad = np.array(tot_rad)
        spec_rad = np.array(spec_rad)
        rs_rad = np.array(rs_rad)

        tot_int = np.array(tot_int)
        spec_int = np.array(spec_int)
        rs_int = np.array(rs_int)

        rad_bin_cen = []

        tot_bin = []
        tot_err_bin = []
        spec_bin = []
        spec_err_bin = []
        rs_bin = []
        rs_err_bin = []

        rad_bin_size = 0.5

        for i in range(0, 4):
            rad_bin_cen.append(rad_bin_size * i + rad_bin_size/2)

            lower = rad_bin_cen[i] - rad_bin_size/2
            upper = rad_bin_cen[i] + rad_bin_size/2

            in_bin_tot = (tot_rad > lower) & (tot_rad < upper)
            in_bin_spec = (spec_rad > lower) & (spec_rad < upper)
            in_bin_rs = (rs_rad > lower) & (rs_rad < upper)

            if i == 0:
                rate1.extend(tot_int[in_bin_tot])
                spec_rate1.extend(spec_int[in_bin_spec])
                rs_rate1.extend(rs_int[in_bin_rs])
            elif i == 1:
                rate2.extend(tot_int[in_bin_tot])
                spec_rate2.extend(spec_int[in_bin_spec])
                rs_rate2.extend(rs_int[in_bin_rs])
            elif i == 2:
                rate3.extend(tot_int[in_bin_tot])
                spec_rate3.extend(spec_int[in_bin_spec])
                rs_rate3.extend(rs_int[in_bin_rs])
            else:
                rate4.extend(tot_int[in_bin_tot])
                spec_rate4.extend(spec_int[in_bin_spec])
                rs_rate4.extend(rs_int[in_bin_rs])

            tot_bin.append(np.sum(tot_int[in_bin_tot])/len(tot_int[in_bin_tot])*100)
            tot_err_bin.append(np.sqrt(np.sum(tot_int[in_bin_tot])) / len(tot_int[in_bin_tot]) * 100)
            spec_bin.append(np.sum(spec_int[in_bin_spec]) / len(spec_int[in_bin_spec]) * 100)
            spec_err_bin.append(np.sqrt(np.sum(spec_int[in_bin_spec])) / len(spec_int[in_bin_spec]) * 100)
            rs_bin.append(np.sum(rs_int[in_bin_rs]) / len(rs_int[in_bin_rs]) * 100)
            rs_err_bin.append(np.sqrt(np.sum(rs_int[in_bin_rs])) / len(rs_int[in_bin_rs]) * 100)

        plt.plot(rad_bin_cen, tot_bin, label='Total', color='g')
        plt.plot(rad_bin_cen, spec_bin, label='Spec', color='b', linestyle='--', alpha=0.8)
        plt.plot(rad_bin_cen, rs_bin, label='RS', color='r', linestyle='--', alpha=0.8)
        # plt.scatter(rad_bin_cen, tot_bin, label='Total', color='g')
        # plt.scatter(rad_bin_cen, spec_bin, label='Spec', color='b')
        # plt.scatter(rad_bin_cen, rs_bin, label='RS', color='r')
        # plt.errorbar(rad_bin_cen, spec_bin, yerr=spec_err_bin)
        plt.errorbar(rad_bin_cen, tot_bin, yerr=tot_err_bin, color='g')
        # plt.errorbar(rad_bin_cen, rs_bin, yerr=rs_err_bin)

plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.savefig(save_dir + f'inter_frac_vs_radius.png')
plt.close(fig)

fig = plt.figure(figsize=(5, 5))
plt.hist(np.array(star_class_star), color='black', histtype='step')
plt.hist(np.array(galaxy_class_star), color='green', histtype='step')
plt.show()

bins = np.array(range(-10, 10)) * 0.005
fig = plt.figure(figsize=(5, 5))
plt.hist(np.array(star_spread_model), color='black', histtype='step', bins=bins)
plt.hist(np.array(galaxy_spread_model), color='green', histtype='step', bins=bins)
plt.show()

fig = plt.figure(figsize=(5, 5))
r1 = np.sum(rate1)/len(rate1)*100
r2 = np.sum(rate2)/len(rate2)*100
r3 = np.sum(rate3)/len(rate3)*100
r4 = np.sum(rate4)/len(rate4)*100
tot_rates = [r1, r2, r3, r4]
tot_ints = [np.sum(rate1), np.sum(rate2), np.sum(rate3), np.sum(rate4)]
tot_lens = [len(rate1), len(rate2), len(rate3), len(rate4)]
sr1 = np.sum(spec_rate1)/len(spec_rate1)*100
sr2 = np.sum(spec_rate2)/len(spec_rate2)*100
sr3 = np.sum(spec_rate3)/len(spec_rate3)*100
sr4 = np.sum(spec_rate4)/len(spec_rate4)*100
spec_rates = [sr1, sr2, sr3, sr4]
spec_ints = [np.sum(spec_rate1), np.sum(spec_rate2), np.sum(spec_rate3), np.sum(spec_rate4)]
spec_lens = [len(spec_rate1), len(spec_rate2), len(spec_rate3), len(spec_rate4)]
rr1 = np.sum(rs_rate1)/len(rs_rate1)*100
rr2 = np.sum(rs_rate2)/len(rs_rate2)*100
rr3 = np.sum(rs_rate3)/len(rs_rate3)*100
rr4 = np.sum(rs_rate4)/len(rs_rate4)*100
rs_rates = [rr1, rr2, rr3, rr4]
rs_ints = [np.sum(rs_rate1), np.sum(rs_rate2), np.sum(rs_rate3), np.sum(rs_rate4)]
rs_lens = [len(rs_rate1), len(rs_rate2), len(rs_rate3), len(rs_rate4)]
re1 = np.sqrt(np.sum(rate1))/len(rate1)*100
re2 = np.sqrt(np.sum(rate2))/len(rate2)*100
re3 = np.sqrt(np.sum(rate3))/len(rate3)*100
re4 = np.sqrt(np.sum(rate4))/len(rate4)*100
sre1 = np.sqrt(np.sum(spec_rate1))/len(spec_rate1)*100
sre2 = np.sqrt(np.sum(spec_rate2))/len(spec_rate2)*100
sre3 = np.sqrt(np.sum(spec_rate3))/len(spec_rate3)*100
sre4 = np.sqrt(np.sum(spec_rate4))/len(spec_rate4)*100
rre1 = np.sqrt(np.sum(rs_rate1))/len(rs_rate1)*100
rre2 = np.sqrt(np.sum(rs_rate2))/len(rs_rate2)*100
rre3 = np.sqrt(np.sum(rs_rate3))/len(rs_rate3)*100
rre4 = np.sqrt(np.sum(rs_rate4))/len(rs_rate4)*100
plt.plot(rad_bin_cen, tot_rates, label='Total', color='g')
plt.plot(np.array(rad_bin_cen)+0.02, spec_rates, label='Spec', color='b')
plt.plot(np.array(rad_bin_cen)-0.02, rs_rates, label='RS', color='r')
plt.errorbar(rad_bin_cen, [r1, r2, r3, r4], yerr=[re1, re2, re3, re4], fmt='o', color='g')
plt.errorbar(np.array(rad_bin_cen)+0.02, [sr1, sr2, sr3, sr4], yerr=[sre1, sre2, sre3, sre4], fmt='o', color='b')
plt.errorbar(np.array(rad_bin_cen)-0.02, [rr1, rr2, rr3, rr4], yerr=[rre1, rre2, rre3, rre4], fmt='o', color='r')
for k in range(0, 4):
    plt.text(rad_bin_cen[k]+0.02, tot_rates[k]-0.5, f'{tot_ints[k]}/{tot_lens[k]}', fontsize=12, color='g')
    plt.text(rad_bin_cen[k]+0.04, spec_rates[k], f'{spec_ints[k]}/{spec_lens[k]}', fontsize=12, color='b')
    plt.text(rad_bin_cen[k], rs_rates[k]+0.5, f'{rs_ints[k]}/{rs_lens[k]}', fontsize=12, color='r')
plt.xlabel(r'R/R$_{200}$', fontsize=16)
plt.ylabel('Fraction of Interacting galaxy [%]', fontsize=16)
plt.legend(fontsize=16, loc='upper left')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlim([0, 2])
plt.savefig(save_dir + f'inter_frac_vs_radius_all.png')
plt.close(fig)


fig = plt.figure(figsize=(5, 5))
plt.scatter(mag_gap, tot_cnt[:,2]/tot_cnt[:,0]*100, label='Total', color='g')
plt.errorbar(mag_gap, tot_cnt[:,2]/tot_cnt[:,0]*100, yerr = np.sqrt(tot_cnt[:,2])/tot_cnt[:,0]*100, fmt='o', color='g')
plt.scatter(mag_gap, spec_cnt[:,2]/spec_cnt[:,0]*100, label='Spec', color='b')
plt.scatter(mag_gap, rs_cnt[:,2]/rs_cnt[:,0]*100, label='Red Sequence', color='r')
for k in range(0, len(clusters)):
    plt.text(mag_gap[k], tot_cnt[k,2]/tot_cnt[k,0]*100, f'{clusters[k]}', fontsize=12)
    plt.text(mag_gap[k], tot_cnt[k, 2] / tot_cnt[k, 0] * 100 - 0.3, f'z={redshifts[k]}', fontsize=12)
plt.xlabel(r'm_{brightest} - m_{2nd brightest}', fontsize=16)
plt.ylabel('Fraction of Interacting galaxy [%]', fontsize=16)
plt.legend(fontsize=12)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.savefig(save_dir + f'inter_frac_vs_mag_gap.png')
plt.close(fig)

fig = plt.figure(figsize=(5, 5))
plt.scatter(m200, tot_cnt[:,2]/tot_cnt[:,0]*100, label='Total', color='g')
plt.errorbar(m200, tot_cnt[:,2]/tot_cnt[:,0]*100, yerr = np.sqrt(tot_cnt[:,2])/tot_cnt[:,0]*100, fmt='o')
plt.scatter(m200, spec_cnt[:,2]/spec_cnt[:,0]*100, label='Spec', color='b')
plt.scatter(m200, rs_cnt[:,2]/rs_cnt[:,0]*100, label='Red Sequence', color='r')
for k in range(0, len(clusters)):
    plt.text(m200[k], tot_cnt[k,2]/tot_cnt[k,0]*100, f'{clusters[k]}', fontsize=12)
    plt.text(m200[k], tot_cnt[k, 2] / tot_cnt[k, 0] * 100 - 0.3, f'z={redshifts[k]}', fontsize=12)
plt.xlabel(r'$m_{200}$ [10$^{14}$ $M\odot$]', fontsize=14)
plt.ylabel('Fraction of Interacting galaxy [%]', fontsize=16)
plt.legend(fontsize=12)
plt.xticks(fontsize=14)
plt.yticks(fontsize=16)
plt.savefig(save_dir + f'inter_frac_vs_m200.png')
plt.close(fig)

fig = plt.figure(figsize=(5, 5))
plt.scatter(redshifts, tot_cnt[:,2]/tot_cnt[:,0]*100, label='Total', color='g')
plt.errorbar(redshifts, tot_cnt[:,2]/tot_cnt[:,0]*100, yerr = np.sqrt(tot_cnt[:,2])/tot_cnt[:,0]*100, fmt='o')
plt.scatter(redshifts, spec_cnt[:,2]/spec_cnt[:,0]*100, label='Spec', color='b')
plt.scatter(redshifts, rs_cnt[:,2]/rs_cnt[:,0]*100, label='Red Sequence', color='r')
for k in range(0, len(clusters)):
    plt.text(redshifts[k], tot_cnt[k,2]/tot_cnt[k,0]*100, f'{clusters[k]}', fontsize=12)
    plt.text(redshifts[k], tot_cnt[k, 2] / tot_cnt[k, 0] * 100 - 0.3, f'z={redshifts[k]}', fontsize=12)
plt.xlabel(r'$z$', fontsize=16)
plt.ylabel('Fraction of Interacting galaxy [%]', fontsize=16)
plt.legend(fontsize=12)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.savefig(save_dir + f'inter_frac_vs_z.png')
plt.close(fig)

fig = plt.figure(figsize=(5, 5))
plt.scatter(pec_vel, tot_cnt[:,2]/tot_cnt[:,0]*100, label='Total', color='g')
plt.errorbar(pec_vel, tot_cnt[:,2]/tot_cnt[:,0]*100, yerr = np.sqrt(tot_cnt[:,2])/tot_cnt[:,0]*100, fmt='o')
plt.scatter(pec_vel, spec_cnt[:,2]/spec_cnt[:,0]*100, label='Spec', color='b')
plt.scatter(pec_vel, rs_cnt[:,2]/rs_cnt[:,0]*100, label='Red Sequence', color='r')
for k in range(0, len(clusters)):
    plt.text(pec_vel[k], tot_cnt[k,2]/tot_cnt[k,0]*100, f'{clusters[k]}', fontsize=12)
    plt.text(pec_vel[k], tot_cnt[k, 2] / tot_cnt[k, 0] * 100 - 0.3, f'z={redshifts[k]}', fontsize=12)
plt.xlabel(r'$Peculiar velocity$', fontsize=16)
plt.ylabel('Fraction of Interacting galaxy [%]', fontsize=16)
plt.legend(fontsize=12)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.savefig(save_dir + f'inter_frac_vs_pec_vel.png')
plt.close(fig)



# fig = plt.figure(figsize=(5, 5))
# plt.scatter(sb, tot_cnt[:,2]/tot_cnt[:,0]*100, label='Total', color='g')
# plt.scatter(sb, spec_cnt[:,2]/spec_cnt[:,0]*100, label='Spec', color='b')
# plt.scatter(sb, rs_cnt[:,2]/rs_cnt[:,0]*100, label='Red Sequence', color='r')
# for k in range(0, len(clusters)):
#     plt.text(sb[k], tot_cnt[k,2]/tot_cnt[k,0]*100, f'{clusters[k]}')
#     plt.text(sb[k], tot_cnt[k, 2] / tot_cnt[k, 0] * 100 - 0.3, f'z={redshifts[k]}')
# plt.xlabel(r'$SB_{100}$')
# plt.ylabel('Fraction of the Interacting galaxy [%]')
# plt.legend()
# plt.savefig(save_dir + f'inter_frac_vs_sb.png')
# plt.close(fig)
#
#
# fig = plt.figure(figsize=(5, 5))
# plt.scatter(mag_gap, tot_cnt[:,1]/tot_cnt[:,0]*100, label='Total', color='g')
# plt.scatter(mag_gap, spec_cnt[:,1]/spec_cnt[:,0]*100, label='Spec', color='b')
# plt.scatter(mag_gap, rs_cnt[:,1]/rs_cnt[:,0]*100, label='Red Sequence', color='r')
# for k in range(0, len(clusters)):
#     plt.text(mag_gap[k], tot_cnt[k,1]/tot_cnt[k,0]*100, f'{clusters[k]}')
#     plt.text(mag_gap[k], tot_cnt[k, 1] / tot_cnt[k, 0] * 100 - 0.3, f'z={redshifts[k]}')
# plt.xlabel(r'm_{brightest} - m_{2nd brightest}')
# plt.ylabel('Fraction of the Post-Merger galaxy [%]')
# plt.legend()
# plt.savefig(save_dir + f'pm_frac_vs_mag_gap.png')
# plt.close(fig)
#
# fig = plt.figure(figsize=(5, 5))
# plt.scatter(mag_gap, tot_cnt[:,3]/tot_cnt[:,0]*100, label='Total', color='g')
# plt.scatter(mag_gap, spec_cnt[:,3]/spec_cnt[:,0]*100, label='Spec', color='b')
# plt.scatter(mag_gap, rs_cnt[:,3]/rs_cnt[:,0]*100, label='Red Sequence', color='r')
# for k in range(0, len(clusters)):
#     plt.text(mag_gap[k], tot_cnt[k,3]/tot_cnt[k,0]*100, f'{clusters[k]}')
#     plt.text(mag_gap[k], tot_cnt[k, 3] / tot_cnt[k, 0] * 100 - 0.3, f'z={redshifts[k]}')
# plt.xlabel(r'm_{brightest} - m_{2nd brightest}')
# plt.ylabel('Fraction of the Paired galaxy [%]')
# plt.legend()
# plt.savefig(save_dir + f'p_frac_vs_mag_gap.png')
# plt.close(fig)
#
# fig = plt.figure(figsize=(5, 5))
# plt.scatter(mag_gap, tot_cnt[:,4]/tot_cnt[:,0]*100, label='Total', color='g')
# plt.scatter(mag_gap, spec_cnt[:,4]/spec_cnt[:,0]*100, label='Spec', color='b')
# plt.scatter(mag_gap, rs_cnt[:,4]/rs_cnt[:,0]*100, label='Red Sequence', color='r')
# for k in range(0, len(clusters)):
#     plt.text(mag_gap[k], tot_cnt[k,4]/tot_cnt[k,0]*100, f'{clusters[k]}')
#     plt.text(mag_gap[k], tot_cnt[k, 4] / tot_cnt[k, 0] * 100 - 0.3, f'z={redshifts[k]}')
# plt.xlabel(r'm_{brightest} - m_{2nd brightest}')
# plt.ylabel('Fraction of the galaxy w/ Faint Companion [%]')
# plt.legend()
# plt.savefig(save_dir + f'fc_frac_vs_mag_gap.png')
# plt.close(fig)
#
#
#
#
# fig = plt.figure(figsize=(5, 5))
# plt.scatter(sigmas, tot_cnt[:,2]/tot_cnt[:,0]*100, label='Total', color='g')
# plt.scatter(sigmas, spec_cnt[:,2]/spec_cnt[:,0]*100, label='Spec', color='b')
# plt.scatter(sigmas, rs_cnt[:,2]/rs_cnt[:,0]*100, label='Red Sequence', color='r')
# for k in range(0, len(clusters)):
#     plt.text(sigmas[k], tot_cnt[k,2]/tot_cnt[k,0]*100, f'{clusters[k]}')
#     plt.text(sigmas[k], tot_cnt[k, 2] / tot_cnt[k, 0] * 100 - 0.3, f'z={redshifts[k]}')
# plt.xlabel(r'$\sigma$ [km s-1]')
# plt.ylabel('Fraction of the Interacting galaxy [%]')
# plt.legend()
# plt.savefig(save_dir + f'inter_frac_vs_sigma.png')
# plt.close(fig)
#
# fig = plt.figure(figsize=(5, 5))
# plt.scatter(sigmas, tot_cnt[:,1]/tot_cnt[:,0]*100, label='Total', color='g')
# plt.scatter(sigmas, spec_cnt[:,1]/spec_cnt[:,0]*100, label='Spec', color='b')
# plt.scatter(sigmas, rs_cnt[:,1]/rs_cnt[:,0]*100, label='Red Sequence', color='r')
# for k in range(0, len(clusters)):
#     plt.text(sigmas[k], tot_cnt[k,1]/tot_cnt[k,0]*100, f'{clusters[k]}')
#     plt.text(sigmas[k], tot_cnt[k, 1] / tot_cnt[k, 0] * 100 - 0.3, f'z={redshifts[k]}')
# plt.xlabel(r'$\sigma$ [km s-1]')
# plt.ylabel('Fraction of the Post-Merger galaxy [%]')
# plt.legend()
# plt.savefig(save_dir + f'pm_frac_vs_sigma.png')
# plt.close(fig)
#
# fig = plt.figure(figsize=(5, 5))
# plt.scatter(sigmas, tot_cnt[:,3]/tot_cnt[:,0]*100, label='Total', color='g')
# plt.scatter(sigmas, spec_cnt[:,3]/spec_cnt[:,0]*100, label='Spec', color='b')
# plt.scatter(sigmas, rs_cnt[:,3]/rs_cnt[:,0]*100, label='Red Sequence', color='r')
# for k in range(0, len(clusters)):
#     plt.text(sigmas[k], tot_cnt[k,3]/tot_cnt[k,0]*100, f'{clusters[k]}')
#     plt.text(sigmas[k], tot_cnt[k, 3] / tot_cnt[k, 0] * 100 - 0.3, f'z={redshifts[k]}')
# plt.xlabel(r'$\sigma$ [km s-1]')
# plt.ylabel('Fraction of the Paired galaxy [%]')
# plt.legend()
# plt.savefig(save_dir + f'p_frac_vs_sigma.png')
# plt.close(fig)
#
# fig = plt.figure(figsize=(5, 5))
# plt.scatter(sigmas, tot_cnt[:,4]/tot_cnt[:,0]*100, label='Total', color='g')
# plt.scatter(sigmas, spec_cnt[:,4]/spec_cnt[:,0]*100, label='Spec', color='b')
# plt.scatter(sigmas, rs_cnt[:,4]/rs_cnt[:,0]*100, label='Red Sequence', color='r')
# for k in range(0, len(clusters)):
#     plt.text(sigmas[k], tot_cnt[k,4]/tot_cnt[k,0]*100, f'{clusters[k]}')
#     plt.text(sigmas[k], tot_cnt[k, 4] / tot_cnt[k, 0] * 100 - 0.3, f'z={redshifts[k]}')
# plt.xlabel(r'$\sigma$ [km s-1]')
# plt.ylabel('Fraction of the galaxy w/ Faint Companion [%]')
# plt.legend()
# plt.savefig(save_dir + f'fc_frac_vs_sigma.png')
# plt.close(fig)
#
#
#
#
#
# fig = plt.figure(figsize=(5, 5))
# plt.scatter(m200, tot_cnt[:,1]/tot_cnt[:,0]*100, label='Total', color='g')
# plt.scatter(m200, spec_cnt[:,1]/spec_cnt[:,0]*100, label='Spec', color='b')
# plt.scatter(m200, rs_cnt[:,1]/rs_cnt[:,0]*100, label='Red Sequence', color='r')
# for k in range(0, len(clusters)):
#     plt.text(m200[k], tot_cnt[k,1]/tot_cnt[k,0]*100, f'{clusters[k]}')
#     plt.text(m200[k], tot_cnt[k, 1] / tot_cnt[k, 0] * 100 - 0.3, f'z={redshifts[k]}')
# plt.xlabel(r'$m_{200}$ [10$^{14}$ $M\odot$]')
# plt.ylabel('Fraction of the Post-Merger galaxy [%]')
# plt.legend()
# plt.savefig(save_dir + f'pm_frac_vs_m200.png')
# plt.close(fig)
#
# fig = plt.figure(figsize=(5, 5))
# plt.scatter(m200, tot_cnt[:,3]/tot_cnt[:,0]*100, label='Total', color='g')
# plt.scatter(m200, spec_cnt[:,3]/spec_cnt[:,0]*100, label='Spec', color='b')
# plt.scatter(m200, rs_cnt[:,3]/rs_cnt[:,0]*100, label='Red Sequence', color='r')
# for k in range(0, len(clusters)):
#     plt.text(m200[k], tot_cnt[k,3]/tot_cnt[k,0]*100, f'{clusters[k]}')
#     plt.text(m200[k], tot_cnt[k, 3] / tot_cnt[k, 0] * 100 - 0.3, f'z={redshifts[k]}')
# plt.xlabel(r'$m_{200}$ [10$^{14}$ $M\odot$]')
# plt.ylabel('Fraction of the Paired galaxy [%]')
# plt.legend()
# plt.savefig(save_dir + f'p_frac_vs_m200.png')
# plt.close(fig)
#
# fig = plt.figure(figsize=(5, 5))
# plt.scatter(m200, tot_cnt[:,4]/tot_cnt[:,0]*100, label='Total', color='g')
# plt.scatter(m200, spec_cnt[:,4]/spec_cnt[:,0]*100, label='Spec', color='b')
# plt.scatter(m200, rs_cnt[:,4]/rs_cnt[:,0]*100, label='Red Sequence', color='r')
# for k in range(0, len(clusters)):
#     plt.text(m200[k], tot_cnt[k,4]/tot_cnt[k,0]*100, f'{clusters[k]}')
#     plt.text(m200[k], tot_cnt[k, 4] / tot_cnt[k, 0] * 100 - 0.3, f'z={redshifts[k]}')
# plt.xlabel(r'$m_{200}$ [10$^{14}$ $M\odot$]')
# plt.ylabel('Fraction of the galaxy w/ Faint Companion [%]')
# plt.legend()
# plt.savefig(save_dir + f'fc_frac_vs_m200.png')
# plt.close(fig)
#
#
#
#
#
# fig = plt.figure(figsize=(5, 5))
# plt.scatter(redshifts, tot_cnt[:,1]/tot_cnt[:,0]*100, label='Total', color='g')
# plt.scatter(redshifts, spec_cnt[:,1]/spec_cnt[:,0]*100, label='Spec', color='b')
# plt.scatter(redshifts, rs_cnt[:,1]/rs_cnt[:,0]*100, label='Red Sequence', color='r')
# for k in range(0, len(clusters)):
#     plt.text(redshifts[k], tot_cnt[k,1]/tot_cnt[k,0]*100, f'{clusters[k]}')
#     plt.text(redshifts[k], tot_cnt[k, 1] / tot_cnt[k, 0] * 100 - 0.3, f'z={redshifts[k]}')
# plt.xlabel(r'$z$')
# plt.ylabel('Fraction of the Post-Merger galaxy [%]')
# plt.legend()
# plt.savefig(save_dir + f'pm_frac_vs_z.png')
# plt.close(fig)
#
# fig = plt.figure(figsize=(5, 5))
# plt.scatter(redshifts, tot_cnt[:,3]/tot_cnt[:,0]*100, label='Total', color='g')
# plt.scatter(redshifts, spec_cnt[:,3]/spec_cnt[:,0]*100, label='Spec', color='b')
# plt.scatter(redshifts, rs_cnt[:,3]/rs_cnt[:,0]*100, label='Red Sequence', color='r')
# for k in range(0, len(clusters)):
#     plt.text(redshifts[k], tot_cnt[k,3]/tot_cnt[k,0]*100, f'{clusters[k]}')
#     plt.text(redshifts[k], tot_cnt[k, 3] / tot_cnt[k, 0] * 100 - 0.3, f'z={redshifts[k]}')
# plt.xlabel(r'$z$')
# plt.ylabel('Fraction of the Paired galaxy [%]')
# plt.legend()
# plt.savefig(save_dir + f'p_frac_vs_z.png')
# plt.close(fig)
#
# fig = plt.figure(figsize=(5, 5))
# plt.scatter(redshifts, tot_cnt[:,4]/tot_cnt[:,0]*100, label='Total', color='g')
# plt.scatter(redshifts, spec_cnt[:,4]/spec_cnt[:,0]*100, label='Spec', color='b')
# plt.scatter(redshifts, rs_cnt[:,4]/rs_cnt[:,0]*100, label='Red Sequence', color='r')
# for k in range(0, len(clusters)):
#     plt.text(redshifts[k], tot_cnt[k,4]/tot_cnt[k,0]*100, f'{clusters[k]}')
#     plt.text(redshifts[k], tot_cnt[k, 4] / tot_cnt[k, 0] * 100 - 0.3, f'z={redshifts[k]}')
# plt.xlabel(r'$z$')
# plt.ylabel('Fraction of the galaxy w/ Faint Companion [%]')
# plt.legend()
# plt.savefig(save_dir + f'fc_frac_vs_z.png')
# plt.close(fig)
#
#
#
#
#
# fig = plt.figure(figsize=(5, 5))
# plt.scatter(sb, tot_cnt[:,1]/tot_cnt[:,0]*100, label='Total', color='g')
# plt.scatter(sb, spec_cnt[:,1]/spec_cnt[:,0]*100, label='Spec', color='b')
# plt.scatter(sb, rs_cnt[:,1]/rs_cnt[:,0]*100, label='Red Sequence', color='r')
# for k in range(0, len(clusters)):
#     plt.text(sb[k], tot_cnt[k,1]/tot_cnt[k,0]*100, f'{clusters[k]}')
#     plt.text(sb[k], tot_cnt[k, 1] / tot_cnt[k, 0] * 100 - 0.3, f'z={redshifts[k]}')
# plt.xlabel(r'$SB_{100}$')
# plt.ylabel('Fraction of the Post-Merger galaxy [%]')
# plt.legend()
# plt.savefig(save_dir + f'pm_frac_vs_sb.png')
# plt.close(fig)