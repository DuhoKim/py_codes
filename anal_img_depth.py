import numpy as np
from astropy import units as u
from astropy.io import ascii
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.modeling.models import Sersic2D, Gaussian2D
import multiprocessing as mp
import abell_cluster_module as ab
import my_module as mm
import random
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import time
import glob
from multiprocessing import Pool
import importlib
importlib.reload(ab)

ANAL_TOT = True

surveys = [
    ['stack', 'short'],
    ['stack', 'short', 'SDSS', 'DECaLS'],
    ['stack', 'short', 'SDSS', 'DECaLS'],
    ['stack', 'short'],
    ['stack', 'short'],
    ['stack', 'short'],
    ['stack', 'short', 'DECaLS']]

# surveys = [
#     [],
#     [],
#     [],
#     [],
#     [],
#     [],
#     ['DECaLS']]

mag_ran = [20, 28]
offset = 3.5

fs = 20     # fontsize
ext = '_2028'

cur_dir = os.getcwd()

NUM_CPUS = None

def find_90_1(mag_arr, comp_arr):
    comp_gt_90 = np.where(comp_arr - 0.9 > 0)[0]    # index where completeness > .9
    comp_lt_1 = np.where(comp_arr - 0.01 < 0)[0]  # index where completeness < .1

    a = comp_arr[comp_gt_90][-1] - 0.9
    b = 0.9 - comp_arr[comp_gt_90[-1]+1]
    e = mag_arr[1] - mag_arr[0]
    dx = a * e / (a + b)

    m90 = mag_arr[comp_gt_90][-1] + dx

    if len(comp_lt_1) == 0: # comps are all greater than 1%
        return m90, 27.5
    else:
        a = comp_arr[comp_lt_1[0]-1] - 0.01
        b = 0.01 - comp_arr[comp_lt_1[0]]
        e = mag_arr[1] - mag_arr[0]
        dx = a * e / (a + b)

        m1 = mag_arr[comp_lt_1[0]-1] + dx

        return m90, m1


def mock_worker(kk):
    mm.print_time()
    print(f'{ab.clusters[l]} {ab.bands[m]} {survey} MOSAIC {kk} tile '
          f'{(mag_ran[1] - mag_ran[0]) / mag_in_step} obj start!')
    fn = f'{ab.clusters[l]}_{ab.bands[m]}si_{kk}'
    fn_psf = f'{ab.clusters[l]}_{ab.bands[m]}_{kk}.psf'
    fn_fits = f'{orig_dir}{fn}_nan.fits' if survey == 'stack' else f'{orig_dir}{fn}.fits'
    fn_check = f'{orig_dir}{fn}_check_nan.fits' if survey == 'stack' else f'{orig_dir}{fn}_check.fits'
    with fits.open(fn_fits) as hdu_orig, \
            fits.open(fn_check) as hdu_ch, \
            open(f'{mock_dir}{fn}_mock_in{ext}.cat', 'w') as fh:

        fh.writelines('# x y mag \n')
        img = hdu_orig[0].data
        x, y = np.meshgrid(np.arange(hdu_orig[0].shape[1]), np.arange(hdu_orig[0].shape[0]))
        for i in np.arange(mag_ran[0] + offset, mag_ran[1] + offset, mag_in_step):
            rand_x = random.randint(50, hdu_orig[0].shape[1] - 50)
            rand_y = random.randint(50, hdu_orig[0].shape[0] - 50)
            while hdu_ch[0].data[rand_y, rand_x] != 0 or np.isnan(hdu_orig[0].data[rand_y, rand_x]):
                rand_x = random.randint(50, hdu_orig[0].shape[1] - 50)
                rand_y = random.randint(50, hdu_orig[0].shape[0] - 50)
            flux = 10 ** ((i - zp) / -2.5)
            mod = Gaussian2D(amplitude=flux, x_mean=rand_x, y_mean=rand_y, x_stddev=2, y_stddev=2)
            img_mod = mod(x, y)
            mod_mag = zp - 2.5 * np.log10(np.sum(img_mod))
            img = img + img_mod
            fh.writelines(f'{rand_x + 1} {rand_y + 1} {mod_mag}\n')
            # reg.writelines(f'circle {rand_x + 1} {rand_y + 1} 5 \n')
        fits.writeto(f'{mock_dir}{fn}_mock{ext}.fits', img, hdu_orig[0].header, overwrite=True)

    os.system(f"sex {mock_dir}{fn}_mock{ext}.fits -PSF_NAME {orig_dir}{fn_psf} "
              f"-CATALOG_NAME {fn}_mock_out{ext}.cat -MAG_ZEROPOINT {zp} -VERBOSE_TYPE QUIET ")

    mock_in = ascii.read(f'{mock_dir}{fn}_mock_in{ext}.cat')
    mock_out = ascii.read(f'{sex_dir}{fn}_mock_out{ext}.cat')

    matched_in = np.ones(len(mock_in), dtype=bool)

    for i in range(0, len(matched_in)):
        for j in range(0, len(mock_out)):
            matched_in[i] = (abs(mock_in['x'][i] - mock_out['X_IMAGE'][j]) < 2) & \
                            (abs(mock_in['y'][i] - mock_out['Y_IMAGE'][j]) < 2)
            if matched_in[i]:
                print(f"{mock_in['x'][i]} {mock_out['X_IMAGE'][j]} {mock_in['y'][i]} {mock_out['Y_IMAGE'][j]} "
                      f"{mock_in['mag'][i]} {mock_out['MAG_BEST'][j]}")
                break

    mags = []
    complete = []

    for i in np.arange(mag_ran[0], mag_ran[1], mag_step):
        mags.append(i)
        num_in = np.sum((mock_in['mag'] > (i - mag_step * .5)) & (mock_in['mag'] < (i + mag_step * .5)))
        num_out = np.sum((mock_in['mag'][matched_in] > (i - mag_step * .5)) &
                         (mock_in['mag'][matched_in] < (i + mag_step * .5)))
        complete.append(num_out / num_in)
        print(f'{i} {num_in} {num_out}')

    fig = plt.figure()
    plt.plot(mags, complete)
    plt.savefig(f'{sex_dir}{ab.clusters[l]}_{ab.bands[m]}_{survey}_{kk}{ext}.png')
    plt.close(fig)

    with open(f'{sex_dir}{ab.clusters[l]}_{ab.bands[m]}_{kk}{ext}.npy', 'wb') as f:
        np.save(f, np.array(mags))
        np.save(f, np.array(complete))

if ANAL_TOT:
    # fig = plt.figure(figsize=(10, 12))
    # gs = gridspec.GridSpec(len(ab.clusters), len(ab.bands), hspace=0, wspace=0)
    fig, axs = plt.subplots(ncols = len(ab.bands), nrows= len(ab.clusters), figsize=(10, 12))
    fig1 = plt.figure(figsize=(10, 3))
    gs = fig1.add_gridspec(1, 3, hspace=0, wspace=0)
    axs1 = gs.subplots(sharey='row')


for m in range(0, len(ab.bands)):
    comb90 = []
    comb01 = []
    sing90 = []
    sing01 = []
    sdss90 = []
    sdss01 = []
    deca90 = []
    deca01 = []
    for l in range(0, len(ab.clusters)):
        for survey in surveys[l]:       # stack n = 0 or single n = 1 or both n = 0, 1 or SDSS n = 2
            #mm.print_time()
            print(f'{ab.clusters[l]} {ab.bands[m]} {survey} MOSAIC start ')
            if survey == 'short':
                orig_dir = f'{ab.work_dir}fits/best_single_extracted/'
                mock_dir = f'{ab.work_dir}fits/best_single_extracted_mock/'
                zp = 25 + ab.short_a[l][m] + ab.short_b[l][m]
                mag_in_step = 0.05
                mag_step = 0.5
                hdu = fits.open(f'/Users/duhokim/work/abell/fits/best_single/{ab.short_sci_fn[l][m]}')
                ntile = len(hdu)
            elif survey == 'stack':
                orig_dir = f'{ab.work_dir}fits/extracted/'
                mock_dir = f'{ab.work_dir}fits/extracted_mock/'
                zp = 25 + ab.stack_a[l][m]
                mag_in_step = 0.05
                mag_step = 0.5
                hdu = fits.open(f'/Users/duhokim/work/abell/fits/stacked/{ab.clusters[l]}_{ab.bands[m]}si.fits')
                ntile = len(hdu)
            elif survey == 'SDSS':
                orig_dir = f'{ab.work_dir}fits/SDSS_extracted/'
                mock_dir = f'{ab.work_dir}fits/SDSS_mock/'
                zp = 22.5
                mag_in_step = 0.05
                mag_step = 0.5
                ntile = 26  # only 1 tile
            else:   # DECaLS
                if m == 0:  # no u band in DECaLS
                    continue
                orig_dir = f'{ab.work_dir}fits/DECaLS/'
                mock_dir = f'{ab.work_dir}fits/DECaLS_mock/'
                zp = 22.5
                mag_in_step = 0.05
                mag_step = 0.5
                ntile = 26  # only 1 tile

            sex_dir = f'{ab.work_dir}sex/run_{survey}_mock/'
            os.chdir(sex_dir)

            if ANAL_TOT:
                mags_all = [[] for x in range(ntile - 1)]
                comp_all = [[] for x in range(ntile - 1)]

                for k in range(1, ntile):
                    with open(f'{sex_dir}{ab.clusters[l]}_{ab.bands[m]}_{k}{ext}.npy', 'rb') as f:
                        mags = np.load(f)
                        comp = np.load(f)
                        mags_all[k - 1] = mags
                        comp_all[k - 1] = comp

                ml90, ml1 = find_90_1(mags_all[0], np.array(comp_all).mean(axis=0))
                print(f'{ml90} {ml1}')

                if survey == 'stack':
                    col = 'green'
                    line1, = axs[l,m].plot(mags_all[1], np.array(comp_all).mean(axis=0), color=col, label='Combined')
                    comb90.append(ml90)
                    comb01.append(ml1)
                elif survey == 'short':
                    col = 'blue'
                    line2, = axs[l, m].plot(mags_all[1], np.array(comp_all).mean(axis=0), '--', color=col, label='Single')
                    sing90.append(ml90)
                    sing01.append(ml1)
                elif survey == 'SDSS':
                    col = 'grey'
                    line3, = axs[l, m].plot(mags_all[1], np.array(comp_all).mean(axis=0), ':', color=col, label='SDSS')
                    sdss90.append(ml90)
                    sdss01.append(ml1)
                else:   # DECaLS
                    col = 'gold'
                    line4, = axs[l, m].plot(mags_all[1], np.array(comp_all).mean(axis=0), '-.', color=col, label='DECaLS')
                    deca90.append(ml90)
                    deca01.append(ml1)


                axs[l, m].fill_between(mags_all[1],
                                np.array(comp_all).mean(axis=0)+np.array(comp_all).std(axis=0),
                                np.array(comp_all).mean(axis=0)-np.array(comp_all).std(axis=0),
                                alpha=0.1,
                                color=col
                                )
                if l == 6:
                    axs[l,m].set_xlabel(rf'$m_{ab.bands[m]}$', fontsize=fs)
                if l < 6:
                    plt.setp(axs[l,m].get_xticklabels(), visible=False)
                if l == 3 and m == 0:
                    axs[l,m].set_ylabel('Completeness', fontsize=fs)
                if m == 2:
                    axs[l,m].yaxis.set_label_position('right')
                    axs[l,m].set_ylabel(f'{ab.clusters[l]}', fontsize=fs, rotation=270, labelpad=20)
                if m > 0:
                    plt.setp(axs[l,m].get_yticklabels(), visible=False)
                # if l == 0 and m == 0:

                axs[l,m].tick_params(direction='in', top=True, right=True)
            else:
                with Pool(NUM_CPUS) as p:
                    p.map(mock_worker, range(1, ntile))

    fill1 = axs1[m].fill_between(range(len(sdss90)), sdss90, sdss01, alpha=0.5, color='grey', label='SDSS')
    fill2 = axs1[m].fill_between(range(len(deca90)), deca90, deca01, alpha=0.5, color='darkgoldenrod', label='DECaLS')
    fill3 = axs1[m].fill_between(range(len(comb90)), [comb90[i] for i in [1, 2, 6, 0, 3, 4, 5]],
                                 [comb01[i] for i in [1, 2, 6, 0, 3, 4, 5]], alpha=0.5, color='green', label='Our combined')
    # fill4 = axs1[m].fill_between(range(len(sing90)), [sing90[i] for i in [1, 2, 6, 0, 3, 4, 5]],
    #                              [sing01[i] for i in [1, 2, 6, 0, 3, 4, 5]], alpha=0.5, color='blue', label='Single')

    axs1[m].tick_params(direction='in', top=True, right=True)
    axs1[m].set_xticks(range(len(comb90)), ['A2399', 'A2670', 'A3716', 'A754', 'A3558', 'A3574', 'A3659'],
                    rotation=70, size=11)

    axs1[m].text(0.8, 0.1, r'$' + ab.bands[m] + '\prime$', transform=axs1[m].transAxes, size=30)
    axs1[m].set_ylim([28, 20])

    box = axs1[m].get_position()
    box.y0 = box.y0 + 0.08
    box.y1 = box.y1 + 0.08
    box.x0 = box.x0 + 0.05
    box.x1 = box.x1 + 0.05
    axs1[m].set_position(box)

if ANAL_TOT:
    axs[0, 1].legend(handles=[line1, line2, line3, line4], loc=3)
    fig.savefig(f'{ab.work_dir}plot/recovered_all{ext}_new_a1.png')

    axs1[0].set_ylabel(r'$m_{lim,90\%}$'+'\n'+r'$\sim$'+'\n'+r'$m_{lim,1\%}$', size=fs, rotation=0, ha='center')
    # axs1[0].yaxis.set_label_position('left')
    axs1[0].yaxis.set_label_coords(-0.3, 0.3)
    # axs1[0].legend(bbox_to_anchor=(2, 24, 5, 3),handles=[fill3, fill4, fill1, fill2], frameon=False)
    # axs1[0].legend(loc='best', handles=[fill3, fill4], frameon=False)
    # axs1[1].legend(loc='best', handles=[fill1, fill2], frameon=False)
    axs1[0].legend(loc='best', handles=[fill1, fill2, fill3], frameon=False)
    fig1.savefig(f'{ab.work_dir}plot/recovered_all{ext}_90_01_one_excl_single.png')

os.chdir(cur_dir)