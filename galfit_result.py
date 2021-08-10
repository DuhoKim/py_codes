from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table, vstack
import numpy as np
import os
import shutil
import glob
from pandas import DataFrame, read_csv
import pandas as pd
from astropy.cosmology import Planck15 as Cosmo
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.backends.backend_pdf as pdf
from matplotlib.colors import LogNorm
import math
from astropy import constants as const

work_dir = '/Users/duhokim/work/abell/galfit/A3558_no_weight/'
cat = ascii.read("/Users/duhokim/work/abell/spec/Shapley/catalog.dat")

c = 299792.5 # speed of light in km/s

fn = sorted(glob.glob(work_dir + 'galfit.*'), key=os.path.getmtime, reverse=True)

btot = np.zeros(len(fn))  # bulge-to-total light ratios
btot[:] = np.nan  # Init w/ nan

btot_disk = np.zeros(len(fn))  # bulge-to-total light ratios
btot_disk[:] = np.nan  # Init w/ nan
btot_bulge = np.zeros(len(fn))  # bulge-to-total light ratios
btot_bulge[:] = np.nan  # Init w/ nan

mtot_disk = np.zeros(len(fn))  # bulge-to-total light ratios
mtot_disk[:] = np.nan  # Init w/ nan
mtot_bulge = np.zeros(len(fn))  # bulge-to-total light ratios
mtot_bulge[:] = np.nan  # Init w/ nan
mtot_disk_err = np.zeros(len(fn))  # bulge-to-total light ratios
mtot_disk_err[:] = np.nan  # Init w/ nan
mtot_bulge_err = np.zeros(len(fn))  # bulge-to-total light ratios
mtot_bulge_err[:] = np.nan  # Init w/ nan

mtot = np.zeros(len(fn))  # Total magnitude
mtot[:] = np.nan  # Init w/ nan
mtot_err = np.zeros(len(fn))  # Total magnitude
mtot_err[:] = np.nan  # Init w/ nan

ttype = np.zeros(len(fn))  # RC3 T-type
ttype[:] = np.nan  # Init w/ nan
ttype_err = np.zeros(len(fn))  # RC3 T-type error
ttype_err[:] = np.nan  # Init w/ nan

chi = np.zeros(len(fn))  # Reduced Chi squre
chi[:] = np.nan  # Init w/ nan
oned = np.zeros(len(fn))
oned[:] = np.nan
names = ["" for i in range(len(fn))]

# subimg = pdf.PdfPages(work_dir+'subs_3rd.pdf')
btot_fn = open(work_dir + 'btot.txt', 'w')
fit_log = open(work_dir + 'fit.log', 'r')
fit_log_lines = fit_log.readlines()

skip = False  # flag for check duplication

for i in range(0, len(fn)):
    res = open(fn[i])
    comp = 0
    disk = np.nan
    bulge = np.nan
    sky = 1050

    for line in res:
        line_spl = line.split(' ')
        if skip:
            break
        if '0) expdisk' in line:
            comp = 1
            continue
        if '0) devauc' in line:
            comp = 2
            continue

        for j in range(0, len(line_spl)):
            if 'feedme' in line_spl[j]:
                name_spl = line_spl[j].split('.')
                name = name_spl[0]

                same_name_idx = np.isin(names, name)
                if sum(same_name_idx):
                    skip = True
                    break
                else:
                    names[i] = name
                z = float(cat['col18'][int(name)-1]) / c

            if 'Chi^2/nu' in line_spl[j]:
                chi[i] = float(line_spl[4][:-1])

            if 'Integrated' in line_spl[j] and comp == 1:
                disk = float(line_spl[2])
                btot_disk[i] = disk
                continue

            if 'Integrated' in line_spl[j] and comp == 2:
                bulge = float(line_spl[2])
                btot_bulge[i] = bulge
                continue

            if 'Sky background' in line_spl[j]:
                sky = float(line_spl[2])
                continue

    if skip:
        skip = False
        continue

    for j in range(len(fit_log_lines) - 1, 0, -1):
        if name+'.feedme' in fit_log_lines[j]:
            for k in range(0, 10):
                if 'expdisk' in fit_log_lines[j + k]:
                    try:
                        mtot_disk[i] = float(fit_log_lines[j + k].split(')')[1].split()[0])
                        mtot_disk_err[i] = float(fit_log_lines[j + k + 1].split(')')[1].split()[0])
                    except ValueError:
                        mtot_disk[i] = float(fit_log_lines[j + k].split(')')[1].split()[0][1:-1])
                        mtot_disk_err[i] = float(fit_log_lines[j + k + 1].split(')')[1].split()[0][1:-1])
                if 'devauc' in fit_log_lines[j + k]:
                    try:
                        mtot_bulge[i] = float(fit_log_lines[j + k].split(')')[1].split()[0])
                        mtot_bulge_err[i] = float(fit_log_lines[j + k + 1].split(')')[1].split()[0])
                    except ValueError:
                        mtot_bulge[i] = float(fit_log_lines[j + k].split(')')[1].split()[0][1:-1])
                        mtot_bulge_err[i] = float(fit_log_lines[j + k + 1].split(')')[1].split()[0][1:-1])
            break

    if np.isnan(bulge) and ~np.isnan(disk):
        btot[i] = 0
        mtot_apparent = disk
        mtot_bulge_err[i] = 99
    else:
        if disk > bulge:
            btot[i] = (10 ** ((disk - bulge) / 2.5) / (10 ** ((disk - bulge) / 2.5) + 1.0))
        if disk == bulge:
            btot[i] = 0.5
        if disk < bulge:
            btot[i] = (1.0 / (10 ** ((bulge - disk) / 2.5) + 1.0))
        mtot_apparent = 31.456 - 2.5 * np.log10(10 ** ((31.456 - disk) / 2.5) + 10 ** ((31.456 - bulge) / 2.5))

    btot_fn.writelines(
        name + ' ' + str(btot[i]) + ' ' + str(chi[i]) + ' ' + str(mtot_apparent) + ' \n')

    ########## ADD SUBTRACT PLOT in PDF ###########
    hdu = fits.open(work_dir + name + '.out.fits')
    img_orig = hdu[1].data
    img_model = hdu[2].data
    img_sub = hdu[3].data

    fig = plt.figure(figsize=(12, 3))

    h1 = fig.add_axes([0, 0, 0.25, 1])
    h2 = fig.add_axes([0.25, 0, 0.25, 1])
    h3 = fig.add_axes([0.5, 0, 0.25, 1])
    h4 = fig.add_axes([0.75, 0, 0.25, 1])

    h1.imshow(img_orig - sky, norm=LogNorm(vmin=0.1, vmax=100), origin='lower')
    h2.imshow(img_orig - sky, norm=LogNorm(vmin=0.1, vmax=1000), origin='lower')
    h3.imshow(img_model - sky, norm=LogNorm(vmin=0.1, vmax=1000), origin='lower')
    h4.imshow(img_sub, norm=LogNorm(vmin=0.1, vmax=100), origin='lower')

    kpc_arcmin = Cosmo.kpc_proper_per_arcmin(z)
    arcmin_5kpc = 5.0 / kpc_arcmin
    frac_5kpc = arcmin_5kpc * 100.0 / img_orig.shape[0]

    props = dict(boxstyle='round', facecolor='wheat')
    h1.text(0.03, 0.92, r'$r \prime$', color='g', size=14, transform=h1.transAxes, bbox=props)
    h2.text(0.03, 0.92, r'$r \prime$', color='g', size=14, transform=h2.transAxes, bbox=props)
    h3.text(0.03, 0.92, 'Model with B/T=' + str('{:5.2}'.format(btot[i])), color='g', size=14, transform=h3.transAxes, bbox=props)
    h4.text(0.03, 0.92, r'Residual ($\chi^{2}$=' + str(chi[i])+')', color='g', size=14, transform=h4.transAxes, bbox=props)
    #h1.plot([0.05, 0.05 + frac_5kpc.value], [0.02, 0.02], color='darkslategray', transform=h1.transAxes, linewidth=5)
    #h1.text(0.05, 0.05, '5kpc, '+'{:4.2f}'.format(float(arcmin_5kpc.value)) + '\'', color='g', fontsize=24,
    #        transform=h1.transAxes, bbox=props)

    h1.get_xaxis().set_visible(False)
    h1.get_yaxis().set_visible(False)
    h2.get_xaxis().set_visible(False)
    h2.get_yaxis().set_visible(False)
    h3.get_xaxis().set_visible(False)
    h3.get_yaxis().set_visible(False)
    h4.get_xaxis().set_visible(False)
    h4.get_yaxis().set_visible(False)

    fig.savefig(work_dir + 'png/' + name + '.png')
    plt.close(fig)

btot_fn.close()
# subimg.close()
