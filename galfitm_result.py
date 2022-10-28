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
import my_module as mm
import math
import copy
import img_scale
import abell_cluster_module as ab

work_dir = '/Users/duhokim/work/abell/galapagos/A3558_rsr_all/t1/galfit/'
cat = ascii.read("/Users/duhokim/work/abell/spec/Shapley/catalog.dat")

n = 7
# cmap = 'PuOr'
cmap = 'viridis'

kron_fac = 2.2
rot_box = 1.5

# for k in range(0, len(ab.clusters)):
for k in range(0, 3):

    ntile = 9 if k > 0 else 8
    # ntile = 1

    for j in range(1, ntile+1):

        fn = sorted(glob.glob(f'/Users/duhokim/work/abell/galapagos/{ab.clusters[k]}_r/t{j}/galfit/*_gf.fits'), key=os.path.getmtime, reverse=True)

        for i in range(0, len(fn)):
            fits_name = fn[i].split('/')[-1]
            name = fits_name.split('_')[0]
            hdu = fits.open(fn[i])
            img_orig = hdu[0].data
            img_model = hdu[1].data
            img_sub = hdu[2].data

            fig = plt.figure(figsize=(12, 4))

            h1 = fig.add_axes([0, 0, 0.33, 1])
            h2 = fig.add_axes([0.33, 0, 0.33, 1])
            h3 = fig.add_axes([0.66, 0, 0.33, 1])

            sky = np.nanmedian(img_orig)
            ser_idx = -1

            res = open(fn[i][:-4] + 'galfit.01')

            for line in res:
                if 'Sky background' in line:
                    sky = float(line.split(' ')[2])
                if 'Position x' in line:
                    pos_x = float(line.split(' ')[2])-1
                if 'Position y' in line:
                    pos_y = float(line.split(' ')[2])-1
                if 'Integrated magnitude' in line:
                    mag = float(line.split(' ')[2])
                if 'R_e' in line:
                    reff = float(line.split(' ')[2])
                if 'Sersic index' in line:
                    ser_idx = float(line.split(' ')[2])
                if 'Axis ratio' in line:
                    boa = float(line.split(' ')[2])
                if 'Position angle' in line:
                    pa = float(line.split(' ')[1])
                if 'Component number: 3' in line:
                    break

            h1.imshow(img_orig - sky, norm=LogNorm(vmin=0.1, vmax=10000), origin='lower', cmap=cmap)

            h2.imshow(img_model - sky, norm=LogNorm(vmin=0.1, vmax=1000), origin='lower', cmap=cmap)
            # h2.text(0.1, 0.9, f'Sersic index = {ser_idx}', color='white', transform=h2.transAxes)

            h3.imshow(img_sub, norm=LogNorm(vmin=0.1, vmax=100), origin='lower', cmap=cmap)

            h1.get_xaxis().set_visible(False)
            h1.get_yaxis().set_visible(False)
            h2.get_xaxis().set_visible(False)
            h2.get_yaxis().set_visible(False)
            h3.get_xaxis().set_visible(False)
            h3.get_yaxis().set_visible(False)

            results_dir = f'/Users/duhokim/work/abell/galapagos/{ab.clusters[k]}_r/t{j}/galfit/png/'

            if not os.path.isdir(results_dir):
                os.makedirs(results_dir)

            fig.savefig(results_dir + name + f'_gala_vir_10000.png')
            plt.close(fig)
