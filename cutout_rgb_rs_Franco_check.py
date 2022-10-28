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
import calc_kcor as kcor


save_dir = ab.work_dir + f'pics/A2670_Franco/'

fig = plt.figure()

r_vis = []
g_vis = []
g_r_vis = []

with open(save_dir+f'A2670_rs_or_spec_tmp.vis', 'rb') as vis:
    cat = ascii.read(save_dir+f'A2670.cat')
    for line in vis:
        ind = np.where(cat['NUMBER'] == int(line))[0][0]
        gmag = cat['MAG_AUTO_g'][ind]
        rmag = cat['MAG_AUTO_r'][ind]

        gmag -= kcor.calc_kcor('g', ab.redshifts[2], 'g - r',
                                                     gmag - rmag)
        rmag -= kcor.calc_kcor('r', ab.redshifts[2], 'g - r',
                                                     gmag - rmag)

        r_vis.append(rmag)
        g_r_vis.append(gmag - rmag)


        # r_vis.append(cat['MAG_AUTO_r'][ind[0][0]])
        # g_r_vis.append(cat['MAG_AUTO_g'][ind[0][0]] - cat['MAG_AUTO_r'][ind[0][0]])


plt.scatter(r_vis,
         g_r_vis)
plt.ylim([-3, 3])
plt.xlim([12, 20])
plt.show()
