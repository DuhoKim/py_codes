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
import xlsxwriter
import statistics
from statistics import mode
from collections import Counter

def most_frequent(List):
    occurence_count = Counter(List)
    return occurence_count.most_common(1)[0][0], occurence_count.most_common(1)[0][1]/len(List)*100

users = ['AL', 'AR', 'FP', 'JC', 'YJ']
morph = ['E', 'S0', 'Spr', 'Irr', 'N/A']
feat = ['none', 'asym', 'warp', 'fan', 'int', 'bridge', 'sheet', 'jelly']

workbook = xlsxwriter.Workbook(ab.work_dir+'vis/results/A3558_all.xlsx')

header_format = workbook.add_format({
    'border': 1,
    'bg_color': '#C6EFCE',
    'bold': True,
    'text_wrap': True,
    'valign': 'vcenter',
    'indent': 1,
})
header_format.set_font_size(30)

id_format = workbook.add_format({
    'bold': True,
    'valign': 'vcenter',
})
id_format.set_font_size(30)

ws = workbook.add_worksheet('A3558 Spec Gal Vis Ins Res')

# ws.set_column(0, 0, 20)  # Set 'ID' column's width [# of character]
# ws.set_column(1, 8, 40)  # Set 'A, B, C' column's width [# of character]
# ws.set_column(6, 6, 140)
# ws.set_column(7, 8, 40, id_format)
# ws.set_column(9, 9, 250, id_format)
# ws.set_default_row(240)  # Set the default row height

ws.write(0, 0, 'ID')

ws.write(0, 2, 'Morph_A')
ws.write(0, 3, 'Morph_B')
ws.write(0, 4, 'Morph_C')
ws.write(0, 5, 'Morph_D')
ws.write(0, 6, 'Morph_E')

ws.write(0, 8, 'Morph_freq')
ws.write(0, 9, 'Morph_%')

ws.write(0, 11, 'Feat_A')
ws.write(0, 12, 'Feat_B')
ws.write(0, 13, 'Feat_C')
ws.write(0, 14, 'Feat_D')
ws.write(0, 15, 'Feat_E')

ws.write(0, 17, 'Feat_freq')
ws.write(0, 18, 'Feat_%')

id = []
row_num = []
id_all = [[], [], [], [], []]
m_all = [[],[],[],[],[]]
f_all = [[],[],[],[],[]]

for i in range(0, len(users)):
    cur_row = 0
    with open(ab.work_dir + f'vis/results/A3558_{users[i]}.vis', 'r') as vis:
        for line in vis:
            item = line.split(' ')
            if len(item) > 2:
                id_all[i].append(item[0])
                m_all[i].append(int(item[1]))
                f_all[i].append(int(item[2]))
                if item[0] in id:
                    id_ind = id.index(item[0])
                    ws.write(row_num[id_ind], i + 2, morph[int(item[1])])
                    ws.write(row_num[id_ind], i + 11, feat[int(item[2])])
                else:
                    id.append(item[0])
                    row_num.append(len(id))
                    ws.write(row_num[-1], 0, id[-1])
                    ws.write(row_num[-1], i + 2, morph[int(item[1])])
                    ws.write(row_num[-1], i + 11, feat[int(item[2])])


for i in range(0, len(id)):
    this_m = []
    this_f = []
    for j in range(0, 5):
        if id[i] in id_all[j]:
            id_ind = id_all[j].index(id[i])
            this_m.append(m_all[j][id_ind])
            this_f.append(f_all[j][id_ind])

    most_m, perc_m = most_frequent(this_m)
    most_f, perc_f = most_frequent(this_f)

    ws.write(row_num[i], 8, f'{morph[most_m]}')
    ws.write(row_num[i], 9, int(perc_m))
    ws.write(row_num[i], 17, f'{feat[most_f]}')
    ws.write(row_num[i], 18, int(perc_f))

workbook.close()