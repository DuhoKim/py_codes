import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import abell_cluster_module as ab
import importlib
importlib.reload(ab)
from scipy.stats import pearsonr

colors = ['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'cyan', 'magenta', 'yellow', 'salmon','skyblue',
          'lightgrey', 'darkgrey', 'black']


# feature type fraction relation with the dynamical parameters for galaxy clusters
clusters =               ['A754',   'A2399',    'A2670', 'A3558','A3562','A3716']
parameters = ['c', 'w', 'p', 'k', 'a', 'd']
xlabels= [r'$\log_{10}(c)$', r'$\log_{10}(w)$', r'$\log_{10}(P_3/P_0)$', r'$\kappa$', r'$\log_{10}(a)$', r'$\delta$']
yuan = np.array([
    [np.nan,   np.nan,     -0.48,  -0.59,  -0.49,  -1.26],
    [np.nan,   np.nan,     -1.92,  -1.93,  -2.29,  -2.83],
    [np.nan,   np.nan,     -7.04,  -6.96,  -6.23,  -4.74],
    [1.7,      np.nan,     1.11,   1.51,   1.42,   1.38],
    [-1.78,    np.nan,     -1.18,  -1.2,   -1.36,  -1.03],
    [0.24,     np.nan,     0.22,   0.5,    0.32,   0.52]])

yuan_err = np.array([
    [np.nan,   np.nan,     0.01,  0.01,  0.01,  0.01],
    [np.nan,   np.nan,     0.01,  0.01,  0.01,  0.01],
    [np.nan,   np.nan,     0.01,  0.01,  0.01,  0.01],
    [0.0,      np.nan,     0.0,   0.0,   0.0,   0.0],
    [0.01,    np.nan,     0.01,   0.01,   0.01,  0.01],
    [0.01,     np.nan,     0.01,   0.01,   0.01,   0.01]])

invert_x = [True, False, False, False, False, False]

tot =   [126,   115,    171,    169,    132,    59]
m   =   [3,     1,      8,      13,     7,      0]
pm   =  [11,    6,      11,     19,     6,      2]
m3   =   [3,     1,      6,      4,     2,      0]
pm3   =  [5,    3,      2,     6,     3,      0]


m_frac = np.zeros((2, len(tot)))
pm_frac = np.zeros((2, len(tot)))
m3_frac = np.zeros((2, len(tot)))
pm3_frac = np.zeros((2, len(tot)))
mpm_frac = np.zeros((2, len(tot)))

for i in range(len(tot)):
    m_frac[0][i] = m[i] / tot[i] * 100
    m_frac[1][i] = m[i] ** 0.5 / tot[i] * 100
    m3_frac[0][i] = m3[i] / tot[i] * 100
    m3_frac[1][i] = m3[i] ** 0.5 / tot[i] * 100
    pm_frac[0][i] = pm[i] / tot[i] * 100
    pm_frac[1][i] = pm[i] ** 0.5 / tot[i] * 100
    pm3_frac[0][i] = pm3[i] / tot[i] * 100
    pm3_frac[1][i] = pm3[i] ** 0.5 / tot[i] * 100
    mpm_frac[0][i] = (m[i] + pm[i]) / tot[i] * 100
    mpm_frac[1][i] = (m[i] + pm[i]) ** 0.5 / tot[i] * 100

fig = plt.figure(figsize=(8, 8))
gs = fig.add_gridspec(2, 3, wspace=0, hspace=0.3)

for i in range(yuan.shape[0]):
    axs = fig.add_subplot(gs[int(i/3), i%3])
    if i % 3 == 0:
        axs.set_ylabel('%', fontsize=15)
    else:
        axs.set_yticks([])
    if i == 0:
        axs.text(0.05, 1.05, 'Relaxed', transform=axs.transAxes, fontsize=15)
    elif i == 2:
        axs.text(0.5, 1.05, 'Disturbed', transform=axs.transAxes, fontsize=15)
    nn = np.where(~np.isnan(yuan[i]))[0]    # index for not nan
    # nn2 = nn[:-1]                           # index for not nan excluding A3716
    axs.errorbar(yuan[i][nn], m_frac[0][nn], xerr=yuan_err[i][nn], yerr=m_frac[1][nn], fmt='*')
    axs.errorbar(yuan[i][nn], pm_frac[0][nn], xerr=yuan_err[i][nn], yerr=pm_frac[1][nn], fmt='P')
    axs.errorbar(yuan[i][nn], m3_frac[0][nn], xerr=yuan_err[i][nn], yerr=m3_frac[1][nn], fmt='*')
    axs.errorbar(yuan[i][nn], pm3_frac[0][nn], xerr=yuan_err[i][nn], yerr=pm3_frac[1][nn], fmt='P')
    axs.set_xlabel(xlabels[i], fontsize=15)

    corr_mat = np.corrcoef(yuan[i][nn], m_frac[0][nn])
    if invert_x[i]:
        corr_xy = corr_mat[0, 1] * -1
    else:
        corr_xy = corr_mat[0, 1]
    axs.text(0.1, 0.9, rf'P_M={"{:.2}".format(corr_xy)}', transform=axs.transAxes, c='blue', fontsize=12)
    corr_mat = np.corrcoef(yuan[i][nn], pm_frac[0][nn])
    if invert_x[i]:
        corr_xy = corr_mat[0, 1] * -1
    else:
        corr_xy = corr_mat[0, 1]
    axs.text(0.1, 0.8, rf'P_PM={"{:.2}".format(corr_xy)}', transform=axs.transAxes, c='orange', fontsize=12)
    corr_mat = np.corrcoef(yuan[i][nn], m3_frac[0][nn])
    if invert_x[i]:
        corr_xy = corr_mat[0, 1] * -1
    else:
        corr_xy = corr_mat[0, 1]
    axs.text(0.1, 0.7, rf'P_M3={"{:.2}".format(corr_xy)}', transform=axs.transAxes, c='green', fontsize=12)
    corr_mat = np.corrcoef(yuan[i][nn], pm3_frac[0][nn])
    if invert_x[i]:
        corr_xy = corr_mat[0, 1] * -1
    else:
        corr_xy = corr_mat[0, 1]
    axs.text(0.1, 0.6, rf'P_PM3={"{:.2}".format(corr_xy)}', transform=axs.transAxes, c='red', fontsize=12)
    axs.tick_params(direction='in', top=True, right=True, labelsize=15)
    if i == 5:
        axs.set_xticks([0.2, 0.4, 0.6])
    elif i == 4:
        axs.set_xticks([-1.8, -1.5, -1.2])
    if invert_x[i]:
        axs.invert_xaxis()

fig.savefig(ab.plot_dir + f'dyn_vs_frac.png')

