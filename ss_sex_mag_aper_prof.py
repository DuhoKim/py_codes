from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
import Funcs



plot_dir=("/Users/dkim108/Documents/work/plot/")
cat_dir=("/Users/dkim108/Documents/work/cat/SouthernStandardStars/")

ss_date = '19aug'
band = 'g'
rad_in_arcsec = np.multiply( [4, 8, 12, 20, 30, 40, 60, 100], 0.263)

if ss_date == '19aug':
    sex_dir=("/Users/dkim108/Documents/work/sex/19aug/")
    ss_name = ['220100', 'E8-A', 'LSE_259', '190000']
    ss_legend_name = ['220100 (1.96)', 'E8-A (1.21)', 'LSE_259 (1.75)', '190000 (1.23)']
    ss_cat_name = ['220100-300000', 'E8_a', 'LSE_259', '190000-295600']
    ss_airmass = [1.96, 1.21, 1.75, 1.23]
    ss_col = ['indianred', 'tan', 'olivedrab', 'mediumpurple']

    if band == 'u':     # col4: u_mag, col7: g_mag, col10: r_mag
        mag_col = 'col4'
        err_col = 'col5'
        a_best = 3.349
        b_best = -0.318
        max_sep = 1.0
        class_star_lim = 0.5
        flag_lim = 5
    elif band == 'g':
        mag_col = 'col7'
        err_col = 'col8'
        a_best = 2.857
        b_best = -0.183
        max_sep = 1.0
        class_star_lim = 0.98
        flag_lim = 0
    elif band == 'r':
        mag_col = 'col10'
        err_col = 'col11'
        a_best = 3.024
        b_best = -0.117
        max_sep = 1.0
        class_star_lim = 0.98
        flag_lim = 0

elif ss_date == '20aug':
    sex_dir=("/Users/dkim108/Documents/work/sex/20aug/")
    ss_name = ['E5-A', 'E6-A', 'LSE_259_1_11', '220100_2_44', '220100_1_01', 'LSE_259_1_71']
    ss_legend_name = ['E5-A (1.7)', 'E6-A (1.12)', 'LSE_259 (1.11)', '220100 (2.44)', '220100 (1.01)', 'LSE_259 (1.71)']
    ss_cat_name = ['E5_a', 'E6_a', 'LSE_259', '220100-300000', '220100-300000', 'LSE_259']
    ss_airmass = [1.7, 1.12, 1.11, 2.44, 1.01, 1.71]

elif ss_date == '21aug':
    sex_dir=("/Users/dkim108/Documents/work/sex/21aug/")
    ss_name = ['E5-A', 'E6-A', 'LSE_259_1_11', 'LSE_259_1_8']
    ss_legend_name = ['E5-A (1.65)', 'E6-A (1.11)', 'LSE_259 (1.11)', 'LSE_259 (1.8)']
    ss_cat_name = ['E5_a', 'E6_a', 'LSE_259', 'LSE_259']
    ss_airmass = [1.65, 1.11, 1.11, 1.8]



plt.close('all')

# Just a figure and one subplot

#ax.title('SS '+ss_date+' '+band+'-band MAG_APER profiles')

sample_numbers = np.zeros(len(ss_name))
avgs = np.zeros((len(ss_name), len(rad_in_arcsec)))
stds = np.zeros((len(ss_name), len(rad_in_arcsec)))
# tot_tmp_diffs = [ [] for i in range(len(rad_in_arcsec))]
# tot_stds = np.zeros(len(rad_in_arcsec))

for i in range(0, len(ss_name)):
    sex_result = ascii.read(sex_dir+ss_name[i]+'_'+band+'a.cat')
    sex_result_good = (sex_result['MAG_ISO'] < 90) & (sex_result['FLAGS'] <= flag_lim) & \
                      (sex_result['CLASS_STAR'] > class_star_lim)
    sex_coords = SkyCoord(sex_result['ALPHA_J2000'][sex_result_good],
                          sex_result['DELTA_J2000'][sex_result_good],
                          unit='deg')

    ss_cat = ascii.read(cat_dir+ss_cat_name[i]+'.dat.trimmed', format='no_header')

    cat_coords = SkyCoord(ss_cat['col2'], ss_cat['col3'], unit=(u.hourangle, u.deg))

    idx, d2d, d3d = sex_coords.match_to_catalog_sky(cat_coords)
    sep_constraint = (d2d.arcsec < max_sep) & (ss_cat[idx][err_col] > 0)
    sex_matches = sex_result[sex_result_good][sep_constraint]
    cat_matches = ss_cat[idx[sep_constraint]]
    sample_numbers[i] = len(sex_matches)
    nrow = int(np.ceil(np.sqrt(len(sex_matches))))

    fig, axarr = plt.subplots(nrow, nrow, sharex='col', sharey='row', gridspec_kw={'hspace':0, 'wspace':0})

    tmp_diffs = np.zeros((len(sex_matches), len(rad_in_arcsec)))
    for l in range(0, len(sex_matches)):
        flux = [sex_matches['FLUX_APER'][l], sex_matches['FLUX_APER_1'][l], sex_matches['FLUX_APER_2'][l],
                sex_matches['FLUX_APER_3'][l], sex_matches['FLUX_APER_4'][l], sex_matches['FLUX_APER_5'][l],
                sex_matches['FLUX_APER_6'][l], sex_matches['FLUX_APER_7'][l]]
        mags = 25.0 - 2.5 * np.log10(flux) + a_best + b_best * ss_airmass[i]
        mag_iso = sex_matches['MAG_ISO'][l] + a_best + b_best * ss_airmass[i]

        tmp_diffs[l] = mags - cat_matches[mag_col][l]

        if l == 0:
            axarr[l/nrow, l%nrow].plot(rad_in_arcsec, mags , alpha=0.5, color=ss_col[i], linewidth = 1,
                     label=ss_legend_name[i]+' ({:d})'.format(int(sample_numbers[i])))
        else:
            axarr[l/nrow, l%nrow].plot(rad_in_arcsec, mags, alpha=0.5, color=ss_col[i], linewidth = 1)

        rad_iso = Funcs.find_x0(rad_in_arcsec, mags, mag_iso)
        axarr[l/nrow, l%nrow].plot(rad_iso, mag_iso, marker='*', alpha=0.5)
        axarr[l/nrow, l%nrow].errorbar(rad_iso, mag_iso, yerr=sex_matches['MAGERR_ISO'][l])

        rad_ss = Funcs.find_x0(rad_in_arcsec, mags, cat_matches[mag_col][l])
        axarr[l/nrow, l%nrow].plot(rad_ss, cat_matches[mag_col][l], marker='D', alpha=0.5)
        axarr[l/nrow, l%nrow].errorbar(rad_ss, cat_matches[mag_col][l], yerr=cat_matches[err_col][l])


        # axarr[l/nrow, l%nrow].plot([24, 24], [8, 20], '--', alpha=0.1)

    plt.savefig(plot_dir + ss_date + '_DECam_' + band + '_MAG_APER_prof_vs_SS_mag_'+ss_name[i]+'.pdf')
    avgs[i] = np.mean(tmp_diffs, axis=0)
    stds[i] = np.std(tmp_diffs, axis=0)

f, ax = plt.subplots()
ax.set_title('Diff b/w MAG_APER and SS trimmed catalog in '+band)
for i in range(0, len(ss_name)):
    ax.plot(rad_in_arcsec, avgs[i], linestyle=':', alpha=0.5, label=ss_legend_name[i])
    ax.errorbar(rad_in_arcsec+0.05*i, avgs[i], yerr=stds[i], alpha=0.5)

ax.plot([min(rad_in_arcsec), max(rad_in_arcsec)],[0, 0], alpha=0.5, linewidth=1)
# plt.xlim([8, 20])
# plt.ylim([8, 20])
ax.set_ylabel('Avg(MAG_APER-trimmed)')
ax.set_xlabel('diameter ["]')
# # plt.gca().invert_xaxis()
# #plt.gca().invert_yaxis()
# ax.text(25, 10, '24\" diameter for SS aper corr')
# ax.text(0.7, 0.9, r'$\star$ : MAG_ISO from SEx', transform=ax.transAxes)
# ax.text(0.7, 0.8, r'$\diamond$ : trimmed catalog', transform=ax.transAxes)
# ax.text(0.1, 0.9, band+' = '+band+'_inst (ZP=25) + {:5.3f} + {:5.3f} * X'.format(a_best, b_best), transform=ax.transAxes)
# ax.text(0.1, 0.8, '# of sample = {:4d}'.format(int(np.sum(sample_numbers))), transform=ax.transAxes)
# ax.text(0.1, 0.7, 'match < {:4.2f}\"'.format(max_sep), transform=ax.transAxes)
# ax.text(0.1, 0.6, 'FLAGS < {:2d}'.format(int(flag_lim)), transform=ax.transAxes)
# ax.text(0.1, 0.5, 'CLASS_STAR > {:4.2f}'.format(class_star_lim), transform=ax.transAxes)
ax.legend(title='Standard Star Field (Airmass)(#)', loc='upper right', prop={'size':7})
plt.savefig(plot_dir+ss_date+'_DECam_'+band+'_MAG_APER-SS_mag.pdf')
plt.close()