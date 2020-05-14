from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt

cluster_name = "A2670"

data_dir=("/Users/dkim108/Documents/work/sex/")
plot_dir=("/Users/dkim108/Documents/work/plot/")

sex_result = ascii.read(data_dir+"a2670_r_nthresh64_back64.cat")
#sex_result_good = (sex_result['FLAGS'] < 4) & (sex_result['MAG_AUTO'] < 90)
sex_result_good = sex_result['MAG_AUTO'] < 90

eff_radius = np.sqrt(sex_result['A_WORLD'] * sex_result['B_WORLD']) * 3600

plt.figure()
plt.title('A2670 r band size vs mag')
plt.scatter(eff_radius[sex_result_good], sex_result['MAG_AUTO'][sex_result_good], alpha=0.1, s=1)
plt.xlabel('effective radius [arcsec]')
plt.ylabel('MAG_AUTO')
plt.gca().invert_yaxis()
#plt.savefig(plot_dir+cluster_name+'_rad_vs_mag.pdf')

plt.show()