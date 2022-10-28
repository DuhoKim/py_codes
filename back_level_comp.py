from astropy.io import ascii
import matplotlib.pyplot as plt

tab = ascii.read('/Users/duhokim/work/abell/galapagos/A2670_test/A2670_back_values.txt')

fig, axs = plt.subplots(1, figsize=(8, 8))

axs.scatter(tab['galapagos_sky'], tab['nanmedian'], label='median', alpha=0.3, marker='o')
axs.scatter(tab['galapagos_sky'], tab['sex'], label='SExtractor', alpha=0.3, marker='^')

axs.plot([1215, 1230], [1215, 1230], linestyle=':', color='black')

axs.set_xlim([1215, 1230])
axs.set_ylim([1215, 1230])

axs.tick_params(direction='in', top=True, right=True, labelsize=15)

axs.set_xlabel('Background level (GALAPAGOS)', fontsize=20)
axs.set_ylabel('Background level (Median, SExtractor)', fontsize=20)

axs.legend(fontsize=20)

fig.savefig('/Users/duhokim/work/abell/plot/back.png')
plt.close(fig)

