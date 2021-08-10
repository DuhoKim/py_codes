from astropy.io import ascii
from astropy.table import Table, vstack
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np
import os
from shutil import copyfile
from pandas import DataFrame, read_csv
import pandas as pd
from astropy.wcs import WCS

ned1=ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_1.txt")
ned2=ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_2.txt")
ned3=ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_3.txt")
ned4=ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_4.txt")
ned5=ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_5.txt")

ned_tot=vstack([ned1,ned2,ned3,ned4,ned5])

sha_cat=ascii.read('/Users/dhk/work/cat/NGC_IC/sha_ned_match_off_lt_1_SDSS_100psfs_b_over_a_gt_05_from_Tom_2_erase_Tamura_match_segmentation.csv') # Sample NGC/IC numbers

rc3=ascii.read("/Users/dhk/work/cat/NGC_IC/RC3/myrc3.dat")
#nair=ascii.read("/Users/dhk/work/cat/NGC_IC/Nair2010/table2.txt")
file=r'/Users/dhk/work/cat/NGC_IC/Steinicke/NI2018.xls'
df = pd.read_excel(file)
file=r'/Users/dhk/work/cat/NGC_IC/ttype.xls'
ttype = pd.read_excel(file)

for x in range(0,len(ned_tot)):
    if ned_tot['col2'][x] in sha_cat['col1']:
        name = ned_tot['col2'][x]
        sha_match = sha_cat[sha_cat['col1']==name]
        boa = "%.2f" % (ned_tot['col12'][x] / ned_tot['col11'][x])

        ######  WCS  read      #######
        w       = WCS('/Users/dhk/work/data/NGC_IC/SHA/SHA_NGC_IC_LONG_mean/'+name+'.fits')
        i,j     = w.all_world2pix(sha_match['col2'][0],sha_match['col3'][0],1)

        semi = ned_tot['col11'][x]                      # read semi-major axis from NED [']
        thumb = int(semi*100)                           # half size of figure size

        if name[0]=='n':
            galnum = name[3:].strip()
            ngc_match = df.loc[(df['N']=='N') & (df['NI']==int(galnum))]
        elif name[0]=='i':
            galnum = name[2:].strip()
            ngc_match = df.loc[(df['N']=='I') & (df['NI']==int(galnum))]

        pa = str(ngc_match.iloc[0,21])
        if pa=='nan':
            pa='0'

		#####  xxx.feedme    ########
        fm=open('/Users/dhk/work/data/NGC_IC/SHA/galfit/10th_const/'+name+'.const','w')

		fm.writelines('A) ../'+name+'.fits \n')		# Input data Image (FITS file)
		fm.writelines('B) '+name+'.out.fits \n')	# Output data image block
		fm.writelines('C) none \n')	# Sigma image name
		fm.writelines('D) ../psf.fits \n')
		fm.writelines('E) 1 \n')
		fm.writelines('F) ../'+name+'.mask.fits \n')	# Bad pixel mask
		fm.writelines('G) none \n')
		fm.writelines('H) '+str(int(i-thumb))+' '+str(int(i+thumb))+' '+str(int(j-thumb))+' '+str(int(j+thumb))+' \n')	# Image region to fit (xmin, xmax, ymin, ymax)
		fm.writelines('I) '+str(thumb*2)+' '+str(thumb*2)+' \n')	# Size of the convolution box (x, y)
		fm.writelines('J) 21.5814 \n')					# Magnitude photometric ZP
		fm.writelines('K) 0.6 0.6 \n')				# Plate scale (dx dy)	[arcsec per pixel]
		fm.writelines('O) regular \n')					# Display type (regular, curses, both)
		fm.writelines('P) 0 \n\n')					# Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps

		#######  Object number : 1  #######
		fm.writelines('0) expdisk \n')				# object type
		fm.writelines('1) '+str(int(i))+' '+str(int(j))+' 1 1 \n')	# position x, y
		fm.writelines('3) '+str(ngc_match.iloc[0,16])+' 1 \n')	# Integrated magnitude
		fm.writelines('4) '+str(thumb/5)+' 1 \n')			# R_e (half-light radius)	[pix]
		fm.writelines('5) 1 1 \n')				# Sersic index n 
		fm.writelines('9) '+str(boa)+' 1 \n') 			# axis ratio (b/a)
		fm.writelines('10) '+pa+' 1 \n')				# position angle (PA) [deg: Up=0, Left=90]
		fm.writelines('Z) 0 \n\n')				# output option (0 = resig., 1 = Don't subtract)

		#######  Object number : 2  #######
        fm.writelines('0) devauc \n')                          # object type
        fm.writelines('1) '+str(int(i))+' '+str(int(j))+' 1 1 \n')        # position x, y
        fm.writelines('3) '+str(ngc_match.iloc[0,16])+' 1 \n')    # Integrated magnitude
        fm.writelines('4) '+str(thumb/5)+' 1 \n')                 # R_e (half-light radius)       [pix]
        fm.writelines('5) 1 1 \n')                              # Sersic index n
        fm.writelines('9) '+str(boa)+' 1 \n')                   # axis ratio (b/a)
        fm.writelines('10) '+pa+' 1 \n')                          # position angle (PA) [deg: Up=0, Left=90]
        fm.writelines('Z) 0 \n\n')                                # output option (0 = resig., 1 = Don't subtract)

        #######  Object number : 3  #######
#                fm.writelines('0) psf \n')                          # object type
#                fm.writelines('1) '+str(int(i))+' '+str(int(j))+' 1 1 \n')        # position x, y
#                fm.writelines('3) '+str(ngc_match.iloc[0,16]+3)+' 1 \n')    # Integrated magnitude
#                fm.writelines('Z) 0 \n\n')                                # output option (0 = resig., 1 = Don't subtract)

        #######  Object number : 4  #######
        fm.writelines('0) sky \n')                          	# object type
        fm.writelines('1) 0 1 \n')        			# sky BG at center of fitting regions [ADUs]
        fm.writelines('2) 0.0 0 \n')    			# dsky/dx (sky gradient in x)
        fm.writelines('3) 0.0 0 \n')    			# dsky/dy (sky gradient in y)
        fm.writelines('Z) 0')                                # output option (0 = resig., 1 = Don't subtract)

        fm.close()
