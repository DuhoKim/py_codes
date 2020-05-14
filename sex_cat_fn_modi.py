#####################################################
# Python script that modify SEx cat file names
# written by Duho Kim (22 Apr 20)
######################################################
import os
from shutil import copyfile

sex_dir = '/Users/dkim108/Documents/work/sex/gals/'

dates = ['19', '20', '21']

for i in range(0, len(dates)):
	for file in os.listdir(sex_dir+dates[i]+'aug'):
		if not file.endswith('aug.cat') and not file.startswith('.'):
			copyfile(sex_dir+dates[i]+'aug/'+file, sex_dir+file[:5]+'/'+file[:-4]+'_'+dates[i]+'aug.cat')