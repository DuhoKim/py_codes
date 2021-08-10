import os
import glob

os.chdir('/Users/duhokim/work/abell/galfit/A3558_no_weight/')
fn = sorted(glob.glob('*.feedme'), key=os.path.getmtime, reverse=True)

for i in range(0, len(fn)):
	os.system(f"galfit {fn[i]}")
