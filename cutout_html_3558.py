import abell_cluster_module as ab
import webbrowser
import glob
import os
from astropy.io import ascii
from os import path

fn = ab.work_dir+'html/A3558_gala/A3558.html'
# fn2 = ab.work_dir+'html/A3558/A3558_id.txt'

f = open(fn, 'w')
# f2 = open(fn2, 'w')

message = """<html>
<head>
<meta name="viewport" content="width=device-width, initial-scale=1">
<style>
* {
  box-sizing: border-box;
}

.row {
  display: flex;
}

/* Create three equal columns that sits next to each other */
.column {
  flex: 25%;
  padding: 5px;
}

img {
  display: block;
  margin-left: auto;
  margin-right: auto;
  width: 40%;
}

</style>
</head>
<body>

<h1>Abell 3558 member (Quintana+2020) galaxies</h1>
<h2>For each ID galaxy,</h2>
<h2>First row is for cutout images in u, g, r, and color (Red: r, Green: g, and Blue: u).</h2>
<h2>Second row is for GALFIT fitting results --- 1st and 2nd: Original r-band cutout images with different stretch, 3rd: 
Two-component model (expdisk + devauc), 4th: Residual image (original - model) </h2>
"""

# os.chdir('/Users/duhokim/work/abell/pics/A3558_cutout_rgb_sqrt/')
# imgs = glob.glob('*_rgb_sqrt_0_1.png')
# thisset = set()
# for i in range(0, len(imgs)):
#     fn_split = imgs[i].split('_')
#     thisset.add(fn_split[0])
# os.chdir('/Users/duhokim/work/abell/galfit/A3558_no_weight/png/')
# galfit = glob.glob('*.png')
# for i in range(0, len(galfit)):
#     fn_split = galfit[i].split('.')
#     thisset.add(fn_split[0])

rgb_dir = ab.work_dir + 'pics/A3558_cutout_rgb_sqrt/'
sample = ascii.read(ab.work_dir+'html/A3558/A3558_id.txt', format='csv')
total_num = len(sample.columns)-1   # last row is N/A

# for x in sorted(thisset):
for i in range(0, total_num):
    id_split = sample.colnames[i].split('\'')
    id_txt = id_split[1]
    if not path.exists(rgb_dir + f'{id_txt}_rgb_sqrt_0_1.png'):  # only for galaxies with rgb image
        continue
    else:
        div = f"""
        <h2>ID: {id_txt} </h2>
    
        <div class="row">
            <div class="column">
            <img src="{id_txt}_u.png" onerror="this.onerror=null; this.src='unnamed.png'" alt="u" style="width:100%">
          </div>
          <div class="column">
            <img src="{id_txt}_g.png" onerror="this.onerror=null; this.src='unnamed.png'" alt="g" style="width:100%">
          </div>
          <div class="column">
            <img src="{id_txt}_r.png" onerror="this.onerror=null; this.src='unnamed.png'" alt="r" style="width:100%">
          </div>
          <div class="column">
            <img src="{id_txt}_rgb.png" onerror="this.onerror=null; this.src='unnamed.png'" alt="rgb" style="width:100%">
          </div>
          <div class="column">
            <img src="{id_txt}_rgb2.png" onerror="this.onerror=null; this.src='unnamed.png'" alt="rgb" style="width:130%">
          </div>
        </div>
        <div class="row">
            <div class="column">
            <img src="{id_txt}_gala.png" onerror="this.onerror=null; this.src='unnamed.png'; this.style='width=1%'" alt="galfit" style="width:70%" 
                margin-left="auto" margin-right="auto" >
          </div>
        </div>
        """
        message+=div

message+="""
</body>
</html>"""

f.write(message)
# f2.write(str(sorted(thisset)))
f.close()
# f2.close()

filename = 'file://'+fn
webbrowser.open_new_tab(filename)
# os.chdir('/Users/duhokim/Documents/Github/py_codes/')