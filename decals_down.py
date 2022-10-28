#####################################################
# Python script that downloads DECaLS tractor catalogs
# written by Duho Kim (10 Jun 22)
######################################################
import os
import abell_cluster_module as ab
import math
import requests
from bs4 import BeautifulSoup

cur_dir = os.getcwd()

for i in range(0, 7):
    cat_dir = f'/Users/duhokim/work/abell/cat/legacy/{ab.clusters[i]}'
    os.chdir(cat_dir)

    ra = ab.coords_cl_cen[i].ra.value
    dec = ab.coords_cl_cen[i].dec.value

    ras = range(math.floor(ra*10)-10, math.floor(ra*10)+10)
    decs = range(math.floor(-dec*10)-10, math.floor(-dec*10)+10)

    def listFD(url, ext=''):
        page = requests.get(url).text
        soup = BeautifulSoup(page, 'html.parser')
        return [url + node.get('href') for node in soup.find_all('a') if node.get('href').endswith(ext)]

    dir_ras = range(math.floor(ra)-1, math.floor(ra)+2)
    for dir_ra in dir_ras:
        url = f"https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/dr9/south/tractor/{dir_ra}/"
        for file in listFD(url, 'fits'):
            file_ra = int(file[-13:-9])
            file_mp = file[-9:-8]
            file_dec = int(file[-8:-5])
            if (file_mp == 'm') and (file_ra in ras) and (file_dec in decs):
                os.system(f"wget {file}")

os.chdir(cur_dir)