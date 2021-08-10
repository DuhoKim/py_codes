import abell_cluster_module as ab
import xlsxwriter
import os.path
from os import path
from astropy.io import ascii
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy import wcs
import shutil

copy_dir = ab.work_dir + 'gui/A3558/'

band_dir = ab.work_dir + 'pics/A3558_cutout_grey_and_Q20_str04_x05/'
rgb_dir = ab.work_dir + 'pics/A3558_cutout_rgb_sqrt/'

sha_cat = ascii.read("/Users/duhokim/work/abell/spec/Shapley/catalog.dat")

hdu_r_stack = fits.open(ab.work_dir + 'fits/stacked/A3558_rsi.fits')

# images list
id_col = 0      # column number for ID number
img_col = 1     # column number for images
cat_col = 2     # column number for categorization

scale = 0.6

sample = ascii.read(ab.work_dir+'html/A3558/A3558_id.txt', format='csv')
total_num = len(sample.columns)-1   # last row is N/A
row_num = 0
rgb_num = 0     # number of galaxies with rgb
gala_num = 0    # to measure GALAPAGOS success rate
wsn = 0         # work sheet number
# comb_num = 0    # to measure combined success rate
#
# tile_no = [0] * total_num
# galfit_no = [0] * total_num

# coord_sha = SkyCoord(f"{sha_cat['col3']}:{sha_cat['col4']}:{sha_cat['col5']} "
#                      f"{sha_cat['col6']}:{sha_cat['col7']}:{sha_cat['col8']}",
#                      unit=(u.hourangle, u.deg))
#
# for i in range(1, 10):
#     sex_cat1 = ascii.read(f'/Users/duhokim/work/abell/galapagos/A3558_rsr_all/t{tile_no}/t{tile_no}r.outcat')


# create an new Excel file and add a worksheet
workbook = xlsxwriter.Workbook(ab.work_dir+'xlsx/A3558_gala_cat_851_zp21_100_fan_3gala_dummy.xlsx')
fn_list = ab.work_dir+'gui/A3558/A3558.list'

# Add a format for the header cells.
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

heading1 = 'ID'
heading2 = 'u'
heading3 = 'g'
heading4 = 'r'
heading5 = 'rgb1'
heading6 = 'rgb2'
heading7 = 'galfit'
heading8 = 'General'
heading9 = 'Features'

exp1 = '''For 'General' classification: \n
“Elliptical“ denotes the normal elliptical and \n
“S0“ lenticular galaxies \n
'''

exp2 = '''“spiral” indicates normal spiral galaxies. \n
“irregular” is for irregular galaxies. 
'''

exp3 = '''For 'Features' classification: \n
“none” shows no peculiar feature. \n
“Post-merger” denotes galaxies showing disturbed features, e.g., asymmetric structures, faint features, discontinuous halo structure, rings, and dust lanes. These are considered galaxy merger remnants. \n
'''

exp4 = ''' “Interacting” galaxies also exhibit disturbed features as well as PM galaxies, but they are close companions and interact with each other, i.e. they are ongoing merger candidates. \n
“Pair” refers to galaxies with a close companion without any disturbed features. \n
“Jelly-fish” denotes galaxies showing ram-pressure stripping.'''

def init_ws(wsn):
    ws = workbook.add_worksheet(f'{wsn * 100}~{(wsn + 1) * 100}')

    # resize cells
    ws.set_column(0, 0, 20)  # Set 'ID' column's width [# of character]
    ws.set_column(1, 8, 40)  # Set 'A, B, C' column's width [# of character]
    ws.set_column(6, 6, 140)
    ws.set_column(7, 8, 40, id_format)
    ws.set_column(9, 9, 250, id_format)
    ws.set_default_row(240)  # Set the default row height

    ws.write('A1', heading1, header_format)
    ws.write('B1', heading2, header_format)
    ws.write('C1', heading3, header_format)
    ws.write('D1', heading4, header_format)
    ws.write('E1', heading5, header_format)
    ws.write('F1', heading6, header_format)
    ws.write('G1', heading7, header_format)
    ws.write('H1', heading8, header_format)
    ws.write('I1', heading9, header_format)

    ws.write('J2', exp1, id_format)
    ws.write('J3', exp2, id_format)
    ws.write('J4', exp3, id_format)
    ws.write('J5', exp3, id_format)

    ws.data_validation(f'H2:H101', {'validate': 'list',
                                      'source': ['elliptical','S0',
                                                 'spiral', 'irregular','N/A']})

    ws.data_validation(f'I2:I101', {'validate': 'list',
                                      'source': ['none', 'asymmetric', 'warped', 'fan', 'interacting',
                                                 'jellyfish']})

    return ws

list_f = open(fn_list, 'w')

for i in range(0, total_num):

    id_split = sample.colnames[i].split('\'')
    id_txt = id_split[1]
    if not path.exists(rgb_dir + f'{id_txt}_rgb_sqrt_0_1.png'):  # only for galaxies with rgb image
        continue
    else:
        row_num = row_num + 1
        rgb_num = rgb_num + 1
        list_f.write(f'{id_txt} \n')
        continue

    if row_num == 101 or i == 0:
        ws = init_ws(wsn)
        row_num = 1
        wsn = wsn + 1

    ws.write(row_num, id_col, f'ID_{id_txt}', id_format)

    fn_img = band_dir + f'{id_txt}_u.png'
    if path.exists(fn_img):
        ws.insert_image(row_num,
                               1,
                               fn_img,
                               {'x_scale': scale, 'y_scale': scale,
            #                    'x_offset': 5, 'y_offset': 5,
                                'object_position': 1})
        shutil.copyfile(fn_img, copy_dir+f'{id_txt}_u.png')
    fn_img = band_dir + f'{id_txt}_g.png'
    if path.exists(fn_img):
        ws.insert_image(row_num,
                               2,
                               fn_img,
                               {'x_scale': scale, 'y_scale': scale,
            #                    'x_offset': 5, 'y_offset': 5,
                                'object_position': 1})
        shutil.copyfile(fn_img, copy_dir + f'{id_txt}_g.png')
    fn_img = band_dir + f'{id_txt}_r.png'
    if path.exists(fn_img):
        ws.insert_image(row_num,
                               3,
                               fn_img,
                               {'x_scale': scale, 'y_scale': scale,
            #                    'x_offset': 5, 'y_offset': 5,
                                'object_position': 1})
        shutil.copyfile(fn_img, copy_dir + f'{id_txt}_r.png')
    fn_img = ab.work_dir+f'html/A3558/{id_txt}_rgb.png'
    if path.exists(fn_img):
        ws.insert_image(row_num,
                               4,
                               fn_img,
                               {'x_scale': scale, 'y_scale': scale,
            #                    'x_offset': 5, 'y_offset': 5,
                                'object_position': 1})
        shutil.copyfile(fn_img, copy_dir + f'{id_txt}_rgb.png')
    fn_img = rgb_dir + f'{id_txt}_rgb_sqrt_0_1.png'
    if path.exists(fn_img):
        ws.insert_image(row_num,
                               5,
                               fn_img,
                               {'x_scale': scale, 'y_scale': scale,
            #                    'x_offset': 5, 'y_offset': 5,
                                'object_position': 1})
        shutil.copyfile(fn_img, copy_dir + f'{id_txt}_rgb2.png')
    # GALFIT result
    # fn_img_gal = ab.work_dir + f'galfit/A3558_no_weight/png/{id_txt}.png'
    # if path.exists(fn_img_gal):
    #     ws.insert_image(row_num + 1,
    #                     6,
    #                     fn_img_gal,
    #                     {'x_scale': 0.855, 'y_scale': 0.855,
    #                      #                    'x_offset': 5, 'y_offset': 5,
    #                      'object_position': 1})
    #     gal_num = gal_num + 1
    #     comb_num = comb_num + 1

    # GALAPAGOS result
    sha_ind = np.where(sha_cat['col1'] == int(id_txt))[0][0]
    coord_sha = SkyCoord(f"{sha_cat['col3'][sha_ind]}:{sha_cat['col4'][sha_ind]}:{sha_cat['col5'][sha_ind]} "
                         f"{sha_cat['col6'][sha_ind]}:{sha_cat['col7'][sha_ind]}:{sha_cat['col8'][sha_ind]}",
                         unit=(u.hourangle, u.deg))
    cel_coord = [[coord_sha.ra.value, coord_sha.dec.value], [0, 0]]
    tile_no = -1
    for jj in range(1, len(hdu_r_stack)):
        # read WCS
        w = wcs.WCS(hdu_r_stack[jj].header)
        pixcrd = w.wcs_world2pix(cel_coord, 1)
        # if w.footprint_contains(sky_coord):
        if (pixcrd[0][0] > 0) & (pixcrd[0][0] < hdu_r_stack[jj].shape[1]) & \
                (pixcrd[0][1] > 0) & (pixcrd[0][1] < hdu_r_stack[jj].shape[0]):
            tile_no = jj
    if tile_no > -1:
        sex_cat = ascii.read(f'/Users/duhokim/work/abell/galapagos/A3558_rsr_cat_all/t{tile_no}/t{tile_no}r.outcat')
        coords_sex = SkyCoord(sex_cat['col13'], sex_cat['col14'], unit='deg')
        d2d = coord_sha.separation(coords_sex)
        matched_sex = (d2d.arcsec < 1)
        if sum(matched_sex):
            galfit_id = sex_cat['col1'][matched_sex][0]
        else:
            continue

        fn_img_gala = f'/Users/duhokim/work/abell/galapagos/A3558_rsr_cat_all/t{tile_no}/galfit/png2/t{tile_no}{galfit_id}.png'
        if path.exists(fn_img_gala):
            gala_num = gala_num + 1
            # if not path.exists(fn_img_gal):
            ws.insert_image(row_num,
                            6,
                            fn_img_gala,
                            {'x_scale': 0.855, 'y_scale': 0.855,
                             #                    'x_offset': 5, 'y_offset': 5,
                             'object_position': 1}
                            )
            shutil.copyfile(fn_img_gala, copy_dir + f'{id_txt}_gala.png')
                # comb_num = comb_num + 1

list_f.close()
workbook.close()


