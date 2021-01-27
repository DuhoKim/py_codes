from astropy.io import ascii

work_dir=("/Users/duhokim/work/abell/cat/")
short_cat_dir=("/Volumes/APPLE SSD/Server/work/sex/best_single/")
out_dir=("/Users/duhokim/work/abell/sex/reg/")

fns = ['n32_m05', 'n32_m0005', 'n64_m005', 'n16_m005', 'default', 'deblend']
color = ['blue', 'cyan', 'magenta', 'yellow', 'green', 'red']
# white
# black
# red
# green
# blue
# cyan
# magenta
# yellow

### from Topcat catalog
# fn = work_dir + 'g-r_blue_subset_A2670'
# with open(fn+'.reg', 'w') as new_rg:
#     cat = ascii.read(fn)
#     reg['col2'] = cat['col2'].astype('U')
#     reg['col7'] = 'dash=9'
#     for i in range(0, len(cat)):
#         if cat[i]['MAG_AUTO_r'] < mag_lim or cat[i]['MAG_AUTO_r_stack'] < mag_lim:
#             reg[i]['col5'] = reg[i]['col5'][:10]
#             new_rg.writelines(' '.join(reg[i]))
#             new_rg.write(' \n')

with open(out_dir + 'g-r_line_A754.reg', 'w') as new_rg:
    new_rg.writelines('# Region file format: DS9 version 4.1\n')
    new_rg.writelines(f'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 \
    highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    new_rg.writelines('fk5\n')

    cat = ascii.read("/Users/duhokim/work/abell/sex/cat/g-r_line_A754")
    for j in range(0, len(cat)):
        new_rg.writelines("circle({:12.7f}, {:12.7f}, {:7.3f}\") # text=\'{}\' \n".format(
            cat[j]['col2'],
            cat[j]['col3'],
            cat[j]['col4'] *  5,
            cat[j]['col1']))

# for i in range(0, len(fns)):
#     with open(out_dir + 'A754_' + fns[i] + '.reg', 'w') as new_rg:
#         new_rg.writelines('# Region file format: DS9 version 4.1\n')
#         new_rg.writelines(f'global color={color[i]} dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 \
#         highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
#         new_rg.writelines('fk5\n')
#
#         cat = ascii.read(short_cat_dir+'A754_rs_300_'+fns[i]+'.cat')
#         for j in range(0, len(cat)):
#             new_rg.writelines("ellipse({:12.7f}, {:12.7f}, {:7.3f}\", {:7.3f}\", {:7.3f}) # \n".format(
#                 cat[j]['ALPHA_J2000'],
#                 cat[j]['DELTA_J2000'],
#                 cat[j]['A_WORLD'] * 3600 * cat[j]['KRON_RADIUS'],
#                 cat[j]['B_WORLD'] * 3600 * cat[j]['KRON_RADIUS'],
#                 180 - cat[j]['THETA_WORLD']))