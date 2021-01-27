from astropy.io import ascii

work_dir=("/Users/duhokim/work/abell/sex/cat/")
clusters=['A2399', 'A2670', 'A3716']
ver = 'v1.2'
mag_lim = 19

for cluster in clusters:
    fn = work_dir + 'DECam_19_21_aug_2014_stacked_SEx_cat_{}_match_rad_1as_Gal_ext_corrected_{}'.format(cluster, ver)
    with open(fn+'_lt_mag19.reg', 'w') as new_rg:
        cat = ascii.read(fn+'.txt')
        reg = ascii.read(fn+'.reg')
        reg['col2'] = reg['col2'].astype('U')
        reg['col7'] = 'dash=9'
        for i in range(0, len(cat)):
            if cat[i]['MAG_AUTO_r'] < mag_lim or cat[i]['MAG_AUTO_r_stack'] < mag_lim:
                reg[i]['col5'] = reg[i]['col5'][:10]
                new_rg.writelines(' '.join(reg[i]))
                new_rg.write(' \n')