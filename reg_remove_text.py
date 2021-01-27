import in_place

work_dir=("/Users/duhokim/work/abell/")

# fn = work_dir + 'sex/cat/DECam_19_21_aug_2014_single_best_exposure_SEx_cat_A2670_match_rad_1as_Gal_ext_corrected_v1.2_no_text.reg'
fn = work_dir + 'sex/cat/stack/A2670_usi_no_text.reg'
#fn = work_dir + 'sex/cat/stack/A2399_gsi_no_text.reg'

with in_place.InPlace(fn) as f:
    for line in f:
        if line == 'fk5\n':
            f.write(line)
            continue
        a = line.split('{')
        b = a[1].split('}')
        new_line = a[0]+'\'\''+b[1]
        f.write(new_line)

        # a = line.split('{')
        # b = a[1].split('}')
        # new_line = a[0]+'\'\''+b[1]
        # f.write(new_line)
