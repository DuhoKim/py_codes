import in_place

cat_dir=("/Users/duhokim/old_macbook_documents/work/cat/SouthernStandardStars/remove_asterisk/")

ss_name = ['220100', 'E8-A', 'LSE_259', '190000', 'E5_a', 'E6_a']
ss_cat_name = ['220100-300000', 'E8_a', 'LSE_259', '190000-295600', 'E5_a', 'E6_a', 'E4_a']

for i in range(0, len(ss_cat_name)):
    with in_place.InPlace(cat_dir+ss_cat_name[i]+'.dat.trimmed') as f:
        for line in f:
            if '*' not in line:
                f.write(line)

            # line = line.replace('*', '')
            # f.write(line)
