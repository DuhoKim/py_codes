import in_place

cat_dir=("/Users/dkim108/Documents/work/cat/SouthernStandardStars/")

ss_name = ['220100', 'E8-A', 'LSE_259', '190000']
#ss_cat_name = ['220100-300000', 'E8_a', 'LSE_259', '190000-295600']
ss_cat_name = ['E5_a', 'E6_a']
for i in range(0, len(ss_cat_name)):
    with in_place.InPlace(cat_dir+ss_cat_name[i]+'.dat.trimmed') as f:
        for line in f:
            line = line.replace('*', '')
            f.write(line)
