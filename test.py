for ii in range(5):
    add_row = ii // 4 if ii else -1
    col444 = (ii % 4) - 1 if ii else 2
    col444 = col444 if col444 >= 0 else 0

    print(f'{3 + add_row} {col444}')
