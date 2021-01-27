import numpy as np
from multiprocessing import Process, Pool, Lock
from astropy.io import fits

# def dqm_check(idx_start, idx_end, is_dqm_check_r):
#     for i in range(idx_start, idx_end):
#         if i % 3:
#             is_dqm_check_r[i] = False
work_dir=("/Users/duhokim/work/abell/")

single_dqm_fn = [['A2399_ud_300_0236_20aug.fits', 'A2399_gd_60_0239_19aug.fits', 'A2399_rd_300_0142_19aug.fits'],
             ['A2670_ud_300_0409_21aug.fits', 'A2670_gd_60_0420_21aug.fits', 'A2670_rd_300_0351_19aug.fits'],
             ['A3716_ud_300_2357_21aug.fits', 'A3716_gd_300_0030_21aug.fits', 'A3716_rd_300_0418_19aug.fits']]

# lock = Lock()
is_dqm_check_r = np.ones(10)

def dqm_check(idx):
    hdu_r = fits.open(work_dir + 'fits/best_single/' + single_dqm_fn[k][2])
    # for i in range(idx_start*10, (idx_start+1)*10):
    #     if i % 3:
    #         is_dqm_check_r[i] = False
    is_dqm_check = np.ones(10, dtype=bool)
    # lock.acquire()
    # try:
    is_dqm_check[idx*ni:(idx+1)*ni] = False
    # finally:
        # lock.release()
    return k



# dqm_check(0, 10)
# print(is_dqm_check_r[0:10])
for k in range(0, 3):
    num_process = 10
    ni = 10
    # hdu_temp = hdu_r
    # procs = []
    # for i in range(0, num_process):
    #     proc = Process(target=dqm_check, args=(i,))
    #     procs.append(proc)
    #     proc.start()
    #
    # for proc in procs:
    #     proc.join()
    with Pool(num_process) as p:
        is_dqm_check_r = p.map(dqm_check, range(num_process))
    # Process(target=dqm_check, args=(0))
    # dqm_check(0)
    print(is_dqm_check_r)
