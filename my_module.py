"""
Functions for personal use
written by https://github.com/DuhoKim
"""

from astropy import coordinates as coords
import astropy.units as u
import numpy as np
from datetime import datetime
from astropy.cosmology import Planck18
from astropy.coordinates import Distance, SkyCoord

def closest(cat1, pos1):
    best_dist = 360.0 * 60.0 * u.arcminute
    best_index = 0
    for i in range(i, len(cat1)):
        dist = pos1.separation(cat1[i]).arcminute
        if dist < best_dist:
            best_dist = dist
            best_index = i
    return best_index, best_dist

def find_nearest(array, value):
    array = np.array(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx

# When there is a increasing or decreasing
# arrays of x corresponding to y,
# find a x0 that corresponds to y0 which is in the middle
def find_x0(x, y, y0):
    grt_y = np.where(y - y0 > 0)[0]  # index where y > y0
    if y[0] < y[-1]:                # increasing
        if len(grt_y) == 0:                     # y0 is smaller than all y
            x0 = x[0]
        elif grt_y[0] == len(x) - 1:       # y0 is greater than all y
            x0 = x[-1]
        else:
            a = y[grt_y][0] - y0
            b = y0 - y[grt_y][-1]
            e = x[grt_y][0] - x[grt_y][-1]
            dx = e * b / (a + b)
            x0 = x[grt_y][-1] + dx
    else:                                   # decreasing
        if grt_y[-1] == len(x)-1:                     # y0 is smaller than all y
            x0 = x[-1]
        elif len(grt_y) == 0:                 # y0 is greater than all y
            x0 = x[0]
        else:
            a = y[grt_y][-1] - y0
            b = y0 - y[grt_y[-1]+1]
            e = x[grt_y[-1]+1] - x[grt_y[-1]]
            dx = e * a / (a + b)
            x0 = x[grt_y[-1]] + dx
    return x0

def print_time():
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("Current Time =", current_time)

####################### Jarrett+03 ########################################################
#  Modified logarithmic visualization method P' = sqrt( log(1+P/(n*sig))  ), where P is the pixel
# intensity value, sig is the image rms "noise",
#  and n is a threshold throttle (w/ satisfactory values between 5 and 10
def jarrett(pix_values, sigma, n):
    return np.sqrt(np.log10(1+pix_values/(n*sigma)))
############################################################################################

def clcendist(z1, z2, coord1, coord2):
    dz1 = Distance(z=z1, cosmology=Planck18)
    dz2 = Distance(z=z2, cosmology=Planck18)
    c1 = SkyCoord(coord1, distance=dz1)
    c2 = SkyCoord(coord2, distance=dz2)
    sep = c1.separation_3d(c2)
    return sep.Mpc

