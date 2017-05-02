import sys
import numpy as np
import Nio
import Ngl
import os
from interpolation import *
# import random as rnd
# from scipy.optimize import leastsq, linprog
# from scipy import stats
# import math
# import pylab
# from mpl_toolkits.mplot3d import Axes3D
# import netCDF4 as nc
import threading as th
from scipy.interpolate import RegularGridInterpolator as RGI
from scipy.interpolate import interp1d

A = np.array
T = np.transpose

class Wind:
    def __init__(self, u, v, lat, lon, time):
        self.u = u
        self.v = v
        self.lat = lat
        self.lon = lon
        self.time = time

################################################################################
######################## Reading from uw #######################################
################################################################################
def read_from_file2(f1, f2):
    ufile = Nio.open_file(f1, "r")
    vfile = Nio.open_file(f2, "r")
    u = ufile.variables["uwnd"][:,:,:]
    v = vfile.variables["vwnd"][:,:,:]
    # print v
    lat = vfile.variables["lat"][:]
    lon = vfile.variables["lon"][:]
    time = vfile.variables["time"][:]
    return Wind(u, v, lat, lon, time)

def read_from_file(f, u_name = 'u', v_name = 'v'):
    input_file = Nio.open_file(f, "r")
    u = input_file.variables[u_name][:,:,:]
    v = input_file.variables[v_name][:,:,:]
    # print v
    lat = input_file.variables["lat"][:]
    lon = input_file.variables["lon"][:]
    time = np.zeros(len(u))
    # time = input_file.variables["time"][:]
    return Wind(u, v, lat, lon, time)


ufile = Nio.open_file("uw2004.nc", "r")
vfile = Nio.open_file("vw2004.nc", "r")
u = ufile.variables["uwnd"]
v = vfile.variables["vwnd"]
# print v
lat = vfile.variables["lat"][35:0:-1]
lon = vfile.variables["lon"][55:100]
print len(lat)
time = vfile.variables["time"]
time = A(list(map(int, time[:]-time[0])))
min_t = np.min(time)/6 + 248*4
max_t = min_t + 7*4
# print (min_t, max_t)
# new_time =
time = np.linspace(min_t, max_t, 28)
u_arr = u[min_t:max_t, 0, 35:0:-1, 55:100]*u.scale_factor + u.add_offset
v_arr = v[min_t:max_t, 0, 35:0:-1, 55:100]*v.scale_factor + v.add_offset
uar = u_arr[0, ::, ::]#*u.scale_factor + u.add_offset
var = v_arr[0, ::, ::]#*v.scale_factor + v.add_offset

# new_lat = lat
# new_lon = lon
# u_new[:, 2, 5] = LAGRANGE1(u_arr[:, 2, 5], time, new_time)
# v_new[:, 2, 5] = LAGRANGE1(v_arr[:, 2, 5], time, new_time)
# print u_new[:, 2, 5]
# print u_arr[:, 2, 5]
# # exit()
# for i in range(u_new.shape[1]):
#     for j in range(u_new.shape[2]):
#         print (i, j)
#         pers += step
#         if (np.abs(last_pers - pers) > 0.1):
#             print pers
#             last_pers = pers
#         u_new[:, i, j] = LAGRANGE1(u_arr[:, i, j], time, new_time)
#         v_new[:, i, j] = LAGRANGE1(v_arr[:, i, j], time, new_time)
# u_new = SPLINE_3D(u_arr, time, lat)
# print (u_new.shape)
# print np.max(np.abs(u_new - u_arr))
