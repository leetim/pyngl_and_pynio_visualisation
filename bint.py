# import matplotlib.pyplot as plt
# import scipy as sp
import sys
import numpy as np
import Nio
import Ngl
import os
from interpolation import *
from draw_field import *
from read_cdf import *
from write_cdf import *
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
# print u

################################################################################
################################################################################

# dirc  = Ngl.pynglpath("data")
# ufile = Nio.open_file(os.path.join(dirc,"cdf","Ustorm.cdf"),"r")  # Open two netCDF files.
# vfile = Nio.open_file(os.path.join(dirc,"cdf","Vstorm.cdf"),"r")
# u = ufile.variables["u"]
# v = vfile.variables["v"]
# ua = u[0,:,:]
# va = v[0,:,:]
# print np.shape(ua)


################################################################################
######################## Reading from interpoled #######################################
################################################################################

# int_file = Nio.open_file("wind2004.nc", "r")
# uarr = int_file.variables["u"][0,:,:]
# varr = int_file.variables["v"][0,:,:]
# lat = int_file.variables["lat"][:]
# lon = int_file.variables["lon"][:]

w = read_from_file("come_temp.nc")
for i in range(40):
    print_vector_field(w.u[i,:,:], w.v[i,:,:], w.lat[:], w.lon[:], "./im/im{}{}.png".format(i/10, i%10), "png")

exit()
new_lat = np.linspace(38.0, 65.0, 200)
new_lon = np.linspace(115.0, 172.0, 200)
new_time = np.linspace(min_t+4, (min_t+4) + 121.0/6, 40)

# new_time = np.linspace(min_t+4, (min_t+4) + 121.0/6, 2)
# new_lat = np.linspace(58.251405, 60.420895, 10)
# new_lon = np.linspace(142.927289, 146.443603, 12)
uarr = np.ones(np.shape(uar))
varr = np.ones(np.shape(var))
for i in range(len(uar)):
    for j in range(len(uar[0])):
        uarr[i][j] = uar[i][j]
        varr[i][j] = var[i][j]
u_new = SPLINE_3D(u_arr, time, lat, lon, new_time, new_lat, new_lon)
v_new = SPLINE_3D(v_arr, time, lat, lon, new_time, new_lat, new_lon)
# print u_new

write_to_cdf("come.nc", u_new, v_new, new_lat, new_lon, new_time)

################################################################################
################################################################################

# u_new = BIG_INTERPOL(u_arr, time, lat, lon, new_time, new_lat, new_lon)
# v_new = BIG_INTERPOL(v_arr, time, lat, lon, new_time, new_lat, new_lon)
# exit()
# return 0
# for i in range(n)

# new_u = np.zeros((24, len(lat), len(lon)))
# new_v = np.zeros((24, len(lat), len(lon)))
# for i in range(len(lat)):
#     for j in range(len(lon)):
#         print (i, j)
#         new_u[:, i, j] = LAGRANGE(u_arr[:, i, j], time, new_time)
#         new_v[:, i, j] = LAGRANGE(v_arr[:, i, j], time, new_time)
# print np.shape(new_u)
# fig = plt.figure()
# plt.streamplot(new_lon, new_lat, new_u[0], new_v[0])
# plt.title('Simple stream plot')
# plt.grid(True)
# plt.show()
