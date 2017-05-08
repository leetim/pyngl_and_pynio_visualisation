import sys
import numpy as np
import Nio
import Ngl
import os
from interpolation import *
from datetime import *
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

class Wave:
    def __init__(self, u, v, lat_u, lon_u, lat_v, lon_v):
        self.u = u
        self.v = v
        self.lat_u = lat_u
        self.lon_u = lon_u
        self.lat_v = lat_v
        self.lon_v = lon_v
        self.time = np.linspace(0, len(u)-1, len(u))
    def __add__(self, other):
        tu = np.zeros((self.u.shape[0] + other.u.shape[0], self.u.shape[1], self.u.shape[2]))
        tu[:len(self.u),:,:] = self.u
        tu[len(self.u):len(self.u)+len(other.u),:,:] = other.u
        del self.u
        self.u = tu
        tv = np.zeros((self.v.shape[0] + other.v.shape[0], self.v.shape[1], self.v.shape[2]))
        tv[:len(self.v),:,:] = self.v
        tv[len(self.v):len(self.v)+len(other.v),:,:] = other.v
        del self.v
        self.v = tv
        self.time = np.linspace(0, len(self.time)+len(other.time)-1, len(self.time)+len(other.time))
        return self

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

def read_wave(f1):
    input_file = Nio.open_file(f1, "r")
    u = input_file.variables["u"][:, 31, :,:]
    v = input_file.variables["v"][:, 31, :,:]
    # print v
    lat_u = input_file.variables["lat_u"][:, 0]
    lon_u = input_file.variables["lon_u"][0, :]
    lat_v = input_file.variables["lat_v"][:, 0]
    lon_v = input_file.variables["lon_v"][0, :]

    return Wave(u, v, lat_u, lon_u, lat_v, lon_v)

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


dt0 = datetime(1800, 1, 1, 0, 0, 0)
dt_start = datetime(2004, 9, 8, 0, 0, 0)
start = (dt_start - dt0).days*24

ufile = Nio.open_file("uw2004.nc", "r")
vfile = Nio.open_file("vw2004.nc", "r")
u = ufile.variables["uwnd"]
v = vfile.variables["vwnd"]

lat = vfile.variables["lat"][35:0:-1]
lon = vfile.variables["lon"][55:100]
time = vfile.variables["time"]

min_t = np.min(A(list(map(int, time[:]))))
ind = (start-min_t)/6


time = A(list(map(int, time[ind-5:ind+9]-time[ind])))
min_t = np.min(time)
max_t = np.max(time)


u_arr = u[ind-5:ind+9, 0, 35:0:-1, 55:100]*u.scale_factor + u.add_offset
v_arr = v[ind-5:ind+9, 0, 35:0:-1, 55:100]*v.scale_factor + v.add_offset
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
