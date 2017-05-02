import sys
import numpy as np
import Nio
import Ngl
import os
from interpolation import *
from draw_field import *
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

def write_to_cdf(f, u, v, lat, lon, time, u_name = 'u', v_name = 'v'):
    nc = Nio.open_file(f, 'w')
    nc.create_dimension('x',len(lon))
    nc.create_dimension('y',len(lat))
    nc.create_dimension('t',len(time))
    nc.create_variable(u_name, 'd', ('t', 'y', 'x'))
    nc.create_variable(v_name, 'd', ('t', 'y', 'x'))
    nc.create_variable('lat', 'd', ('y',))
    nc.create_variable('lon', 'd', ('x',))
    nc.create_variable('time', 'd', ('t',))
    # nc.create_variable('time', 'd', ('time',))
    # print np.shape(u_new)
    # print len(lon)
    # print len(lat)
    # print np.max(np.abs(u_arr))
    # print nc.variables['u'].shape
    print time
    nc.variables[u_name][:,:,:] = u[:,:,:]
    nc.variables[v_name][:,:,:] = v[:,:,:]
    nc.variables['lat'][:] = lat[:]
    nc.variables['lon'][:] = lon[:]
    nc.variables['time'][:] = time[:]
    # nc.variables['time'] = new_time
    nc.close()

def write_to_grb_like(f, u, v, lat, lon, u_name = 'u', v_name = 'v'):
    nc = Nio.open_file(f, 'w')
    nc.create_dimension('x',len(lon))
    nc.create_dimension('y',len(lat))
    nc.create_variable(u_name, 'd', ('y', 'x'))
    nc.create_variable(v_name, 'd', ('y', 'x'))
    nc.create_variable('lat', 'd', ('y',))
    nc.create_variable('lon', 'd', ('x',))
    # nc.create_variable('time', 'd', ('time',))
    # print np.shape(u_new)
    # print len(lon)
    # print len(lat)
    # print np.max(np.abs(u_arr))
    # print nc.variables['u'].shape
    nc.variables[u_name][:,:] = u[:,:]
    nc.variables[v_name][:,:] = v[:,:]
    nc.variables['lat'][:] = lat[:]
    nc.variables['lon'][:] = lon[:]
    # nc.variables['time'] = new_time
    nc.close()
