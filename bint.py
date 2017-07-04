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
######################### Write GribLike #######################################
################################################################################

def write_for_spillmod():
    w = read_from_file("come1.nc")
    print len(w.time)
    # exit()
    for i in range(len(w.time)):
        write_to_grb_like("./grb/wrf_data0{}{}.nc".format(i/10, i%10), w.u[i,:,:], w.v[i,:,:],
                          w.lat, w.lon, "u-component_of_wind_height_above_ground",
                          "v-component_of_wind_height_above_ground")

# int_file = Nio.open_file("wind2004.nc", "r")
# uarr = int_file.variables["u"][0,:,:]
# varr = int_file.variables["v"][0,:,:]
# lat = int_file.variables["lat"][:]
# lon = int_file.variables["lon"][:]
################################################################################
############################# Draw Field #######################################
################################################################################

def Draw():
    w = read_from_file("sah_wave.nc")
    max_lim = np.max(np.sqrt(w.u**2 + w.v**2))*1.1
    for i in range(w.u.shape[0]):
        print_vector_field(w.u[i,:,:], w.v[i,:,:], w.lat[:], w.lon[:], "./sah_wave_im/im00{}{}.png".format(i/10, i%10), "png",
            title = "Wave Sakhalin {} September 2004 {}{}:{}0".format(8 + i/24/2, (i/2%24)/10, (i/2%24)%10, 3*(i%2)), max_limit = max_lim)
    w = read_from_file("sah_wind.nc")
    max_lim = np.max(np.sqrt(w.u**2 + w.v**2))*1.1
    for i in range(w.u.shape[0]):
        print_vector_field(w.u[i,:,:], w.v[i,:,:], w.lat[:], w.lon[:], "./sah_wind_im/im00{}{}.png".format(i/10, i%10), "png",
            title = "Wind Sakhalin {} September 2004 {}{}:{}0".format(8 + i/24/2, (i/2%24)/10, (i/2%24)%10, 3*(i%2)), max_limit = max_lim)

def Print_wave():
    w1 = read_wave("SPILLMOD_DATA/roms_his_0303.nc")
    w2 = read_wave("SPILLMOD_DATA/roms_his_0304.nc")
    w1 = w1+w2
    # new_time = np.linspace(0, 47, 96)
    # new_lat = np.linspace(45.77, 47.81, 30)
    # new_lon = np.linspace(140.62, 144.45, 30)
    new_time = np.linspace(0, 47, 48)
    new_lat = np.linspace(43.6, 65.0, 50)
    new_lon = np.linspace(130.1, 165.0, 50)
    # print "________"
    # print w1.time
    # print w2.time
    # print new_time
    # print "________"
    u_new = SPLINE_3D(w1.u, w1.time, w1.lat_u, w1.lon_u, new_time, new_lat, new_lon)
    v_new = SPLINE_3D(w1.v, w1.time, w1.lat_v, w1.lon_v, new_time, new_lat, new_lon)
    # u_new = np.ones((2, 20, 20))*1E37
    # v_new = np.ones((2, 20, 20))
    # i = 4
    # print_vector_field(u_new[i,:,:], v_new[i,:,:], new_lat[:], new_lon[:], "./im{}{}.png".format(i/10, i%10), "png",
    #     title = "Wave Okhotsk sea {} May 2017 {}{}:{}0".format(10 + i/24/2, (i/2%24)/10, (i/2%24)%10, 3*(i%2)), max_limit = max_lim)
    # return
    for i in range(u_new.shape[0]):
        for j in range(u_new.shape[1]):
            for k in range(u_new.shape[2]):
                if u_new[i,j,k] > 1E5:
                    u_new[i,j,k] = 0.0
    for i in range(v_new.shape[0]):
        for j in range(v_new.shape[1]):
            for k in range(v_new.shape[2]):
                if v_new[i,j,k] > 1E5:
                    v_new[i,j,k] = 0.0
    max_lim = np.max(np.sqrt(u_new**2 + v_new**2))*1.1
    for i in range(len(u_new)):
        print u_new[i,:,:]
        print v_new[i,:,:]
        print_vector_field(u_new[i,:,:], v_new[i,:,:], new_lat[:], new_lon[:], "./okhotsk_wave_im/im{}{}.png".format(i/10, i%10), "png",
            title = "Wave Okhotsk sea {} May 2017 {}{}:{}0".format(8 + i/24/2, (i/2%24)/10, (i/2%24)%10, 3*(i%2)), max_limit = max_lim)
    # write_to_cdf("okhotsk_wave_20170510.nc", u_new, v_new, new_lat, new_lon, new_time)



# exit()
################################################################################
########################## Interpolation #######################################
################################################################################
def Interpolation_write():
    # new_lat = np.linspace(46.88, 47.17, 20)
    # new_lon = np.linspace(141.7, 142.4, 20)
    # new_lat = np.linspace(43.6, 65.0, 20)
    # new_lon = np.linspace(130.1, 165.0, 20)
    new_time = np.linspace(0, 47, 96)
    new_lat = np.linspace(45.77, 47.81, 30)
    new_lon = np.linspace(140.62, 144.45, 30)
    # new_time = np.linspace(min_t+4, (min_t+4) + 121.0/6, 40)

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
    max_lim = np.max(np.sqrt(u_new**2 + v_new**2))*1.1
    for i in range(len(u_new)):
        print_vector_field(u_new[i], v_new[i], new_lat, new_lon, "./okhotsk_wave_im/00{}{}.png".format(i/10, i%10), "png",
            title = "Wind Sakhalin {} September 2004 {}{}:{}0".format(8 + i/24/2, (i/2%24)/10, (i/2%24)%10, 3*(i%2)), max_limit = max_lim)
    # write_to_cdf("sah_wind.nc", u_new, v_new, new_lat, new_lon, new_time)
    # print u_new

def read_wind_from_grib():
    time = np.linspace(-72, 48, 41)
    # print time
    Z = np.zeros
    new_time = np.linspace(0, 47, 96)
    new_lat = np.linspace(43.6, 65.0, 30)
    new_lon = np.linspace(130.1, 165.0, 30)
    f = ["%03d"%(i,) for i in range(41)]
    w = read_grib("SPILLMOD_DATA/wrfprs.%s.grb"%(f[0],))
    shape = (41, w.u.shape[0], w.u.shape[1])
    wind = Wind(Z(shape), Z(shape), None, None)
    wind.u[0,:,:] = w.u
    wind.v[0,:,:] = w.v
    wind.lat = w.lat
    wind.lon = w.lon
    for i in range(1, len(f)):
        w = read_grib("SPILLMOD_DATA/wrfprs.%s.grb"%(f[i],))
        wind.u[i,:,:] = w.u
        wind.v[i,:,:] = w.v
    wind.time = time
    w = wind
    # print time
    print "_________________________________________"
    print w.u.shape
    print w.v.shape
    print w.time.shape
    print w.lat
    print w.lon
    print new_time.shape
    print new_lat
    print new_lon
    print "_________________________________________"
    u = SPLINE_3D(w.u, w.time, w.lat, w.lon, new_time, new_lat, new_lon)
    v = SPLINE_3D(w.v, w.time, w.lat, w.lon, new_time, new_lat, new_lon)
    # write_to_cdf("okh_wind_20170510.nc", u, v, new_lat, new_lon, new_time)
        # print w.u.shape

    max_lim = np.max(np.sqrt(u**2 + v**2))
    for  i in range(u.shape[0]):
        print_vector_field(u[i,:,:], v[i,:,:], new_lat[:], new_lon[:], "./okhotsk_wind_im/im{}{}.png".format(i/10, i%10), "png",
            title = "Wind Okhotsk sea {} May 2017 {}{}:{}0".format(10 + i/24/2, (i/2%24)/10, (i/2%24)%10, 3*(i%2)), max_limit = max_lim)
        # print_vector_field(w.u[::10, ::10], w.v[::10, ::10], w.lat[::10], w.lon[::10], "wrf_20170510/arim%s.png"%(f[i],), "png",
        #     title = "Wind", max_limit = max_lim)

################################################################################
################################################################################
# read_wind_from_grib()
# Interpolation_write()
Print_wave()
# Draw()
# write_for_spillmod()

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
