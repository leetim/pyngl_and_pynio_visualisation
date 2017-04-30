# import matplotlib.pyplot as plt
# import scipy as sp
import sys
import numpy as np
import Nio
import Ngl
import os
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
F = lambda x, y, z: (x + y**2 + z**3)/2

n = 5
ax, ay, az = [0.0 for i in range(3)]
bx, by, bz = [1.0 for i in range(3)]
x_base = np.linspace(ax, bx, n)
y_base = np.linspace(ay, by, n)
z_base = np.linspace(az, bz, n)
x_new = np.linspace(ax, bx, 5*n)
y_new = np.linspace(ay, by, 5*n)
z_new = np.linspace(az, bz, 5*n)
y, x, z = np.meshgrid(x_base, y_base, z_base) #important!
yn, xn, zn = np.meshgrid(x_new, y_new, z_new) #important!
u = F(x, y, z)
# print u
################################################################################
################################################################################

def LAGRANGE(u, x, new_x):
    n = len(u)
    # print x
    new_u = np.zeros(len(new_x))
    FUNC = lambda x, x1, x2: (x - x2)/(x1 - x2)
    for k in range(len(new_x)):
        for i in range(n):
            temp = u[i]
            for j in range(n):
                if j != i:
                    temp *= FUNC(new_x[k], x[i], x[j])
            new_u[k] += temp
    return new_u

################################################################################
################################################################################

def spline(u, arx, x):
    n = len(u)
    a = np.zeros((n))
    b = np.zeros((n))
    c = np.zeros((n))
    d = np.zeros((n))

def lagrange_interpol(u, arx, ary, arz, x, y, z):
    n = np.shape(u)
    inds = [n[0]/2, n[1]/2, n[2]/2]
    args = [arx, ary, arz]
    vals = A([x, y, z])
    for i in range(3):
        j1 = 0
        j2 = n[i] - 1
        while j1 != j2 - 1:
            k = (j1 + j2)/2
            if (args[i][k] > vals[i]):
                j2 = k
            else:
                j1 = k
        inds[i] = (j1, j2)
    res = 0.0
    def MULT(x, i, arx):
        res = 1.0
        for j in range(len(arx)):
            if (i != j):
                res *= (x - arx[j])/(arx[i] - arx[j])
        return res
    for  i in range(len(arx)):
        for j in range(len(ary)):
            for k in range(len(arz)):
                mx = MULT(x, i, arx)
                my = MULT(y, j, ary)
                mz = MULT(z, k, arz)
                res += u[i][j][k]*mx*my*mz
    return res


def bliner_interpol(u, arx, ary, arz, x, y, z):
    n = np.shape(u)
    inds = [n[0]/2, n[1]/2, n[2]/2]
    args = [arx, ary, arz]
    vals = A([x, y, z])
    for i in range(3):
        j1 = 0
        j2 = n[i] - 1
        while j1 != j2 - 1:
            k = (j1 + j2)/2
            if (args[i][k] > vals[i]):
                j2 = k
            else:
                j1 = k
        inds[i] = (j1, j2)
    res = 0.0
    def MULT(x, i, arx, lim):
        res = 1.0
        for j in range(lim[0], lim[1]):
            if (i != j):
                res *= (x - arx[j])/(arx[i] - arx[j])
        return res
    for  i in range(inds[0][0], inds[0][1]+1):
        for j in range(inds[1][0], inds[1][1]+1):
            for k in range(inds[2][0], inds[2][1]+1):
                mx = MULT(x, i, arx, inds[0])
                my = MULT(y, j, ary, inds[1])
                mz = MULT(z, k, arz, inds[2])
                res += u[i][j][k]*mx*my*mz
    return res

def bcub_interpol(u, arx, ary, arz, x, y, z):
    n = np.shape(u)
    inds = [n[0]/2, n[1]/2, n[2]/2]
    args = [arx, ary, arz]
    vals = A([x, y, z])
    for i in range(3):
        j1 = 0
        j2 = n[i] - 1
        while j1 != j2 - 1:
            k = (j1 + j2)/2
            if (args[i][k] > vals[i]):
                j2 = k
            else:
                j1 = k
        if j1-1 >= 0 and j2+1 < n[i]:
            j1 -= 1
            j1 += 1
        else:
            if j1-2 >= 0:
                j1 -= 2
            if j2+2 < n[i]:
                j2 += 2
        inds[i] = (j1, j2)
    res = 0.0
    def MULT(x, i, arx, lim):
        res = 1.0
        for j in range(lim[0], lim[1]):
            if (i != j):
                res *= (x - arx[j])/(arx[i] - arx[j])
        return res
    for  i in range(inds[0][0], inds[0][1]+1):
        for j in range(inds[1][0], inds[1][1]+1):
            for k in range(inds[2][0], inds[2][1]+1):
                mx = MULT(x, i, arx, inds[0])
                my = MULT(y, j, ary, inds[1])
                mz = MULT(z, k, arz, inds[2])
                res += u[i][j][k]*mx*my*mz
    return res

def LAGRANGE1(u, x_cur, x_new):
    n = len(u)
    F = interp1d(x_cur, u, kind='linear')
    res = np.zeros(len(x_new))
    for k in range(len(res)):
        res[k] = F(x_new[k])
    return res

def SPLINE_3D(u, x_cur, y_cur, z_cur, x, y, z):
    print u.shape
    print (x.shape, y.shape, z.shape)

    u_new1 = np.zeros((len(x), len(y_cur), len(z_cur)))
    print u_new1.shape
    for i in range(u_new1.shape[1]):
        for j in range(u_new1.shape[2]):
            # print (i, j)
            u_new1[:, i, j] = LAGRANGE1(u[:, i, j], x_cur, x)
    print "x done"

    u_new2 = np.zeros((len(x), len(y), len(z_cur)))
    print u_new2.shape
    for i in range(u_new2.shape[0]):
        for j in range(u_new2.shape[2]):
            # print (i, j)
            u_new2[i, :, j] = LAGRANGE1(u_new1[i, :, j], y_cur, y)
    print "y done"

    u_new3 = np.zeros((len(x), len(y), len(z)))
    print u_new3.shape
    for i in range(u_new3.shape[0]):
        for j in range(u_new3.shape[1]):
            print (i, j)
            u_new3[i, j, :] = LAGRANGE1(u_new2[i, j, :], z_cur, z)
    print "z done"
    return u_new3
################################################################################
################################################################################


def BIG_INTERPOL(u, arx, ary, arz, x, y, z):
    print np.shape(u)
    print np.shape(arx)
    print np.shape(ary)
    print np.shape(arz)
    print np.shape(x)
    print np.shape(y)
    print np.shape(z)
    result = np.zeros((len(x), len(y), len(z)))
    iterations = 0
    oldp = 0.0
    # F = RGI(u, )
    for i in range(len(x)):
        for j in range(len(y)):
            for k in range(len(z)):
                iterations += 1
                # result[i][j][k] = lagrange_interpol(u, arx, ary, arz, x[i], y[j], z[k])
                # result[i][j][k] = bliner_interpol(u, arx, ary, arz, x[i], y[j], z[k])
                result[i][j][k] = bcub_interpol(u, arx, ary, arz, x[i], y[j], z[k])
                newp = 100.0*iterations/(len(x)*len(y)*len(z))
                if newp - oldp > 0.001:
                    sys.stdout.write("\r                                     " +
                    "\r{}% {}".format(round(newp, 3), i))
                    sys.stdout.flush()
                    oldp = newp

    return result
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
######################## Reading from uw #######################################
################################################################################

ufile = Nio.open_file("uw2004.nc", "r")
vfile = Nio.open_file("vw2004.nc", "r")
u = ufile.variables["uwnd"]
v = vfile.variables["vwnd"]
# print v
lat = vfile.variables["lat"][30:0:-1]
lon = vfile.variables["lon"][60:95]
print len(lat)
time = vfile.variables["time"]
time = A(list(map(int, time[:]-time[0])))
min_t = np.min(time)/6 + 248*4
max_t = min_t + 7*4
# print (min_t, max_t)
# new_time =
time = np.linspace(min_t, max_t, 28)
u_arr = u[min_t:max_t, 0, 30:0:-1, 60:95]*u.scale_factor + u.add_offset
v_arr = v[min_t:max_t, 0, 30:0:-1, 60:95]*v.scale_factor + v.add_offset
uar = u_arr[0, ::, ::]#*u.scale_factor + u.add_offset
var = v_arr[0, ::, ::]#*v.scale_factor + v.add_offset

# new_time = np.linspace(min_t+4, (min_t+4) + 121.0/6, 2)
new_lat = np.linspace(38.0, 65.0, 250)
new_lon = np.linspace(115.0, 172.0, 250)

new_time = np.linspace(min_t+4, (min_t+4) + 121.0/6, 40)
# new_lat = np.linspace(58.251405, 60.420895, 10)
# new_lon = np.linspace(142.927289, 146.443603, 12)

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
################################################################################
######################## Reading from interpoled #######################################
################################################################################

# int_file = Nio.open_file("wind2004.nc", "r")
# uarr = int_file.variables["u"][0,:,:]
# varr = int_file.variables["v"][0,:,:]
# lat = int_file.variables["lat"][:]
# lon = int_file.variables["lon"][:]

uarr = np.ones(np.shape(uar))
varr = np.ones(np.shape(var))
for i in range(len(uar)):
    for j in range(len(uar[0])):
        uarr[i][j] = uar[i][j]
        varr[i][j] = var[i][j]
u_new = SPLINE_3D(u_arr, time, lat, lon, new_time, new_lat, new_lon)
v_new = SPLINE_3D(v_arr, time, lat, lon, new_time, new_lat, new_lon)
# print u_new




wks_type = "ps"
resources = Ngl.Resources()
resources.mpLimitMode = "LatLon"
# resources.mpMinLonF   = float(142.927289)
# resources.mpMaxLonF   = float(146.443603)
# resources.mpMinLatF   = float(58.251405)
# resources.mpMaxLatF   = float(60.420895)

# resources.mpMinLonF   = np.min(lon)
# resources.mpMaxLonF   = np.max(lon)
# resources.mpMinLatF   = np.min(lat)
# resources.mpMaxLatF   = np.max(lat)

cmap = np.array([[1.00, 0.00, 0.00], [1.00, 0.00, 0.40], \
                [1.00, 0.00, 0.80], [1.00, 0.20, 1.00], \
                [1.00, 0.60, 1.00], [0.60, 0.80, 1.00], \
                [0.20, 0.80, 1.00], [0.20, 0.80, 0.60], \
                [0.20, 0.80, 0.00], [0.20, 0.40, 0.00], \
                [0.20, 0.45, 0.40], [0.20, 0.40, 0.80], \
                [0.60, 0.40, 0.80], [0.60, 0.80, 0.80], \
                [0.60, 0.80, 0.40], [1.00, 0.60, 0.80]],'f')
rlist = Ngl.Resources()
rlist.wkColorMap         = "BlueRedGray"     # Set colors for contours.

resources.mpMinLonF   = np.min(new_lon)
resources.mpMaxLonF   = np.max(new_lon)
resources.mpMinLatF   = np.min(new_lat)
resources.mpMaxLatF   = np.max(new_lat)
resources.mpPerimOn   = True

resources.vcMonoLineArrowColor      =  True
resources.vcMinFracLengthF = 0.33
resources.vcRefMagnitudeF  = 20.0
resources.vcRefLengthF     = 0.045
resources.vcGlyphStyle              = "CurlyVector"
resources.vcLineArrowThicknessF     =   6.0

resources.mpFillOn              = True         # Turn on map fill.
resources.mpFillAreaSpecifiers  = ["Water","Land"]
resources.mpSpecifiedFillColors = [0,90]
resources.vcMonoLineArrowColor  = False
resources.mpAreaMaskingOn       = True # Draw vectors in color.

resources.tiMainString  = "wind at ohotsk sea - September 2004"
resources.tiXAxisString = "longitude"
resources.tiYAxisString = "latitude"
resources.mpPerimOn             = True            # Turn on a perimeter.
resources.mpGridMaskMode        = "MaskLand"      # Mask grid over land.
# resources.vfXCStartV  = float(lon[0])             # Define X/Y axes range
# resources.vfXCEndV    = float(lon[len(lon[:])-1]) # for vector plot.
# resources.vfYCStartV  = float(lat[0])
# resources.vfYCEndV    = float(lat[len(lat[:])-1])

#59.251405, 142.927289
#59.420895, 143.443603
# # resources.cnFillDrawOrder       = "Predraw"
resources.vfXCStartV  = float(new_lon[0])             # Define X/Y axes range
resources.vfXCEndV    = float(new_lon[len(new_lon[:])-1]) # for vector plot.
resources.vfYCStartV  = float(new_lat[0])
resources.vfYCEndV    = float(new_lat[len(new_lat[:])-1])

resources.mpGridMaskMode            = "MaskNotOcean"                #-- draw grid over ocean, not land
resources.mpGridLineDashPattern     =   2                           #-- grid dash pattern
resources.pmLabelBarDisplayMode     = "Always"                      #-- turn on a labelbar
resources.lbOrientation             = "Horizontal"                  #-- labelbar orientation
resources.lbLabelFontHeightF        =  0.008                        #-- labelbar label font size
resources.lbBoxMinorExtentF         =  0.22                         #-- decrease height of labelbar boxes
resources.lbTitleString             = "TEMPERATURE (~S~o~N~F)"      #-- labelbar title string
resources.lbTitleFontHeightF        =  0.010                        #-- labelbar title font size
resources.lbBoxMinorExtentF         =  0.18                         #-- decrease height of labelbar boxes
# resources.lbFillColor         =  False                         #-- decrease height of labelbar boxes


i = 3
# vc = Ngl.vector_map(wks,u_new,v_new,resources)
# for i in range(len(u_new)):
wks = Ngl.open_wks(wks_type,"./some.ps".format(i))
Ngl.set_values(wks,rlist)
vc = Ngl.vector_map(wks,u_new[i],v_new[i],resources)
Ngl.destroy(wks)
Ngl.end()
# for i in u_arr:
    # for j in i:
    #     for k in j:
    #         k = float(k)

# print type(u_arr[0][0][0])

################################################################################
################################################################################

# u_new = BIG_INTERPOL(u_arr, time, lat, lon, new_time, new_lat, new_lon)
# v_new = BIG_INTERPOL(v_arr, time, lat, lon, new_time, new_lat, new_lon)
# exit()
# return 0
# for i in range(n)
nc = Nio.open_file('windl2004_000.nc','w')
nc.create_dimension('x',len(new_lon))
nc.create_dimension('y',len(new_lat))
nc.create_dimension('time',len(new_time))
nc.create_variable('u', 'd', ('time', 'y', 'x'))
nc.create_variable('v', 'd', ('time', 'y', 'x'))
nc.create_variable('lat', 'd', ('y',))
nc.create_variable('lon', 'd', ('x',))
# nc.create_variable('time', 'd', ('time',))
print np.shape(u_new)
print len(new_lon)
print len(new_lat)
print np.max(np.abs(u_arr))
print nc.variables['u'].shape
nc.variables['u'][:,:,:] = u_new[:,:,:]
nc.variables['v'][:,:,:] = v_new[:,:,:]
nc.variables['lat'][:] = new_lat[:]
nc.variables['lon'][:] = new_lon[:]
# nc.variables['time'] = new_time
nc.close()

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
