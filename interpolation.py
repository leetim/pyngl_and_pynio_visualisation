import sys
import numpy as np
import Nio
import Ngl
import os
import threading as th
from scipy.interpolate import RegularGridInterpolator as RGI
from scipy.interpolate import interp1d

A = np.array
T = np.transpose


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
            # print (i, j)
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


print "Interpolation is imported"
