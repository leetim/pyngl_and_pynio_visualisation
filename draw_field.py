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

################################################################################
######################## Reading from uw #######################################
################################################################################

cmap = np.array([[1.00, 0.00, 0.00], [1.00, 0.00, 0.40], \
                [1.00, 0.00, 0.80], [1.00, 0.20, 1.00], \
                [1.00, 0.60, 1.00], [0.60, 0.80, 1.00], \
                [0.20, 0.80, 1.00], [0.20, 0.80, 0.60], \
                [0.20, 0.80, 0.00], [0.20, 0.40, 0.00], \
                [0.20, 0.45, 0.40], [0.20, 0.40, 0.80], \
                [0.60, 0.40, 0.80], [0.60, 0.80, 0.80], \
                [0.60, 0.80, 0.40], [1.00, 0.60, 0.80]],'f')
# rlist = Ngl.Resources()
# rlist.wkColorMap         = "BlueRedGray"     # Set colors for contours.


#convert -delay 15 -trim +repage *.png kholmsk1.gif
def print_vector_field(u, v, lat, lon, out_file = "/some.ps", wks_type = "ps", title = "Okhotsk sea"):
    wks_type = wks_type
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
    resources.mpMinLonF   = np.min(lon)
    resources.mpMaxLonF   = np.max(lon)
    resources.mpMinLatF   = np.min(lat)
    resources.mpMaxLatF   = np.max(lat)
    resources.mpPerimOn   = True

    resources.vcMonoLineArrowColor      =  True
    resources.vcMinFracLengthF = 0.33
    resources.vcRefMagnitudeF  = 20.0
    resources.vcRefLengthF     = 0.045
    resources.vcGlyphStyle              = "CurlyVector"
    # resources.vcLineArrowThicknessF     =   6.0

    resources.mpFillOn              = True         # Turn on map fill.
    resources.mpFillAreaSpecifiers  = ["Water","Land"]
    resources.mpSpecifiedFillColors = [0,90]
    resources.vcMonoLineArrowColor  = False
    resources.vcLabelsUseVectorColor = True
    resources.vcLevelColors = A([[1.0*(10 + i)/20, 0.0, 1.0*(10 - i)/10] for i in range(11)])
    resources.mpAreaMaskingOn       = True # Draw vectors in color.

    resources.tiMainString  = title
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
    resources.vfXCStartV  = float(lon[0])             # Define X/Y axes range
    resources.vfXCEndV    = float(lon[len(lon[:])-1]) # for vector plot.
    resources.vfYCStartV  = float(lat[0])
    resources.vfYCEndV    = float(lat[len(lat[:])-1])

    # resources.gsnSpreadColors = True
    resources.vcLevelSelectionMode = "ManualLevels"
    resources.vcLevelSpacingF = 28.0/9
    resources.vcMinLevelValF = 0.0
    resources.vcMaxLevelValF = 28.0

    # resources.cnLevelSelectionMode = "ManualLevels"
    # resources.cnLevelSpacingF = 28.0/11
    # resources.cnMinLevelValF = 0.0
    # resources.cnMaxLevelValF = 28.0
    #
    # resources.mpLevelSpacingF = 28.0/11
    # resources.mpMinLevelValF = 0.0
    # resources.mpMaxLevelValF = 28.0

    resources.mpGridMaskMode            = "MaskNotOcean"                #-- draw grid over ocean, not land
    resources.mpGridLineDashPattern     =   2                           #-- grid dash pattern
    resources.pmLabelBarDisplayMode     = "Always"                      #-- turn on a labelbar
    resources.lbOrientation             = "Horizontal"                  #-- labelbar orientation
    resources.lbLabelFontHeightF        =  0.008                        #-- labelbar label font size
    resources.lbBoxMinorExtentF         =  0.22                         #-- decrease height of labelbar boxes
    resources.lbTitleString             = "WIND SPEED (m/c~S~2~N~)"      #-- labelbar title string
    resources.lbTitleFontHeightF        =  0.010                        #-- labelbar title font size
    resources.lbBoxMinorExtentF         =  0.18                         #-- decrease height of labelbar boxes
    # resources.lbFillColor         =  False                         #-- decrease height of labelbar boxes


    i = 3
    # vc = Ngl.vector_map(wks,u_new,v_new,resources)
    # for i in range(len(u_new)):
    wks = Ngl.open_wks(wks_type, out_file)
    # Ngl.set_values(wks,rlist)
    vc = Ngl.vector_map(wks,u,v,resources)
    Ngl.destroy(wks)
    # Ngl.end()
# for i in u_arr:
    # for j in i:
    #     for k in j:
    #         k = float(k)

# print type(u_arr[0][0][0])
