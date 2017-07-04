from PIL import Image, ImageDraw
from interpolation import *
# draw.point
from scipy.misc import imsave
import numpy as np
import subprocess
from sys import argv
# print argv
res_name = argv[4]
# exit()
A = np.array
def get_im_arr(image):
    pix = image.load()
    z = np.ones((image.size[1], image.size[0], 3), dtype=np.int16)
    # print pix[1, 1]
    # print z.shape
    for i in range(image.size[0]):
        for j in range(image.size[1]):
            z[j][i][0] = pix[i, j][0]
            z[j][i][1] = pix[i, j][1]
            z[j][i][2] = pix[i, j][2]
            # print z[i, j]
            # print pix[i, j]
    return z

image = Image.open(argv[1])
image1 = Image.open(argv[2])
image2 = Image.open(argv[3])
pix = get_im_arr(image)
pix1 = get_im_arr(image1)
pix2 = get_im_arr(image2)

# print pix
arr_x = np.linspace(0.0, 100.0, pix2.shape[0])
arr_y = np.linspace(0.0, 100.0, pix2.shape[1])
# print arr_y.shape
# print arr_x.shape
# print pix1.shape
# exit()
# new_x = np.linspace(0.0, 100.0, n_height)
# new_y = np.linspace(0.0, 100.0, width)
# new_z = np.zeros((width, n_height, 3))
k1 = (1.0 - float(pix2.shape[0])/pix2.shape[1])*100.0
new_x = np.linspace(0.0, 100.0, pix.shape[0])
new_y = np.linspace(k1/2, 100.0-k1/2, pix.shape[0])
# print (k1/2, 100.0-k1/2)
# print pix2.shape
new_z1 = pix1
# new_z2 = pix2
# print pix.shape
# print pix1.shape
print pix2.shape
# new_z1 = np.zeros((len(new_x), len(new_y), 3))
# new_z2 = np.zeros((len(new_x), len(new_y), 3))
# for i in range(3):
#     new_z2[:,:,i] = SPLINE_2D(pix2[:,:,i], arr_x, arr_y, new_x, new_y)
new_z2 = pix2[:,pix.shape[0]/2:pix.shape[0]/2+pix.shape[0],:]
print new_z2.shape
# k1 = (new_z2.shape[1]-1024)/2
# k2 = k1 + 1024
# new_z2 = new_z2[:,k1:k2,:]
# np.round(new_z1)
# for i in range(3):
#     new_z2[:,:,i] = SPLINE_2D(pix2[:,:,i], arr_x, arr_y, new_x, new_y)
# np.round(new_z2)
del pix1
del pix2
print (pix.shape[0], pix.shape[1] + new_z1.shape[1], new_z2.shape[1], 3)
res = np.zeros((max(pix.shape[0], new_z1.shape[0], new_z2.shape[0]), pix.shape[1] + new_z1.shape[1] + new_z2.shape[1], 3), dtype=np.int16)
# print res.shape
# print pix.shape
for i in range(new_z1.shape[0]):
    for j in range(new_z1.shape[1]):
        res[i][j] = new_z1[i][j]
    for j in range(new_z2.shape[1]):
        if new_z2.shape[0] <= i:
            break
        res[i][new_z1.shape[1]+j] = new_z2[i][j]
    for j in range(pix.shape[1]):
        res[i][new_z1.shape[1] + new_z2.shape[1] + j] = pix[i][j]
# print pix1
imsave(res_name, res)
# print res.shape

# arr = 256.0*np.abs(np.sin(arr))+0.1
# print arr
# print z[0]
# imsave("old_some.png", pix1)
# subprocess.call('convert -delay 5 +repage zalip/*.png 5.gif', shell=True)
