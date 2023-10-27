import sys

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy import wcs

file = sys.argv[1]
file1 = sys.argv[2]
hdu = fits.open(file)
data = hdu[0].data
data = data
crpix = [350, 350]  # Original CRPIX values  # 232.1  350
crval = [0, 0]  # Original CRVAL values
cdelt = [50, -50]  # Original CDELT values  # 75.4 50
w1 = wcs.WCS(naxis=2)
w1.wcs.crpix = crpix
w1.wcs.crval = crval
w1.wcs.cdelt = cdelt
w1.wcs.ctype = ["LINEAR", "LINEAR"]
w1.wcs.cunit = ["parsec", "parsec"]

data_cut=[]
fig, ax = plt.subplots(1, 4, sharex=False, sharey=False, figsize=(24, 6))
i =0
for j in [2500,4500,6500,8500]:
    x, y = w1.all_world2pix(j,13000, 0)  # gives pixel
    x=int(x)
    data_cut.append(data[:,x]) # cut a vertical line
    Y = np.arange(data_cut[0].shape[0])

    wcs_coords = [w1.wcs_pix2world(x, y, 1) for y  in Y]  # converting cut pixel back to radius
    radius = [wcs_coords[i][1] for i in Y]
    ax[i].set_title(j)
    ax[i].plot(radius,data_cut[i],label ='model')
    ax[i].set_xlabel('Radius')
    ax[i].set_ylabel('$M_\odot\: pc^{-2}$')
    ax[i].legend()
    i = i+1
#plt.show()

hdu = fits.open(file1)
data1 = hdu[0].data
data1 = data1
crpix = [34, 33]  # Original CRPIX values
crval = [0, 0]  # Original CRVAL values
cdelt = [339.291, -339.291]  # Original CDELT values
w1 = wcs.WCS(naxis=2)
w1.wcs.crpix = crpix
w1.wcs.crval = crval
w1.wcs.cdelt = cdelt
w1.wcs.ctype = ["LINEAR", "LINEAR"]
w1.wcs.cunit = ["parsec", "parsec"]

data_cut1=[]
data_cut2 =[]
min_data=[]
norm_data=[]
i =0
for j in [2500,4500,6500,8500]:
    x, y = w1.all_world2pix(j,13000, 0)  # gives pixel
    x=int(x)
    print(x)
    data_cut1.append(data1[:,x]) # cut a vertical line
    data_cut2 = [np.where(arr == 0, np.nan, arr) for arr in data_cut1]
    min_data.append(np.nanmin(data_cut2[i]))
    norm_data.append(data_cut2[i]-min_data[i])

    #print(norm_data)
    #print(data_cut[0].shape[0])
    Y = np.arange(data_cut1[0].shape[0])

    wcs_coords = [w1.wcs_pix2world(x, y, 1) for y  in Y]  # converting cut pixel back to radius
    radius = [wcs_coords[i][1] for i in Y]

    #plt.plot(radius,data_cut1[i],label =j)
    ax[i].plot(radius, norm_data[i], label='obs')

    #plt.xlabel('Radius')
    #plt.ylabel('$M_\odot\: pc^{-2}$')
    ax[i].legend()
    i = i+1
plt.show()
