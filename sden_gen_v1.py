#!/usr/bin/env python
'''
To generate a face-on galaxy with some exponential/Gaussian
distribution.
'''
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pylab
from scipy.interpolate import griddata
import shutil
if len(sys.argv) != 8:
    print('<infile> <dim> <Ngauss_to_consider> <intrinsic_q> <prefix> <main file> <index>')
    print(sys.argv)
    sys.exit()

infile = sys.argv[1]
dim = int(sys.argv[2])
ngauss = int(sys.argv[3])
aq = float(sys.argv[4])  # intrinsic axial ratio # .2
prefix = sys.argv[5]


# ------------------------------
# Reading the input file

data = np.genfromtxt(infile)
L = data[0, :]
sig = data[1, :]
q = data[2, :]

cx = dim / 2.
cy = dim / 2.
# Calculating pseudo inclination
incl = np.arccos(np.sqrt((q ** 2 - aq ** 2) / (1 - aq ** 2)))
print('The angle is ',incl)# inclinations of all the Gaussians in radian
#print(incl)
#sys.exit()
# The deprojection of the surface densities
# has to be done with different inclination

# Constructing the surface density map and preparing for deprojection

x = np.arange(dim)
y = np.arange(dim)
xx, yy = np.meshgrid(x, y)  # Actual projected grid

xx = xx - cx
yy = yy - cy

xn = np.arange(dim)
yn = np.arange(dim)
xxn, yyn = np.meshgrid(xn, yn)
xxn = xxn - cx
yyn = yyn - cy

sden = np.zeros((len(xx),len(xx)))
zz = np.zeros((len(xx),len(xx)))
for i in range(ngauss):
    x = np.arange(dim)
    y = np.arange(dim)
    xx, yy = np.meshgrid(x, y)  # Actual projected grid

    xx = xx - cx
    yy = yy - cy
    xn = np.arange(dim)
    yn = np.arange(dim)
    xxn, yyn = np.meshgrid(xn, yn)
    xxn = xxn - cx
    yyn = yyn - cy

    tsden = (L[i] / (2 * np.pi * (sig[i] ** 2) * q[i])) * (
        np.exp(-((1 / (2 * (sig[i] ** 2))) * ((xx) ** 2 + (yy) ** 2 / q[i] ** 2)))) # eq 33 of capillari 2020

    # Now deprojection has to be done
    yy = yy / np.cos(incl[i])
    xxx = xx.flatten()
    yyy = yy.flatten()
    z = tsden
    #print('Light of individual gaussian',np.sum(tsden),'component',L[i],sig[i],q[i])
    zzz = tsden * np.cos(incl[i])

    zz = zz+z
    sden_incl = griddata((xxx, yyy), zzz.flatten(), (xxn, yyn), method='linear', fill_value=0.0) # interpolate the data zzz from the meshgrid on the point xxn
    sden = sden + sden_incl
    #print('Light of indiidual deprojected gaussian', np.sum(zzz), np.sum(sden_incl))
"""
height, width = sden.shape
center_x, center_y = width // 2, height // 2

# Step 3 and Step 4: Calculate the radial profile
max_radius = min(center_x, center_y)
radii = np.arange(max_radius + 1)
profile = []

for radius in radii:
    y, x = np.ogrid[-center_y:height - center_y, -center_x:width - center_x]
    mask = x ** 2 + y ** 2 <= radius ** 2
    pixels_in_ring = np.sum(mask)
    sum_pixel_values = np.sum(sden[mask])
    average_pixel_value = sum_pixel_values / pixels_in_ring
    profile.append(average_pixel_value)
print(profile)
sys.exit()
tot = [radii,profile]
"""
# Extracting the surface density profile
rad = np.sqrt(xxn ** 2 + yyn ** 2)

tdim = int(dim / 2) - 1
rad_bin = np.arange(tdim)

prof = np.zeros(tdim)
for i in range(tdim - 1):
    wr = np.where((rad >= rad_bin[i]) & (rad < rad_bin[i + 1]))
    prof[i] = np.mean(sden[wr])

tot = [rad_bin, prof]

outfile = prefix + '_sden_prof_before'
np.save(outfile, tot)

# Making a numpy file folder and storing all galaxy .npy files
directory = 'npy_files'
if (os.path.exists('npy_files')):
    src = os.getcwd() + "\\"+outfile+".npy"
    dst = os.getcwd() + "/npy_files\\"+outfile+".npy"
    shutil.copyfile(src, dst)
    os.remove(src)
else:
    os.mkdir(directory)
    src = os.getcwd() + "\\"+outfile+".npy"
    dst = os.getcwd() + "/npy_files\\"+outfile+".npy"
    shutil.copyfile(src, dst)
    os.remove(src)


pylab.figure()
pylab.plot(rad_bin, prof, 'r-')
pylab.xlim(0, cx)
pylab.ylim(0, np.max(prof))
pylab.figure()
pylab.imshow(sden, origin='lower')
pylab.colorbar()
pylab.plot(cx, cy, marker='o', ms=10, color='r')
pylab.contour(sden, levels=10, colors='w', linewidths=2)
pylab.show()

main_file = sys.argv[6]
index = sys.argv[7]
os.system('python sden_unit_conversion.py '+main_file+' '+index)
