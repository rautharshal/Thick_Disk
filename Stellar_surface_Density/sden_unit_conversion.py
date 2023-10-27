import os
import sys
import numpy as np
from astropy.io import fits
import shutil

file=sys.argv[1]
data = np.genfromtxt(file,delimiter=" ",dtype='str')
G_name = data[:,0]
F_file = data[:,1]
data1 = np.genfromtxt(file,delimiter=" ",dtype='float')

dist = data1[:,4]
M_l = data1[:,5]
dist_lum = data1[:,6]
#zp = data1[:,5]
#fo = data1[:,7]
i = int(sys.argv[2])
mu=[]
infile=G_name[i]
#B= data1[:,6]

#print(a)
hdu = fits.open(F_file[i])
#exp_time=hdu[0].header['EXPTIME']
exp_time=53.907456
char='spitzer'
if F_file[i].find(char) == -1:
    #running r_band
    data_l_pc2 = np.load('npy_files/' + infile + '_r_sden_prof_before.npy')
    scale = 0.39598
    count = np.delete(data_l_pc2[1], np.where(data_l_pc2[1] == 0))
    """
    #------------------sdss flux to Mag/pix^2-------------------
    a = a * 3.631 * 10 ** (-6) # nanomaggy to jansky
    #M = -2.5*np.log10(a)+zp[i] # jansky to Mag
    f = np.full(a.shape,fo[i])
    b = np.full(a.shape,B[i])
    print(a.shape,f.shape)
    M = (-2.5 / np.log(10))*(np.arcsinh((a/f) / (2*b)) + np.log(b)) + 5*np.log(scale)
    print(M)# jansky to mag
    """
    #--------------------------R band flux to mass---------------
    """
    #--------------Old method-------------------
    # Dl = 74.5 * 10**6 *3.08*10**16 # Luminosity distance in metre for ngc551
    Dl = 67.5 * 10 ** 6 * 3.08 * 10 ** 16  # Luminosity distance in metre for ngc2410
    a = a * 3.631 * 10 ** (-6)  # nanomaggy to Jansky
    surf_d_L = 4 * np.pi * Dl ** 2 * a * 10 ** (-26) * 4.81 * 10 ** 14  # Flux to Luminosity
    surf_d_L = surf_d_L / (3.823 * 10 ** 26)  # luminosity in L_sun
    surf_d_L = surf_d_L * (206265 / (scale * 10 ** 6 * dist)) ** 2  # L/pix^2 to L/pc^2
    """
    L = count * 4 * np.pi * (dist_lum[i] ** 2) * (3.086 ** 2) * 4.86 * 3.631 / 3.823  # From sdss nanomaggy to L_sun/pix^2
    surf_d_L = L * (206265 / (scale * 10 ** 6 * dist[i])) ** 2 # l_sun/pix^2 to L_sun/pc^2
    surf_d_M = surf_d_L * M_l[i]
    out_file = G_name[i] + '_r_sden_prof_M_pc2_final'
    radbin = data_l_pc2[0]
    b = np.delete(data_l_pc2[0], -1)
    tot = [b, surf_d_M]
    np.save(out_file, tot)
else:
    data_l_pc2 = np.load('npy_files/' + infile + '_s_sden_prof_before.npy')
    a = np.delete(data_l_pc2[1], np.where(data_l_pc2[1] == 0))
    print('not found')
   #A=a*((10**(-.03))/(.05*.05))
    A = a * 373.3017  # Mjy/str to M_sun/pc^2
    out_file = G_name[i] + '_s_sden_prof_M_pc2_final'
    radbin = data_l_pc2[0]
    b = np.delete(data_l_pc2[0], -1)
    tot = [b, A]
    np.save(out_file, tot)
#print(a)
directory = 'npy_files_sden'
if (os.path.exists('npy_files_sden')):
    src = os.getcwd() + "\\"+out_file+".npy"
    dst = os.getcwd() + "/npy_files_sden\\"+out_file+".npy"
    shutil.copyfile(src, dst)
    os.remove(src)
else:
    os.mkdir(directory)
    src = os.getcwd() + "\\"+out_file+".npy"
    dst = os.getcwd() + "/npy_files_sden\\"+out_file+".npy"
    shutil.copyfile(src, dst)
    os.remove(src)