import numpy as np
import matplotlib.pyplot as plt 
from astropy.io import fits
from mgefit.find_galaxy import find_galaxy
from mgefit.sectors_photometry import sectors_photometry
from mgefit.mge_fit_sectors import mge_fit_sectors
from mgefit.mge_print_contours import mge_print_contours
import pandas as pd
import os
import shutil
import sys
from  sigma_logic import SIG  # for minimum noise
from astropy.wcs import WCS



file = sys.argv[1]
file1 = sys.argv[2]
file2 = sys.argv[3]

data = np.genfromtxt(file1,delimiter=" ",dtype='str')
F_file = data[:,1]

char='spitzer'
if F_file[int(file2)].find(char) ==-1:
    # Running for r_band
    print('In r_band')
    I=int(file2)
    #os.system('python sigma_logic.py '+file1+' '+I)
    ml=np.abs(float(SIG(file1,I)))

    print(file)
    hdu = fits.open(file)
    header=hdu[0].header
    img = hdu[0].data
    img -= 0.003 # Noise subtraction
    img[np.isnan(img)] = 0
    scale =0.39598
    print(scale,ml)
    sigma = 1.3/scale
    image_wcs = WCS(header)
    ngauss = 20


    f = find_galaxy(img,plot=1,nblob=1,level=ml) # change fraction to shift the contour

    s = sectors_photometry(img, f.eps, f.theta, f.xpeak, f.ypeak,n_sectors=19,plot=1,minlevel=ml)

    #m = mge_fit_sectors(s.radius, s.angle, s.counts, f.eps, ngauss=ngauss, plot=1,sigmapsf=sigma,qbounds=[1-f.eps+.06,1],outer_slope=2)
    m = mge_fit_sectors(s.radius, s.angle, s.counts, f.eps, ngauss=ngauss, plot=1, sigmapsf=sigma,bulge_disk=True,qbounds=[1-f.eps+.06,1], outer_slope=2)

    plt.clf()

    """
     Cropping function
    """
    n = 90
    img = img[f.xpeak-n:f.xpeak+n, f.ypeak-n:f.ypeak+n]
    xc, yc = n-f.xpeak + f.xmed, n -f.ypeak+ f.ymed

    #plt.figure(1)
    mge_print_contours(img, f.theta, xc, yc, m.sol,scale=scale,sigmapsf=sigma)

    plt.show()
    f_sigma=str(int(5*m.sol[1][-1]))
    t_gauss = str(len(m.sol[0]))

    #Adding Galaxy name to mgefitfile

    data = np.genfromtxt(file1,delimiter=" ",dtype='str')
    G_name = data[:,0]
    F_name = data[:,1]
    i_index=int(file2)
    e=G_name[i_index]
    E=e+'_r'


    # Storing the mgefit parameters in a text file
    with open("mge_file_" + e + "_rband.txt", mode="w+") as f:
        for j in range(len(m.sol)):
            content = ' '.join(map(str, m.sol[j]))
            f.write("%s\n" % content)
    f.close()

    # Calling the surface Density program
    print(f_sigma,t_gauss)
    os.system('python sden_gen_v1.py' ' mge_file_' + e + '_rband.txt '+f_sigma+' '+t_gauss+ ' 0.2 ' + E+' '+file1+' '+file2)# need to take 5 sigma

    # Making a mgefit folder and storing all galaxy mgefit text files

    directory = 'mge_fit_file'
    if (os.path.exists('mge_fit_file')):
        src = os.getcwd() + "\mge_file_" + e + "_rband.txt"
        dst = os.getcwd() + "\mge_fit_file\mge_file_" + e + "_rband.txt"
        shutil.copyfile(src, dst)
        os.remove(src)
    else:
        os.mkdir(directory)
        src = os.getcwd() + "\mge_file_" + e + "_rband.txt"
        dst = os.getcwd() + "\mge_fit_file\mge_file_" + e + "_rband.txt"
        shutil.copyfile(src, dst)
        os.remove(src)
else:
    # Running for Spitzer
    I = int(file2)
    #os.system('python sigma_logic.py ' + file1+' '+I)
    print('In s')
    ml = np.abs(float(SIG(file1, I)))
    hdu = fits.open(file)
    header = hdu[0].header
    img = hdu[0].data
    scale = .6  # to convert from pix to arcsec
    sigma = 1.6/scale  # psf in pix
    image_wcs = WCS(header)
    ngauss = 20
    df = pd.DataFrame(img)  ##for spitzer
    df = df.replace(np.nan, 0, regex=True)  ## spitzer
    img = df.to_numpy()  ## spitzer

    f = find_galaxy(img, plot=1, nblob=1, level=ml)  # change fraction to shift the contour

    s = sectors_photometry(img, f.eps, f.theta, f.xpeak, f.ypeak, n_sectors=19, plot=1, minlevel=ml)

    #m = mge_fit_sectors(s.radius, s.angle, s.counts, f.eps, ngauss=ngauss, plot=1, sigmapsf=sigma,qbounds=[1-f.eps+.06,1],outer_slope=2)
    m = mge_fit_sectors(s.radius, s.angle, s.counts, f.eps, ngauss=ngauss, plot=1, sigmapsf=sigma,bulge_disk=True,qbounds=[1 - f.eps + .06, 1], outer_slope=2)
    plt.clf()
    print(ml)
    print('The inclination is ',np.arccos(1-f.eps))
    print('position angle',f.pa)

    """
     Cropping function
    """
    n = 60
    img = img[f.xpeak - n:f.xpeak + n, f.ypeak - n:f.ypeak + n]
    xc, yc = n - f.xpeak + f.xmed, n - f.ypeak + f.ymed

    plt.figure(1)
    mge_print_contours(img, f.theta, xc, yc, m.sol, scale=.6, sigmapsf=sigma)

    plt.show()
    f_sigma = str(int(5 * m.sol[1][-1]))
    t_gauss = str(len(m.sol))

    # Adding Galaxy name to mgefitfile

    data = np.genfromtxt(file1, delimiter=" ", dtype='str')
    G_name = data[:, 0]
    F_name = data[:, 1]
    i_index = int(file2)
    e = G_name[i_index]
    E = e+'_s'

    # Storing the mgefit parameters in a text file
    with open("mge_file_" + e + "_3.6um.txt", mode="w+") as f:
        for j in range(len(m.sol)):
            content = ' '.join(map(str, m.sol[j]))
            f.write("%s\n" % content)
    f.close()

    # Calling the surface Density program

    os.system('python sden_gen_v1.py' ' mge_file_' + e + '_3.6um.txt ' + f_sigma + ' ' + t_gauss + ' 0.2 ' + E+' '+file1+' '+file2)  # need to take 5 sigma

    # Making a mgefit folder and storing all galaxy mgefit text files

    directory = 'mge_fit_file'
    if (os.path.exists('mge_fit_file')):
        src = os.getcwd() + "\mge_file_" + e + "_3.6um.txt"
        dst = os.getcwd() + "\mge_fit_file\mge_file_" + e + "_3.6um.txt"
        shutil.copyfile(src, dst)
        os.remove(src)
    else:
        os.mkdir(directory)
        src = os.getcwd() + "\mge_file_" + e + "_3.6um.txt"
        dst = os.getcwd() + "\mge_fit_file\mge_file_" + e + "_3.6um.txt"
        shutil.copyfile(src, dst)
        os.remove(src)

