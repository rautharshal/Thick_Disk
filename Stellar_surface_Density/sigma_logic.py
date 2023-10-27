import sys
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from lmfit import Model
from matplotlib import pyplot as plt
import pandas as pd
def SIG(a,b):

    #Take the file of galaxy name,co-ordinates and fits file
    file=a
    i = b
    # Extracting the data from the file
    data = np.genfromtxt(file,delimiter=" ",dtype='str')
    G_name = data[:,0]
    F_file = data[:,1]
    RA_list = data[:,2]
    DEC_list= data[:,3]
    c=(SkyCoord(RA_list,DEC_list, unit="deg"))


    # Function for fitting
    def Gauss(x,amp,mu,sigma):
        return amp*np.exp(-(x-mu)**2/(2*sigma**2))

    # Cropping the fits file

    hdu = fits.open(F_file[i])
    img = hdu[0].data
    df = pd.DataFrame(img)  ##for spitzer
    df = df.replace(np.nan, 0, regex=True)  ## spitzer
    img = df.to_numpy()
    img = np.squeeze(img)
    med = np.median(img)
    img -= med  # Noise subtraction
    image_wcs = WCS(hdu[0].header)

    # Extracting noise from the image
    fimg = img.flatten()
    wr = np.where(fimg < 0.0)
    timg = fimg[wr]

    # Producing Histogram of the noise
    h,b =np.histogram(timg,bins=10,density=True)
    db = b[1]-b[0]
    b = b[:-1]+db/2

    # lmfit to the noise histogram
    mod = Model(Gauss)
    result=mod.fit(h,x=b,amp=1,mu=1,sigma=1)
    #result.plot()
    #plt.show()
    sig=result.best_values['sigma']
    return sig
    #print(sig)
    """

    # Finding the best crop of the image
    x1=img.flatten().std()
    crop = 50
    a=[]
    n = 5
    size = (crop, crop)
    d = image_wcs.world_to_pixel(c[i])

    while x1 > np.abs(sig):
        size = (crop,crop)
        cutout = Cutout2D(img, d, size,wcs=image_wcs)

        new_cutout=cutout.data[0:n,crop-n:crop]
        new_cutout1=cutout.data[0:n,:]
        new_cutout2=cutout.data[crop-n,:]
        new_cutout3=cutout.data[:,0:n]
        new_cutout4=cutout.data[:,crop-n]
        x2=new_cutout.std()
        x2_1=new_cutout1.std()
        x2_2 = new_cutout2.std()
        x2_3 = new_cutout3.std()
        x2_4 = new_cutout4.std()
        #print(sig,x2_1,x2_2,x2_3,x2_4,crop)
        if (x2_1<np.abs(sig) and x2_2<np.abs(sig) and x2_3<np.abs(sig) and x2_4<np.abs(sig)):
            plt.imshow(cutout.data, origin='lower', vmin=2e-5, vmax=1, cmap='gray')
            plt.show()
            break
        crop=crop+10
        Fcrop=crop
    # Adding extra 20 by 20 for precaution and writing newfits file
    Final_size=(Fcrop+20,Fcrop+20)

    cutout1 = Cutout2D(img, d, Final_size, wcs=image_wcs)
    """
