import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import  os
import sys
from astropy.coordinates import SkyCoord
import shutil

import redshift_logic as rl


file =sys.argv[1]
data = np.genfromtxt(file,delimiter=" ",dtype='str')

G_name = data[:,0]
F_file = data[:,1]
RA_list = data[:,2]
DEC_list= data[:,3]

c=(SkyCoord(RA_list,DEC_list, unit="deg"))


for i in range(len(F_file)):

    hdu = fits.open(F_file[i])
    Final_image = rl.f_data[i]

    # Put the cutout image in the FITS HDU
    hdu[0].data = Final_image
    hdu[0].header.update(rl.cutout.wcs.to_header())
    #hdu[0].header.update(sl.cutout1.wcs.to_header())

    cutout_filename = G_name[i]+'_cutout.fits'

    hdu[0].writeto(cutout_filename, overwrite=True)
    hdu.close()

    # Making a fits folder and storing all galaxy fits files
    directory = 'dot_croped_fits'
    if (os.path.exists('dot_croped_fits')):
        src = os.getcwd() + "\\" +cutout_filename
        dst = os.getcwd() + "\dot_croped_fits\\"+cutout_filename
        shutil.copyfile(src, dst)
        print('done1')

    else:
        os.mkdir(directory)
        src = os.getcwd() + "\\" +cutout_filename
        dst = os.getcwd() + "\dot_croped_fits\\"+cutout_filename
        shutil.copyfile(src, dst)
        print('done')


    #call the mgefit program
    I=str(i)

    os.system('python mgefit_parameter.py ' +cutout_filename +' '+file +' '+I) # add loop I and see
    #os.system('python mgefit_parameter.py ' + F_file[i] +' '+ file + ' ' + I)  # add loop I and see
    #os.remove(src)