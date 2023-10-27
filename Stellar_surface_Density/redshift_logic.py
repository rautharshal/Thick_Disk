import sys
from astropy.wcs import WCS
from astroquery.ipac.ned import Ned
import numpy as np
from astropy.nddata import Cutout2D
from astropy.io import fits
from astropy import units as u
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import redshift_calc as rc
import pandas as pd
#Take the file of galaxy name,co-ordinates and fits file
file=sys.argv[1]

# Extracting the data from the file
data = np.genfromtxt(file,delimiter=" ",dtype='str')
data1 =np.genfromtxt(file,delimiter=" ",dtype='float')
G_name = data[:,0]
F_file = data[:,1]
RA_list = data[:,2]
DEC_list= data[:,3]
lum_dist = data1[:,4]
c=(SkyCoord(RA_list,DEC_list, unit="deg"))
f_data=[]
arcsec=[]
# Finding redshifts

for galaxy in range(len(G_name)):
    #r_z=result_table['Published Redshift']
    #b=min(r_z)
    #dist=(b*300000)/65
    #rc.z=b
    dist1=lum_dist[galaxy]
    theta=.05/dist1
    arcsec.append(((3600*180)/np.pi)*theta)

# Cropping the fits file
for i in range(len(G_name)):

    hdu = fits.open(F_file[i])
    hdu.info()
    img = hdu[0].data
    df = pd.DataFrame(img)  ##for spitzer
    df = df.replace(np.nan, 0, regex=True)  ## spitzer
    img = df.to_numpy()
    img = np.squeeze(img)
    med = np.median(img)
    img -= med  # Noise subtraction

    print(img)
    image_wcs = WCS(hdu[0].header)
    size = u.Quantity((arcsec[i]), u.arcsec)
    d = image_wcs.world_to_pixel(c[i])
    print(size,d)
    cutout = Cutout2D(img, d, size,wcs=image_wcs)
    print(cutout)
    f_data.append(cutout.data)
    plt.imshow(cutout.data, origin='lower', vmin=2e-5, vmax=1, cmap='gray')
    plt.show()