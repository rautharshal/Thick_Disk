#!/usr/bin/env python
'''
To calculate the MOMNT maps from the data cube
'''
import os
import sys
import numpy as np
import astropy.io.fits as p
import pylab

if len(sys.argv) != 3:
	print '<incube> <prefix>'
	sys.exit()

incube = sys.argv[1]
prefix = sys.argv[2]

cube = p.getdata(incube)
hdr = p.getheader(incube)

nch = cube.shape[0]
crpix3 = hdr['crpix3']
crval3 = hdr['crval3']
cdelt3 = hdr['cdelt3']
crpix1 = hdr['crpix1']
crval1 = hdr['crval1']
cdelt1 = hdr['cdelt1']
crpix2 = hdr['crpix2']
crval2 = hdr['crval2']
cdelt2 = hdr['cdelt2']
incl = hdr['incl']
incl = incl*np.pi/180. # in radian
print crpix1,crval1,cdelt1,crpix2,crval2,cdelt2,crpix3,crval3,cdelt3
'''
dimx = cube.shape[1]
dimy = cube.shape[2]

x = (np.arange(dimx)-crpix1)*cdelt1 + crval1
y = (np.arange(dimy)-crpix2)*cdelt2 + crval2

xx, yy = np.meshgrid(x,y)
a = tradii
b = a*np.cos(incl)

el = xx**2/a**2 + yy**2/b**2
wr = np.where(el > 1.)
'''
vel = (np.arange(nch) - crpix3)*cdelt3 + crval3

m0 = np.nansum(cube, axis=0)

tcube = np.zeros(cube.shape)

for i in range(nch):
	tcube[i] = vel[i]*cube[i]

m1 = np.nansum(tcube, axis=0)
m1 = m1/m0

for i in range(nch):
	tcube[i] = ((vel[i] - m1)**2)*cube[i]

m2 = np.nansum(tcube, axis=0)/m0
m2 = np.sqrt(m2)

m0 = m0*abs(cdelt3)

#m0[wr] = np.nan
#m1[wr] = np.nan
#m2[wr] = np.nan

outfile = prefix + '_mom0.fits'
primaryhdu = p.PrimaryHDU(data=m0, header=hdr)
hdulist = p.HDUList([primaryhdu])
cmd = 'rm ' + outfile
os.system(cmd)
hdulist.writeto(outfile)

outfile = prefix + '_mom1.fits'
primaryhdu = p.PrimaryHDU(data=m1, header=hdr)
hdulist = p.HDUList([primaryhdu])
cmd = 'rm ' + outfile
os.system(cmd)
hdulist.writeto(outfile)

outfile = prefix + '_mom2.fits'
primaryhdu = p.PrimaryHDU(data=m2, header=hdr)
hdulist = p.HDUList([primaryhdu])
cmd = 'rm ' + outfile
os.system(cmd)
hdulist.writeto(outfile)

