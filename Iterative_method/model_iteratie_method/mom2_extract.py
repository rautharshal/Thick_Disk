#!/usr/bin/env python
'''
To extract the MOM2 profile from the MOM2 image produced 
by the momgen code.
'''
import os
import sys
import numpy as np
import astropy.io.fits as p
import pylab

if len(sys.argv) != 5:
	print '<in_mom2> <in_obs_sig_file> <inner_radius_to_blank> <prefix>'
	sys.exit()

inmom2 = sys.argv[1]
insig = sys.argv[2]
r_blnk = float(sys.argv[3])
prefix = sys.argv[4]

mom2 = p.getdata(inmom2)
hdr = p.getheader(inmom2)

data_sig = np.load(insig)

dimx = mom2.shape[0]
dimy = mom2.shape[1]

cx = hdr['crpix1']
cy = hdr['crpix2']
crval1 = hdr['crval1']
crval2 = hdr['crval2']
cdelt1 = hdr['cdelt1']
cdelt2 = hdr['cdelt2']
incl = hdr['incl']
incl = incl*np.pi/180.
# ------------------------------

r_ag = data_sig[0]
sig_ag = data_sig[1]

dr = r_ag[1]-r_ag[0] 
r_bound = r_ag - dr/2.
a_pc = np.concatenate([r_bound, [(r_ag[-1]+dr/2.)]])
b_pc = a_pc*np.cos(incl)
a_axis = a_pc/cdelt1  # in pixel
b_axis = b_pc/cdelt1 

nr = a_pc.size-1
# ------------------------------

tx = np.arange(dimx)
ty = np.arange(dimy)

xx, yy = np.meshgrid(tx, ty)

# ------------------------------

# blanking the central region 
wr = np.where(abs((xx-cx)*cdelt1+crval1) < r_blnk)
mom2[wr] = np.nan
# ------------------------------

vdisp = np.zeros(nr)
vdisp_err = np.zeros(nr)
pa = 0.0
for i in range(nr):
	sys.stdout.write('\rRunning for radi %d/%d pc'%(a_pc[i],a_pc[-1]))
	sys.stdout.flush()

	tidx1 = np.zeros(xx.shape)
	tidx2 = np.zeros(yy.shape)

	aa1 = a_axis[i]
	bb1 = b_axis[i]

	aa2 = a_axis[i+1]
	bb2 = b_axis[i+1]

	elip1 = (((xx-cx)*np.cos(pa) + (yy-cy)*np.sin(pa))**2)/(aa1**2) + (((xx-cx)*np.sin(pa) - (yy-cy)*np.cos(pa))**2)/(bb1**2)
	elip2 = (((xx-cx)*np.cos(pa) + (yy-cy)*np.sin(pa))**2)/(aa2**2) + (((xx-cx)*np.sin(pa) - (yy-cy)*np.cos(pa))**2)/(bb2**2)

	wr = np.where(elip1 >= 1.)
	tidx1[wr] = 1.

	wr = np.where(elip2 >= 1.)
	tidx2[wr] = 1.

	tidx = tidx1 + tidx2

	wr = np.where(tidx == 1.)
	vdisp[i] = np.nanmean(mom2[wr])
	vdisp_err[i] = np.nanstd(mom2[wr])

outfile = prefix + '_sig'
tot = [r_ag, vdisp]
np.save(outfile, tot)


print
pylab.figure()
pylab.plot(r_ag, vdisp/1000., 'ro')
pylab.plot(r_ag, sig_ag/1000., marker='^', color='none', mec='b')
pylab.axhline(y=8., ls='--',c='k')
#pylab.ylim(4,10)
#pylab.show()

