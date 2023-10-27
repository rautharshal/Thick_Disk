#!/usr/bin/env python
'''
Generate the data cube of HI

To minimize loading, we do not incline the whole velocity
slice or plane, instead, we only incline 1/4 th of the plane
and use the symmetry criteria to fill the rest of it.

As the memory of Jeeves is quite high, instead of doing one
velocity slice at a time, we implement it like all velocity 
slice at the same time.
'''
import os
import sys
import numpy as np
import astropy.io.fits as p
from astropy.io import fits
from scipy import interpolate
from astropy.convolution import convolve_fft
import matplotlib
from mpi4py import MPI
#matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
#import pylab
#import resource
from scipy.interpolate import griddata
#import gc

if len(sys.argv) != 14:
	print '<in_density_file> <in_sig_file> <dist> <pix_res(pc)> <dimension> <fwhm1 (arcsec)> <fwhm2 (arcsec)> <bpa> <v_sys(km/s)> <vres (km/s)> <incl (deg)> <nproc> <out_prefix>'
	sys.exit()

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

infile = sys.argv[1]
sig_file = sys.argv[2]
dist = float(sys.argv[3])
pixres = float(sys.argv[4])
dim = int(sys.argv[5])
fwhm1 = float(sys.argv[6])
fwhm2 = float(sys.argv[7])
bpa = float(sys.argv[8])
v_sys = float(sys.argv[9])*1000. # m/s
vres = float(sys.argv[10])*1000. # m/s
incl = float(sys.argv[11]) # degree
nproc = int(sys.argv[12])
prefix = sys.argv[13]

rho_data = p.getdata(infile)
hdr = p.getheader(infile)


crpix = hdr['crpix1'] #v gets data from sden fits file
crval = hdr['crval1']
cdelt = hdr['cdelt1']

n_r = rho_data.shape[0]
radii = (np.arange(n_r) - crpix)*cdelt + crval

fwhm1_pc = (np.pi/180.)*(fwhm1/3600.)*dist*1e6 # in pc 
fwhm2_pc = (np.pi/180.)*(fwhm2/3600.)*dist*1e6 # in pc 
beam_area = fwhm1_pc*fwhm2_pc*1.13 # A Gaussian beam

r_rho = radii*1
phi = bpa
# ------------------------------

# generating the velocity profile
# For ngc 551
bf_vmax = 193.121 # km/s +/- 0.1
bf_rmax = 11392 # pc +/- 100
bf_n = .25 # +/- 0.07

rv = (bf_vmax*(radii/bf_rmax))/( (1./3.) + (2.0/3.0)*(radii/bf_rmax)**bf_n)**(1.5/bf_n) # km/s
rv = rv*1000. # m/s
# ------------------------------
# Generating dispersion profile

data_sig_ag = np.load(sig_file)
r_ag = data_sig_ag[0] # pc
disp_ag = data_sig_ag[1] # m/s
sigag = np.nanmax(disp_ag) # m/s

func_sigag = interpolate.interp1d(r_ag, disp_ag, kind='linear', bounds_error=False, fill_value='extrapolate')
#sig_ag = func_sig_ag(radii)

# Generating cube
# PACKING THE SUBROUTINES FOR RUNNING IN PARALLEL
# ------------------------------
#    Finding the kernal parameter

'''
ba1,bb1,bp1 -- Output beam parameter
ba2,bb2,bp2 -- Input beam parameter
ba3,bb3,bp3 -- Kernel beam parameter
'''
def kernpar(ba1,bb1,bp1,ba2,bb2,bp2):

	d0 = ba1**2 - bb1**2
	d2 = ba2**2 - bb2**2
	d1 = max(0.0,d0**2+d2**2-2.0*d0*d2*np.cos(2.0*(bp1-bp2)))
	d1 = np.sqrt(d1)
	bp3 = 0.0
	
	a = d0*np.sin(2.0*bp1) - d2*np.sin(2.0*bp2)
	b = d0*np.cos(2.0*bp1) - d2*np.cos(2.0*bp2)

	if ( a != 0.0 ) or (b != 0.0) :
		bp3 = 0.5*np.arctan2(a,b)
	a = ba1**2 + bb1**2 - ba2**2 - bb2**2
	ba3 = max(0.0,0.5*(a+d1))
	ba3 = np.sqrt(ba3)
	bb3 = max(0.0,0.5*(a-d1))
	bb3 = np.sqrt(bb3)

	if bb3 < 0.0 :
		bp3 = bp3 + np.pi
	return ba3,bb3,bp3
# ------------------------------
def lfill(x0,y0,x):
	f = interpolate.interp1d(x0, y0,kind='linear',bounds_error=False,fill_value=0.0)
	return f(x)
# ------------------------------
def fill(x0,y0,x):
	f = interpolate.interp1d(x0, y0,kind='cubic',bounds_error=False,fill_value=0.0)
	return f(x)
# ------------------------------

def plane(rho_plane,h_plane,th,a_r_pix,tr):
#	global d_pix

#	vertical direction
	nr = rho_plane.shape[0]

#	results = pprocess.Map(limit=nproc,reuse=1)
#	pfill = results.manage(pprocess.MakeReusable(fill))

	tplane = np.zeros((nr,th.size))
	print 'Interpolating in vertical direction.'
	for j in range(nr): # for all r
		sys.stdout.write('\rRunning %d out of %d'%(j+1,nr))
		sys.stdout.flush()
		a_h_pix = h_plane[j,:]/d_pix_x
		y0 = rho_plane[j,:]
		tplane[j,:] = lfill(a_h_pix,y0,th)

#	print 'Interpolating in vertical direction.'
#	for j in range(nr):
#		sys.stdout.write('\rRunning %d out of %d'%(i+1,nr))
#		sys.stdout.flush()
#		tplane[j,:] = results[j]

#	horizontal direction

#	results = pprocess.Map(limit=nproc,reuse=1)
#	pfill = results.manage(pprocess.MakeReusable(fill))
	
	iplane = np.zeros((tr.size,th.size))
	print 'Interpolating in radial direction.'
	for j in range(th.size):
		sys.stdout.write('\rRunning %d out of %d'%(j+1,th.size))
		sys.stdout.flush()
		y0 = tplane[:,j]
		iplane[:,j] = fill(a_r_pix,y0,tr)

#	print
#	for i in range(th.size):
#		sys.stdout.write('\rRunning %d out of %d'%(i+1,th.size))
#		sys.stdout.flush()
#		iplane[:,i] = results[i]
	print

	return iplane
# ------------------------------

def allot2(y,z,cy,cz):
	r_r = int(np.sqrt((y-cy)**2. + (z-cz)**2 ))+1
	if r_r < viplane.shape[0] :
		spec = viplane[r_r,:]
	else:
		spec = np.zeros(viplane.shape[1])
	return spec
# ------------------------------
def tconvl(image,kernel):
	conv = convolve_fft(image,kernel,interpolate_nan=True,normalize_kernel=True)
	ttt = conv.flatten()
	rslice = griddata((yyy,zzz),ttt,(ryy,rzz), method='nearest')
	return rslice
# ------------------------------

def incline_rh(i,j):
	tx = np.arange(dimx)
	ty = np.repeat(yy[i,j],dimy)
	tz = np.repeat(zz[i,j],dimz)

	xp = (tx-float(cx))*np.cos(ra) - (tz-float(cz))*np.sin(ra) + cx
	zp = (tx-float(cx))*np.sin(ra) + (tz-float(cz))*np.cos(ra) + cz
	yp = ty*1

#	Converting to r,h unit
	h_incl = xp*1
	r_incl = np.sqrt((yp-cy)**2 + (zp-cz)**2)

	h_incl = h_incl.astype(int)
	r_incl = r_incl.astype(int)

	wr =  np.where( (0 <= h_incl) & ( h_incl < dimx) & ( 0 <= r_incl) & (r_incl < dimy/2.))
	flx = (np.nansum(viplane[r_incl[wr],h_incl[wr]]))*d_pix_x/5.5 # K*km/s

	return flx
# ------------------------------
def incline_rh_cube(i,j,vspec):
	tx = np.arange(dimx)
	ty = np.repeat(yy[i,j],dimy)
	tz = np.repeat(zz[i,j],dimz)

	xp = (tx-float(cx))*np.cos(ra) - (tz-float(cz))*np.sin(ra) + cx
	zp = (tx-float(cx))*np.sin(ra) + (tz-float(cz))*np.cos(ra) + cz
	yp = ty*1

#	Converting to r,h unit
	h_incl = xp*1
	r_incl = np.sqrt((yp-cy)**2 + (zp-cz)**2)

	r_incl_pc = r_incl/d_pix_x # in parsec
	rv_incl = (bf_vmax*(r_incl_pc/bf_rmax))/( (1./3.) + (2.0/3.0)*(r_incl_pc/bf_rmax)**bf_n)**(1.5/bf_n) # km/s
	rv_incl = rv_incl + v_sys # km/s

#	rv_incl = rv_incl*1000. # m/s

	h_incl = h_incl.astype(int)
	r_incl = r_incl.astype(int)

	wr =  np.where( (0 <= h_incl) & ( h_incl < dimx) & ( 0 <= r_incl) & (r_incl < dimy/2.))
	int_incl = (viplane[r_incl[wr],h_incl[wr]])*d_pix_x/5.5 # K*km/s
	amp_incl = int_incl/(sigmg*np.sqrt(2*np.pi)) # K
	v_incl = rv_incl[wr] # km/s

	spec = np.nansum(amp_incl*np.exp(-(((vspec-v_incl)**2.)/(2.*sigmg**2.))))

#	spec = np.zeros(nch)
#	for k in range(nch):
#		spec[k] = np.nansum(amp_incl*np.exp(-(((v[k]-v_incl)**2.)/(2.*sigmg**2.))))

#	flx = (np.nansum(viplane[r_incl[wr],h_incl[wr]]))*d_pix_x/5.5 # K*km/s

	return spec
# ------------------------------
def make_slice(vspec):

	y = np.arange(dimy)
	z = np.arange(dimz)
	yy, zz = np.meshgrid(y,z)

	tslice = np.zeros((dimy,dimz))
	for i in range(dimy):
#		sys.stdout.write('\rRunning %d out of %d'%(i+1,dimy))
#		sys.stdout.flush()
		for j in range(dimz/2):

			tx = np.arange(dimx)
			ty = np.repeat(yy[j,i],dimy)
			tz = np.repeat(zz[j,i],dimz)

			xp = (tx-float(cx))*np.cos(ra) - (tz-float(cz))*np.sin(ra) + cx
			zp = (tx-float(cx))*np.sin(ra) + (tz-float(cz))*np.cos(ra) + cz
			yp = ty*1

#			Converting to r,h unit
			h_incl = xp*1
			r_incl = np.sqrt((yp-cy)**2 + (zp-cz)**2)

# 			Calculating systemic velocity
			tty = yp + 0.5 - cy
			ttz = zp + 0.5 - cz

			sgn = np.sign(5.0*np.sign(tty)+np.sign(ttz))
			psi = np.arctan((zp-cz)/(yp-cy+0.000001)) # in radian

			r_incl_pc = r_incl*d_pix_x # in parsec
			rv_incl = ((bf_vmax*(r_incl_pc/bf_rmax))/( (1./3.) + (2.0/3.0)*(r_incl_pc/bf_rmax)**bf_n)**(1.5/bf_n))*np.cos(psi)*np.sin(incl)*sgn # km/s
			rv_incl = rv_incl + v_sys # km/s

#	rv_incl = rv_incl*1000. # m/s

			h_incl = h_incl.astype(int)
			r_incl = r_incl.astype(int)

			wr =  np.where( (0 <= h_incl) & ( h_incl < dimx) & ( 0 <= r_incl) & (r_incl < dimy/2.))
			int_incl = (viplane[r_incl[wr],h_incl[wr]])*d_pix_x/5.5 # K*km/s
			amp_incl = int_incl/(sigmg*np.sqrt(2*np.pi)) # K
			v_incl = rv_incl[wr] # km/s

			spec = np.nansum(amp_incl*np.exp(-(((vspec-v_incl)**2.)/(2.*sigmg**2.))))
			tslice[j,i] = spec
			tslice[dimz-j-1,i] = spec
#			tslice[j,dimz-i-1] = spec
#			tslice[dimy-i-1,dimz-j-1] = spec
	tcslice = convl(tslice,kernel) # convlved
	
	yyy = yy.flatten()
	zzz = zz.flatten()
	ttt = tcslice.flatten()

	ry = np.arange(0.,dimy, factor)
	rz = np.arange(0.,dimz, factor)
	ryy, rzz = np.meshgrid(ry, rz)
	rslice = griddata((yyy,zzz),ttt,(ryy,rzz), method='nearest')

#	del tcslice, yy,zz,yyy,zzz,ryy,rzz	

#	spec = np.zeros(nch)
#	for k in range(nch):
#		spec[k] = np.nansum(amp_incl*np.exp(-(((v[k]-v_incl)**2.)/(2.*sigmg**2.))))

#	flx = (np.nansum(viplane[r_incl[wr],h_incl[wr]]))*d_pix_x/5.5 # K*km/s

	return rslice
# ------------------------------

def free_test():
	return np.arange(nch)


def make_vslice(i,j):

#	y = np.arange(dimy)
#	z = np.arange(dimz)
#	yy, zz = np.meshgrid(y,z)
#	tx = np.arange(dimx)

	ty = np.repeat(yy[j,i],dimy)
	tz = np.repeat(zz[j,i],dimz)

	xp = (tx-float(cx))*np.cos(ra) - (tz-float(cz))*np.sin(ra) + cx
	zp = (tx-float(cx))*np.sin(ra) + (tz-float(cz))*np.cos(ra) + cz
	yp = ty*1

#	Converting to r,h unit
	h_incl = xp*1
	r_incl = np.sqrt((yp-cy)**2 + (zp-cz)**2)

# 	Calculating systemic velocity
	tty = yp + 0.5 - cy
	ttz = zp + 0.5 - cz

	sgn = np.sign(5.0*np.sign(tty)+np.sign(ttz))
	psi = np.arctan((zp-cz)/(yp-cy+0.000001)) # in radian

	r_incl_pc = r_incl*d_pix_x # in parsec
	rv_incl = ((bf_vmax*(r_incl_pc/bf_rmax))/( (1./3.) + (2.0/3.0)*(r_incl_pc/bf_rmax)**bf_n)**(1.5/bf_n))*np.cos(psi)*np.sin(incl)*sgn # km/s

	rv_incl = rv_incl*1000. # m/s
	rv_incl = rv_incl + v_sys # m/s

	h_incl = h_incl.astype(int)
	r_incl = r_incl.astype(int)

	wr =  np.where( (0 <= h_incl) & ( h_incl < dimx) & ( 0 <= r_incl) & (r_incl < dimy/2.))
	tflx = viplane_ag[r_incl[wr],h_incl[wr]] # M_sun/pc^3
	"""
	#--------------------Radio----------------
	tflx = tflx*d_pix_x**3 # M_sun
	tflx = tflx/(1e5*2.356*dist**2.0) # Jy*km/s
	tflx = tflx*1000. # Jy*m/s
	md = tflx/(sig_ag*np.sqrt(2.0*np.pi)) # Jy, sigag unit m/s
	amp_incl = md*(beam_area/(d_pix_x**2)) # Jy/Bm
	"""

	#-----------------------3.6um-----------------
	tflx = tflx*d_pix_x**3 # M_sun
	#tflx = tflx *(1/373.3017) # MJy/str
	
	"""
	#--------------------R band-----------
	m_l = .585
	scale = .39598
	d =72
	Dl = 74.5  # Luminosity Dist in Mpc
	t_sflx = tflx *d_pix_x # M_sun/pc^2
	#t_slx = t_sflx/m_l # L_sun/pc^2
	#t_sflx_pix = t_sflx *((scale*10**6*d)/206265)**2 # L_sun/pix^2
	#t_count = (t_sflx_pix*3.823)/(4*np.pi*(Dl*3.086)**2*4.86*3.631) # nanomagy/pix^2
	amp_incl = t_sflx
	"""
	
	tsigag = func_sigag(r_incl_pc)
	sig_ag = tsigag[wr] # m/s

	amp_incl = tflx*1000./(sig_ag*np.sqrt(2.0*np.pi)) # consider f = 1, so M_sun= f*Mjy* km/s  is divided by sig_ag which is in m/s
	v_incl = rv_incl[wr] # m/s
	
	spec = np.zeros(nch)
	for k in range(nch):
		spec[k] = np.nansum(amp_incl*np.exp(-(((v[k]-v_incl)**2.)/(2.*sig_ag**2.))))
	
	return spec
# ------------------------------


n_r = r_rho.size
itr = rho_data.shape[2]

beam_min = min(fwhm1,fwhm2)
beam_max = max(fwhm1,fwhm2)

d_pix_x = pixres # parsec
d_pix_y = pixres # parsec
d_pix_z = pixres # parsec

pix_x = (d_pix_x/(dist*1e6))*180.*3600/np.pi
pix_y = (d_pix_y/(dist*1e6))*180.*3600/np.pi
pix_z = (d_pix_z/(dist*1e6))*180.*3600/np.pi

dimx = dim
dimy = dim
dimz = dim
cx,cy,cz = dimx/2.0,dimy/2.0,dimz/2.0 # center of the galaxy

r_max = np.max(r_rho)
r_min = np.min(r_rho)

a_r = np.linspace(r_min,r_max,n_r)
a_r_pix = a_r/d_pix_y

beam_max_pc = (beam_max*np.pi/(180.*3600.))*dist*1e6 # pc
reres = beam_max_pc/10. # observing resolution in pc
#factor = reres/d_pix_y
factor = 1

#	rdimx = int((dimx*d_pix_x)/reres)
#	rdimy = rdimx*1
#	rdimz = rdimx*1

#	rcx = rdimx/2
#	rcy = rdimy/2
#	rcz = rdimz/2

#	Preparing kernel
ba1 = (beam_max/pix_y)/2.35482 # sigma
bb1 = (beam_min/pix_y)/2.35482 # sigma
bp1 = np.pi/2.0 - phi

ba2 = 1.0/2.35482
bb2 = 1.0/2.35482
bp2 = 0.0

sig1,sig2,bpa = kernpar(ba1,bb1,bp1,ba2,bb2,bp2)

gx=np.arange(0,2.0*dimz,1.0)
gy=gx*1.0
gxx,gyy=np.meshgrid(gx,gy)

loc1=float(dimz)
loc2=loc1

theta = (np.pi/2.0 - bpa) # PA of the output beam, bpa is in radian
a=(np.cos(theta)**2.0/(2.0*sig1**2))+(np.sin(theta)**2.0/(2.0*sig2**2))
b=-(np.sin(2*theta)/(4.0*sig1**2))+(np.sin(2*theta)/(4.0*sig2**2))
c=(np.sin(theta)**2/(2.0*sig1**2))+(np.cos(theta)**2/(2.0*sig2**2))
kernel = np.exp(-((a*(gxx-loc1)**2.0) + 2.0*b*(gxx-loc1)*(gyy-loc2) + c*(gyy-loc2)**2.0) )

# 	------------------------------
# 	Making the plane

h_plane = rho_data[:,0,:] # pc
h_plane[:,0]=0.0
rho_plane = rho_data[:,1,:] # M_sun/pc^3

incl = incl*np.pi/180. # in radian
ra = -incl

tr = np.linspace(0.0,dimz/2-1,dimz/2)
th = np.arange(dimx)
th = abs(th - cx) # pixel unit

#print "rho ",rho_plane,"height",h_plane,"height pixel",th,"radius to pixel",a_r_pix,"radius pixel",tr

viplane_ag = plane(rho_plane,h_plane,th,a_r_pix,tr)
# viplane ready: Now we have to incline and get the spectrum

dimx = dim
dimy = dim
dimz = dim
cx, cy, cz = dimx/2., dimy/2., dimz/2.

# Defining spectral cube
v_l = v_sys - np.max(rv) - 10.*sigag # m/s
v_u = v_sys + np.max(rv) + 10.*sigag # m/s
v = np.arange(v_l,v_u,vres)
nch = v.size

y = np.arange(dimy)
z = np.arange(dimz)
yy, zz = np.meshgrid(y,z)
tx = np.arange(dimx)

idxy = []
idxz = []

for i in range(dimy):
	for j in range(dimz/2):
		idxy.append(i)	
		idxz.append(j)

ncomb = len(idxy)
mitr = ncomb/size + 1
scount = np.repeat(mitr, size)
extra = mitr*size - ncomb

sp = np.empty(size)
ep = np.empty(size)
sp[0] = 0
ep[0] = scount[0]
for i in range(1,size,1):
	sp[i] = sp[i-1] + scount[i-1]
	ep[i] = ep[i-1] + scount[i-1]

sp = sp.astype(int)
ep = ep.astype(int)

dspec = ep[0]-sp[0]
tspec = np.zeros((dspec,nch))

result = []
for i in range(sp[rank], ep[rank], 1):
	if i < ncomb:
#		sys.stdout.write('\rRunning %d out of %d'%(i+1,ncomb))
#		sys.stdout.flush()
		tspec[i-sp[rank]] = make_vslice(idxy[i],idxz[i])

#		result.append(make_vslice(idxy[i],idxz[i]))
#		result.append(free_test())
	else:
		tspec[i-sp[rank]] = np.zeros(nch)

print 'Done computation for proc %d, communicating.'%rank
comm.Barrier()
# Gathering data
# Writing to output fits file
hdr = p.Header()
outlabel = prefix + '_cube'
hdr.set('marker', 'Cubelets')
hdr.set('label', outlabel)
hdr.set('dspec', dspec)
hdr.set('nch', nch)
hdr.set('dimx', dimx)
hdr.set('dimy', dimy)
hdr.set('dimz', dimz)
hdr.set('factor', factor)
hdr.set('nproc', size)
hdr.set('sp', sp[rank])
hdr.set('ep', ep[rank])
hdr.set('crpix1',cy/factor,'reference pix y dir') # /factor was present
hdr.set('crpix2',cz/factor,'reference pix z dir') # /factor was present
hdr.set('crpix3',0.0,'reference pix vel dir')
hdr.set('crval1', 0.0)
hdr.set('crval2', 0.0)
hdr.set('crval3', v[0])
hdr.set('cdelt1',d_pix_y*factor,'in pc') # *factor was present
hdr.set('cdelt2',d_pix_y*factor,'in pc') # *factor was present
hdr.set('cdelt3',vres,'in km/s')
hdr.set('bunit','K','in pc')
#hdr.set('datamax',np.nanmax(rcube))
#hdr.set('datamin',np.nanmin(rcube))
hdr.set('incl',incl*180./np.pi, 'degree')
hdr.set('bmaj',fwhm1/3600.,'degree')
hdr.set('bmin',fwhm2/3600.,'degree')
hdr.set('bpa',bpa,'degree')
hdr.set('dist',dist,'Mpc')
hdr.set('pixres',pixres,'actual pix res')
hdr.set('infile',infile,'input file')

primaryhdu = p.PrimaryHDU(data=tspec, header=hdr)
hdulist = p.HDUList([primaryhdu])
outfile = prefix + '_cube'+ str(rank)+'.fits'
cmd = 'rm ' + outfile
os.system(cmd)
hdulist.writeto(outfile)

# Writing kernel
if rank == 0:
	hdr.set('marker', 'Kernel')
	primaryhdu = p.PrimaryHDU(data=kernel, header=hdr)
	hdulist = p.HDUList([primaryhdu])
	outfile = prefix + '_kernel.fits'
	cmd = 'rm ' + outfile
	os.system(cmd)
	hdulist.writeto(outfile)
