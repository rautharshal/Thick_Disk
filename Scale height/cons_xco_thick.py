#!/usr/bin/env python
'''
To solve the differential eqn of Hydrostatic Equilibrium 
in presence of molecular Hydrogen to find best X_CO
'''
import numpy as np
import sys
import os
import dsolve_v4
import pprocess
import astropy.io.fits as p
import pylab
from scipy import interpolate

if len(sys.argv) != 17:
	print '<in_h1_den_file> <H_2 den file> <stellar_den_file> <stellar_vfile> <h2_vfile> <dist(Mpc)> <r1(pc)> <r2(pc)> <dr(pc)> <sig_h1(km/s)> <serial/parallel(s/p)> <nproc> <outden_rho_star> <out_rho_file_atomic_gas> <out_rho_file_mg> <out_sig_file>'
	sys.exit()

h1_file = sys.argv[1]
h2_file = sys.argv[2]
stellar_file = sys.argv[3]
stellar_vfile = sys.argv[4]
h2_vfile = sys.argv[5]
dist = float(sys.argv[6])
r1 = float(sys.argv[7])
r2 = float(sys.argv[8])
dr = float(sys.argv[9])
sig_ag = float(sys.argv[10]) # km/s
proc = sys.argv[11]
nproc = int(sys.argv[12])
outdenfile_s = sys.argv[13]
outdenfile_ag = sys.argv[14]
outdenfile = sys.argv[15]
outsigfile = sys.argv[16]
# ------------------------------

radii = np.arange(r1, r2, dr)
# ------------------------------
# HI surface density 
data = np.genfromtxt(h1_file) # M_sun/pc^2
r_h1 = data[:,0] # in arcsec
r_h1 = r_h1*1e3 # in pc
sden_th1 = data[:,1] # M_sun/pc^2
func_sden = interpolate.interp1d(r_h1, sden_th1, kind='linear', bounds_error=False, fill_value=0.0)
sden_h1 = func_sden(radii)
# Correcting for Helium
sden_ag = sden_h1*1.0
#sig_ag = 7.0 # km/s
# ------------------------------

# H_2 surface density
data = np.genfromtxt(h2_file)
r_h2 = data[:,0] # arcsec
r_h2 = 1e6*dist*np.pi*(r_h2/3600.)/180. # in pc
sden_th2 = data[:,1] # M_sun/pc^2
func_sden = interpolate.interp1d(r_h2, sden_th2, kind='linear', bounds_error=False, fill_value=0.0)
sden_mg = func_sden(radii) # M_sun/pc^2

#sig_mg = 5.0 # km/s

# -----------------------------

# H_2 velocity dispersion

data_h2 = np.genfromtxt(h2_vfile)
sig_th2 = data_h2[:,1]
r_sig_h2 = data_h2[:,0] # in arcsec
r_sig_h2 = 1e6*dist*np.pi*(r_sig_h2/3600.)/180. # in pc
func_sigh2 = interpolate.interp1d(r_sig_h2, sig_th2, kind='linear', bounds_error=False, fill_value=0.0)

#sig_h2 = func_sigh2(radii)
#sig_mg = sig_h2
sig_mg = 6.

#sig_mg=[6.]
#sden_th2 = data[:,1] # K*km/s
#func_sden = interpolate.interp1d(r_h2, sden_th2, kind='linear', bounds_error=False, fill_value=0.0)
#sden_h2 = func_sden(radii) # K*km/s
# Taking care of X factor
#sden_mg = sden_h2*4.04 # M_sun/pc^2
# ------------------------------

# Stellar surface density

data = np.genfromtxt(stellar_file)
#r0_s = 3.3*1e3 # in pc
l_s = 5.421 # in kpc
r_s = data[:,0] # in pc
#r_s = 1e6*dist*np.pi*(r_s/3600.)/180. # in pc
sden_ts = data[:,1] # M_sun/pc^2
func_sden = interpolate.interp1d(r_s, sden_ts, kind='linear', bounds_error=False, fill_value=0.0)
sden_s = func_sden(radii) # M_sun/pc^2

# ------------------------------------

# Stellar vvelocity dispersion
#sig_s = (60.85*np.sqrt(r0_s*sden_s))/1000.0 # km/s
#sig_s = (1.879*np.sqrt(l_s*sden_s))  # km/s


data_s = np.genfromtxt(stellar_vfile)
sig_ts = data_s[:,1]
r_sig_s = data_s[:,0] # in arcsec
r_sig_s = 1e6*dist*np.pi*(r_sig_s/3600.)/180. # in pc
func_sigs = interpolate.interp1d(r_sig_s, sig_ts, kind='linear', bounds_error=False, fill_value=0.0)
sig_s = func_sigs(radii)

# ------------------------------

'''
# H2 and stellar surface densities are to be calculated from 
# the exponential fit given Leroy et al. 2009

# Stellar surface density
mstar = 10**10.9
r0_star = 3.3*1e3 # pc
a0_star = mstar/(2*np.pi*r0_star**2) # M_sun/pc^2

sden_star = a0_star*np.exp(-radii/r0_star) # M_sun/pc^2
# ------------------------------

# Molecular surface density
mh2 = 10**9.7
r0_h2 = 3.1*1e3 # pc
a0_h2 = mh2/(2*np.pi*r0_h2**2) # M_sun/pc^2

sden_h2 = a0_h2*np.exp(-radii/r0_h2) # M_sun/pc^2
# ------------------------------
'''
# N7331 isothermal Halo property from deBlok 2008 (ISO, cons gamma, Kruppa IMF, pg: 2668 Table-4)

# ISO
#rho_0 = 2.44e-2
#rc_kpc = 5.40*1e3 # in pc

# NFW
rho_0 = 28*1e-3 # == rho_s  # for ngc551
rc_kpc = 10.214*1e3 # in pc == r_s # for ngc551

# Rotation curve

# ngc551 
bf_vmax = 193.121 # km/s +/- 0.1
bf_rmax = 11392 # pc +/- 100
bf_n = .25 # +/- 0.07

rv = (bf_vmax*(radii/bf_rmax))/( (1./3.) + (2.0/3.0)*(radii/bf_rmax)**bf_n)**(1.5/bf_n)

dvdr = (bf_vmax/bf_rmax)* ( ( (2./3.)*  (radii/bf_rmax)**bf_n + 1./3. )**(-1.5/bf_n)) * ( 0.5 - 0.5*(radii/bf_rmax)**bf_n ) / ( (radii/bf_rmax)**bf_n + 0.5 )

rf = 2.0*rv*dvdr/radii
# ------------------------------

# Start sampling radius
'''
UNITS
dist in pc
vel in km/s
density output will be in M_sun/pc^3 
'''

# ------------------------------

n = radii.size
percentage_err = 0.1

rho_mg = np.zeros((radii.size,2,500))
rho_g = np.zeros((radii.size,2,500))
rho_s = np.zeros((radii.size,2,500))
# ------------------------------

def loop(i):
	hs, rhos, hg, rhog, hmg, rhomg= dsolve_v4.control(rho_0, rc_kpc, radii[i], rf[i], sden_s[i], sden_ag[i], sden_mg[i],sig_s[i], sig_ag, sig_mg, percentage_err)
	return hs, rhos, hg, rhog, hmg, rhomg

# ------------------------------
'''
i = 40
hs, rhos, hg, rhog, hmg, rhomg = dsolve_v3.control(rho_0, rc_kpc, radii[i], rf[i], sden_s[i], sden_ag[i], sden_mg[i], sig_s[i], sig_ag, sig_mg, percentage_err)

saveout = sys.stdout
f = open('rho_den.txt', 'w')
sys.stdout = f
for i in range(500):
	print hs[i], rhos[i], hg[i], rhog[i], hmg[i], rhomg[i]
f.close()
sys.stdout = saveout

sys.exit()
'''
# Serial processing 
if proc == 's':
	for i in range(n):
		print 'Running for radii %d pc'%radii[i]
		print 'hello'
		rho_s[i,0,:], rho_s[i,1,:], rho_g[i,0,:], rho_g[i,1,:], rho_mg[i,0,:], rho_mg[i,1,:]= loop(i)
		#raw_input('Halt\n')

elif proc == 'p':

	results = pprocess.Map(limit=nproc, reuse=1)
	ploop = results.manage(pprocess.MakeReusable(loop))
	for i in range(n):
		ploop(i)

	for i in range(n):
		sys.stdout.write('\rRunning for radius %d pc'%radii[i])
		sys.stdout.flush()
		rho_s[i,0,:], rho_s[i,1,:], rho_g[i,0,:], rho_g[i,1,:], rho_mg[i,0,:], rho_mg[i,1,:]= results[i]
else:
	print 'Please specify p or s for parallel or serial processing.'
	sys.exit()
print

hdr = p.Header()
hdr.set('sigg', sig_ag,'sig of atomic gas')

hdr.set('ctype1','radius in pc')
hdr.set('ctype2','z height in pc')
hdr.set('label','Atomic Gas')
hdr.set('bunit','M_sun/pc^3 for atomic gas')
primaryhdu = p.PrimaryHDU(data=rho_g,header=hdr)
hdulist = p.HDUList([primaryhdu])
cmd = 'rm '+outdenfile_ag
os.system(cmd)
hdulist.writeto(outdenfile_ag)

hdr = p.Header()
hdr.set('sigg', sig_ag,'sig of thick comp')

hdr.set('ctype1','radius in pc')
hdr.set('crpix1',0.0,'reference pix of radius')
hdr.set('crval1',r1,'in parsec')
hdr.set('cdelt1',dr,'in pc')
hdr.set('ctype2','z height in pc')
hdr.set('label','Star')
hdr.set('bunit','M_sun/pc^3 for star')
primaryhdu = p.PrimaryHDU(data=rho_s,header=hdr)
hdulist = p.HDUList([primaryhdu])
cmd = 'rm '+outdenfile_s
os.system(cmd)
hdulist.writeto(outdenfile_s)

# Constructing sigma file

vc = np.zeros((1,2,radii.size))

vc[0,0,:] = radii[:]
vc[0,1,:] = np.repeat(sig_mg,1)


hdr = p.Header()
#hdr.set('fthin', pc_mg,'fraction of thin component')
#hdr.set('fthick', pc_mgt,'fraction of thick component')
#hdr.set('sigtn', sig_mg,'sig of thin comp')
#hdr.set('sigtk', sig_mgt,'sig of thick comp')

hdr.set('ctype1','radius in pc')
hdr.set('ctype2','H2 dispersion in km/s')
hdr.set('ctype3','z height in pc')

hdr.set('crpix1',0.0,'reference pix of radius')
hdr.set('crpix2',0.0,'ref pix of height')
hdr.set('crpix3',0.0,'ref pix of height')

hdr.set('crval1',r1,'in parsec')
hdr.set('crval3',0.0,'pc')

hdr.set('cdelt1',dr,'in pc')
hdr.set('label','Molecular Gas')

hdr.set('bunit','M_sun/pc^3 for molecular gas')
primaryhdu = p.PrimaryHDU(data=rho_mg,header=hdr)
hdulist = p.HDUList([primaryhdu])
cmd = 'rm '+outdenfile
os.system(cmd)
hdulist.writeto(outdenfile)

hdr.set('label','Thin component velocity disp')
primaryhdu = p.PrimaryHDU(data=vc, header=hdr)
hdulist = p.HDUList([primaryhdu])
cmd = 'rm '+outsigfile
os.system(cmd)
hdulist.writeto(outsigfile)
