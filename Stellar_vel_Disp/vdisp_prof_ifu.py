#!/usr/bin/env python
'''
To estimate the surface density profile from the MOM0 maps.
We also incorporate variable PA and incl.
'''
import os

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits as p
import sys
import pylab
from pylab import rc
from scipy import interpolate
from astropy.wcs import WCS
from astropy import wcs
from astropy import units as u

from matplotlib.colors import LinearSegmentedColormap



rc('font', family='serif', size=18)
if len(sys.argv) != 11:
    print('<infits_file> <cx> <cy> <incl> <pa> <start_radius(arcsec)> <end_radius> <bin_width(-1==1/3 beam)> <name> <prefix>')
    print(len(sys.argv))
    print(sys.argv[9])
    sys.exit()

infile = sys.argv[1]
wcx = float(sys.argv[2])
wcy = float(sys.argv[3])
incl = float(sys.argv[4])
pa = float(sys.argv[5])
st_rad = float(sys.argv[6]) # arcsec
en_rad = float(sys.argv[7]) # arcsec
d_rad = float(sys.argv[8]) # pc
name = sys.argv[9]
prefix = sys.argv[10]

# ------------------------------
#rc('text', usetex=True)
#rc('font', family='serif', size=18)
# ------------------------------

hdul = p.open(infile)
phdr = hdul[0].header
ehdr = hdul[1].header
vdisp = hdul[1].data[15]
mass_d = hdul[1].data[19]
flux = hdul[1].data[3]

# ------------------------------
ctype1 = phdr['ctype1']
ctype2 = phdr['ctype2']
crval1 = phdr['crval1']
crval2 = phdr['crval2']
crpix1 = phdr['crpix1']
crpix2 = phdr['crpix2']
cd1_1 = phdr['cd1_1']
cd1_2 = phdr['cd1_2']
cd2_1 = phdr['cd2_1']
cd2_2 = phdr['cd2_2']
cunit1 = phdr['cunit1']
cunit2 = phdr['cunit2']
# ------------------------------

w = wcs.WCS(naxis=2)

w.wcs.crpix = [crpix1, crpix2]
w.wcs.cdelt = np.array([cd1_1, cd2_2])
w.wcs.crval = [crval1, crval2]
w.wcs.ctype = [ctype1, ctype2]
w.wcs.cunit = [cunit1, cunit2]
#w.wcs.set_pv([(2, 1, 45.0)])

worldcrd = np.array([[wcx, wcy]], dtype=np.float64)
#worldcrd = np.array([[crval1, crval2]], dtype=np.float64)
pix = w.wcs_world2pix(worldcrd, 0)
cx, cy = pix[0] # This is in python unit 0 --> 1
# ------------------------------

# Creating new header
hdr = p.Header()
hdr.set('ctype1', ctype1)
hdr.set('ctype2', ctype2)
hdr.set('crval1', crval1)
hdr.set('crval2', crval2)
hdr.set('crpix1', crpix1)
hdr.set('crpix2', crpix2)
hdr.set('cd1_1', cd1_1)
hdr.set('cd1_2', cd1_2)
hdr.set('cd2_1', cd2_1)
hdr.set('cd2_2', cd2_2)
hdr.set('cunit1', cunit1)
hdr.set('cunit2', cunit2)

#hdr.set('ctype1', 'RA--SIN')
#hdr.set('ctype2', 'DEC--SIN')

primaryhdu = p.PrimaryHDU(data=vdisp, header=hdr)
primaryhdu.writeto('vdisp.fits', overwrite=True)

primaryhdu = p.PrimaryHDU(data=flux, header=hdr)
primaryhdu.writeto('flux.fits', overwrite=True)

mass_d = 10**mass_d
mass_d = mass_d/(((.972*10**6*72)/206265)**2)
# ------------------------------

cdelt1 = cd1_1*3600. # in arcsec
cdelt2 = cd2_2*3600. # in arcsec
#dimx = vdisp.shape[1]
#dimy = vdisp.shape[0]

dimx = mass_d.shape[1]
dimy = mass_d.shape[0]

incl = incl*(np.pi/180.) # in radian
pa = pa*(np.pi/180.) + np.pi/2.

'''
# ------------------------------
hist = hdr['history*']
# Determining the beam
nhist = len(hist)
for i in range(nhist):
	hline = str(hist[i])
	tval = hline.find('CLEAN BMAJ')
	if tval > 0 : # Phrase found
		bm1 = float(hline.split('=')[1].split()[0])*3600. # arcs
		bm2 = float(hline.split('=')[2].split()[0])*3600.
		bpa = float(hline.split('=')[3].split()[0])
		break;

#beam_pix = abs(cdelt1*cdelt2)/(1.13*bm1*bm2) # beam per pixel
#bm1_kpc = (bm1*np.pi/180.)*dist*1e3 # kpc
#bm2_kpc = (bm2*np.pi/180.)*dist*1e3 # kpc
beam_area = bm1*bm2*1.13 # arcs^2
# ------------------------------
# Reading PA and INCL from the FAT output file

fat_data = np.genfromtxt(infat_file)
fat_rad = fat_data[:,0] # in arcsec
fat_incl = fat_data[:,5] # in deg
fat_pa = fat_data[:,7] # in deg

fat_incl = fat_incl*(np.pi/180.) # in radian
fat_pa = fat_pa*(np.pi/180.) + np.pi/2.

# ------------------------------
'''
# Defining the radial bins
radii = np.arange(st_rad, en_rad, d_rad) # in arcsec
nr = radii.size-1
a_arcs = radii + 0.01
#+ (radii[1]-radii[0])/2. # in arcsec
b_arcs = a_arcs*np.cos(incl) # in arcsec

a_axis = a_arcs/abs(cdelt1) # in pixel
b_axis = b_arcs/abs(cdelt2) # in pix

# ------------------------------
'''
fat_rad = fat_rad/abs(cdelt1) # in pixel
func_incl = interpolate.interp1d(fat_rad, fat_incl, kind='linear', bounds_error=fat_incl[-1])
func_pa = interpolate.interp1d(fat_rad, fat_pa, kind='linear', bounds_error=fat_incl[-1])
# ------------------------------
'''
tx = np.arange(dimx)
ty = np.arange(dimy)

xx, yy = np.meshgrid(tx,ty)
sden = np.zeros(nr)
sden2 = np.zeros(nr)
sden3 = np.zeros(nr) # Jy*km/s/arcsec^2

for i in range(nr):
#    sys.stdout.write('\rRunning for radi %d/%d pc'%(a_pc[i],en_rad))
#    sys.stdout.flush()

    tidx1 = np.zeros(xx.shape)
    tidx2 = np.zeros(xx.shape)

    aa1 = a_axis[i]
    bb1 = b_axis[i]

    aa1 = a_axis[i] # in pixel
    tincl = incl*1
    tpa = pa*1
    bb1 = aa1*np.cos(tincl)

    aa2 = a_axis[i+1]
    bb2 = aa2*np.cos(tincl)

#    elip1 = (((xx-cx)*np.cos(pa) + (yy-cy)*np.sin(pa))**2)/(aa1**2) + (((xx-cx)*np.sin(pa) - (yy-cy)*np.cos(pa))**2)/(bb1**2)
#    elip2 = (((xx-cx)*np.cos(pa) + (yy-cy)*np.sin(pa))**2)/(aa2**2) + (((xx-cx)*np.sin(pa) - (yy-cy)*np.cos(pa))**2)/(bb2**2)

    elip1 = (((xx-cx)*np.cos(tpa) + (yy-cy)*np.sin(tpa))**2)/(aa1**2) + (((xx-cx)*np.sin(tpa) - (yy-cy)*np.cos(tpa))**2)/(bb1**2)
    elip2 = (((xx-cx)*np.cos(tpa) + (yy-cy)*np.sin(tpa))**2)/(aa2**2) + (((xx-cx)*np.sin(tpa) - (yy-cy)*np.cos(tpa))**2)/(bb2**2)


    wr = np.where(elip1 >= 1.)
    tidx1[wr] = 1.
    wr = np.where(elip2 >= 1.)
    tidx2[wr] = 1.

    tidx = tidx1 + tidx2
    wr = np.where(tidx == 1.0)

    sx = np.arange(dimx)
    ttheta = tpa*1
    sy = np.tan(ttheta)*(sx-cx) + cy
    #print(sx,sy)
    tidx[tidx==0] = np.nan
    tidx[tidx==2] = np.nan
    #    input('Halt\n')
    fig =plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111, projection=w)
    # Plot the data using WCS coordinates

    viridis = plt.cm.get_cmap('viridis')
    colors = viridis(np.linspace(0, 1, 256))
    colors[0] = [1, 1, 1, 1]  # Set color for 0 intensity to white
    cmap = LinearSegmentedColormap.from_list('custom_cmap', colors)
    im = ax.imshow(mass_d, cmap=cmap, origin='lower')
    pylab.imshow(tidx, origin='lower',alpha=.8,cmap='binary')
    pylab.plot(sx, sy, ls='--', color='r', lw=2)
    # Add a colorbar
    cbar = plt.colorbar(im,fraction=.05, pad=.04)
    ax.tick_params(axis='both', labelsize=20)
    cbar.ax.tick_params(labelsize=20)
    plt.text(1.4, .5, '$\sigma_{\star} (km s^{-1})$', va='center', rotation='vertical',fontsize=25,transform=ax.transAxes)
    pylab.xlabel('RA(J2000)',fontsize=23)
    pylab.ylabel('DEC(J2000)',fontsize=25)
    xaxis = ax.coords[0]
    xaxis.set_major_formatter('hh:mm:ss')
    xaxis.set_ticks(number=4,size =4)
#    pylab.xlim(100,150)
#    pylab.ylim(100,150)
    outfig = prefix + '_rings' + str(i+1) + '.pdf'
    #print(outfig)
    pylab.savefig(outfig, bbox_inches='tight')
    #pylab.show()
    #sys.exit()
#    input('Halt\n')

    #tdisp = vdisp[wr] # unit Jy/bm*m/s
    tdisp = mass_d[wr]
    tsden = np.nanmean(tdisp) # km/s
#    tsden = np.nanmedian(tmom0) # Jy/bm*m/s
#    sden3[i] = (tsden/beam_area)/1000. # jy/arcsec^2*km/s

#    sden2[i] = (2.343*tsden*(0.18**2))/(1.13*(bm1/3600.)*(bm2/3600.)*(np.pi**2)*1e4) # M_sun/pc^2 Starting from HI mass within a beam
#    sden2[i] = sden2[i]*np.cos(tincl) # deprojected surface density
#    sden[i] = tsden*8.85/(bm1*bm2) # M_sun/pc^2 Starting form T_B.
#    sden[i] = sden[i]*np.cos(tincl) # deprojected surface density
    sden[i] = tsden # km/s

rad_bin = a_arcs + d_rad/2. 
rad_bin = rad_bin[:-1]
dist = (rad_bin*10**6*72)/206265

outfile = prefix + '551_mass.txt'
saveout = sys.stdout
f = open(outfile, 'w')
sys.stdout = f
print('# Radius(arcs)  Velocity dispersion(km/s)')
for i in range(nr):
    print(dist[i], sden[i])
sys.stdout = saveout
f.close()


# ------------------------------
fig = pylab.figure(figsize=(6,6))
ax = fig.add_subplot(111)
pylab.plot(rad_bin, sden, marker='', ls='-', color='r', lw=2)
#pylab.plot(rad_bin, sden2, marker='', ls='--', color='b', lw=2)
pylab.text(0.65, 0.85, r'%s'%name, transform=ax.transAxes)
pylab.xlabel(r'Radius (arcsec)')
pylab.ylabel(r'$\rm \sigma_* \ (km \thinspace s^{-1})$')
outfile = prefix + '_plot.pdf'
pylab.savefig(outfile, bbox_inches='tight')
#plt.show()