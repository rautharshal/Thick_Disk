#!/usr/bin/env python
'''
To accumulate all the fits files and stitch them together
and generatet the final cube
'''
import os
import sys
import numpy as np
import astropy.io.fits as p
from astropy.convolution import convolve_fft
from scipy.interpolate import griddata
from mpi4py import MPI

if len(sys.argv) != 4:
	print '<Any sample file> <kernel_file> <out_label>'
	sys.exit()

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
rsize = comm.Get_size()

infile = sys.argv[1]
kfile = sys.argv[2]
prefix = sys.argv[3]

def tconvl(image,kernel):
	conv = convolve_fft(image,kernel,nan_treatment='interpolate',normalize_kernel=True, allow_huge=True)
	ttt = conv.flatten()
	rslice = griddata((yyy,zzz),ttt,(ryy,rzz), method='nearest')
	return conv
# ------------------------------

hdr = p.getheader(infile)
label = hdr['label']
size = hdr['nproc']
dspec = hdr['dspec']
nch = hdr['nch']
dimx = hdr['dimx']
dimy = hdr['dimy']
dimz = hdr['dimz']
factor = hdr['factor']

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

tcomb = ep[-1]

res = np.zeros((tcomb, nch))
for i in range(size):
	tout = label + str(i) + '.fits'
	tdata = p.getdata(tout)
	res[sp[i]:ep[i]] = tdata

fres = res[0:ncomb]
rcube = np.zeros((nch,dimz,dimx))
for i in range(ncomb):
	rcube[:,idxz[i],idxy[i]] = fres[i]
	rcube[:,dimz-idxz[i]-1,idxy[i]] = fres[i]

# Cube restored proceeding to convolution
kernel = p.getdata(kfile)

y = np.arange(dimy)
z = np.arange(dimz)
yy, zz = np.meshgrid(y,z)

yyy = yy.flatten()
zzz = zz.flatten()
ry = np.arange(0.,dimy, factor)
rz = np.arange(0.,dimz, factor)
ryy, rzz = np.meshgrid(ry, rz)

idxch = []
for i in range(nch):
	idxch.append(i)

ncomb = len(idxch)
mitr = ncomb/rsize + 1
scount = np.repeat(mitr, rsize)
extra = mitr*rsize - ncomb

sp = np.empty(rsize)
ep = np.empty(rsize)
sp[0] = 0
ep[0] = scount[0]
for i in range(1,rsize,1):
	sp[i] = sp[i-1] + scount[i-1]
	ep[i] = ep[i-1] + scount[i-1]

sp = sp.astype(int)
ep = ep.astype(int)

yyy = yy.flatten()
zzz = zz.flatten()
ry = np.arange(0.,dimy, factor)
rz = np.arange(0.,dimz, factor)
ryy, rzz = np.meshgrid(ry, rz)

print 'Convolving. Time consuming part. Be patient.'

rdim = np.arange(0.,dimz,factor).size
dchan = ep[0] - sp[0]
tcube = np.zeros((dchan, rdim, rdim))

for i in range(sp[rank], ep[rank], 1):
	if i < ncomb:
#		sys.stdout.write('\rRunning %d out of %d'%(i+1,ncomb))
#		sys.stdout.flush()
		tcube[i-sp[rank]] = tconvl(rcube[idxch[i]], kernel)
	else:
		tcube[i-sp[rank]] = np.zeros((rdim, rdim))

print 'Done convolution for proc %d, communicating.'%rank
comm.Barrier()

# Will be written different files and then will be collected by
# a serial program. The file sizes are so heavy that commumication
# overhead is much much higher in parallel processing than the 
# computation time. Hence, each processor will simply write the files
# into the disk, and later will be collected by a different code.

outlabel = prefix + '_ccube'
hdr.set('marker', 'channelets')
hdr.set('dchan', dchan)
hdr.set('label', outlabel)
hdr.set('nproc', rsize)
hdr.set('sp', sp[rank])
hdr.set('ep', ep[rank])
hdr.set('rdim',rdim)

primaryhdu = p.PrimaryHDU(data=tcube, header=hdr)
hdulist = p.HDUList([primaryhdu])
outfile = prefix + '_ccube'+ str(rank)+'.fits'
cmd = 'rm ' + outfile
os.system(cmd)
hdulist.writeto(outfile)
