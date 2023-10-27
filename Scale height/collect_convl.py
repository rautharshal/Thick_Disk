#!/usr/bin/env python
'''
To collect the convolved and reshampled cubes.
'''
import os
import sys
import numpy as np
import astropy.io.fits as p

if len(sys.argv) != 3:
	print '<any_sample_file> <outfile>'
	sys.exit()

infile = sys.argv[1]
outfile = sys.argv[2]

hdr = p.getheader(infile)

label = hdr['label']
size = hdr['nproc']
dchan = hdr['dchan']
nch = hdr['nch']
dimx = hdr['dimx']
dimy = hdr['dimy']
dimz = hdr['dimz']
factor = hdr['factor']
rdim = hdr['rdim']

idxch = []
for i in range(nch):
	idxch.append(i)

ncomb = len(idxch)
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
res = np.zeros((tcomb, rdim, rdim)) # the cube
for i in range(size):
	tout = label + str(i) + '.fits'
	tdata = p.getdata(tout)
	res[sp[i]:ep[i]] = tdata

rcube = res[0:ncomb]

hdr.set('marker', 'Convld reshampled cube')
primaryhdu = p.PrimaryHDU(data=rcube, header=hdr)
hdulist = p.HDUList([primaryhdu])
cmd = 'rm ' + outfile
os.system(cmd)
hdulist.writeto(outfile)
