#!/usr/bin/env python
'''
To plot the sig and the MOM2 values after convergence.
'''
import os
import sys
import numpy as np
import pylab
import astropy.io.fits as p
from pylab import rc

if len(sys.argv) != 7:
    print '<mom_extract_out> <m2_original> <sig_last_file> <optimization_l(pc)> <optimization_r(pc)> <prefix>'
    sys.exit()

infile_now = sys.argv[1]
infile_last = sys.argv[2]
sig_last = sys.argv[3]
opt_l = float(sys.argv[4])/1000. # kpc
opt_r = float(sys.argv[5])/1000. # kpc
prefix = sys.argv[6]

#rc('text', usetex=True)
rc('font', family='serif', size=20)

data1 = np.load(infile_now)
data2 = np.load(infile_last)
data3 = np.load(sig_last)

r_now = data1[0]/1000. # kpc
m2_now = data1[1]/1000. # km/s

r_org = data2[0]/1000. # kpc
m2_org = data2[1]/1000. # km/s
#m2_err_org = data2[2]/1000. # km/s

r_sig = data3[0]/1000. # kpc
sig_sig = data3[1]/1000. # km/s

wr = np.where(r_org > opt_r)
r_org[wr] = np.nan
m2_org[wr] = np.nan
#m2_err_org[wr] = np.nan

wr = np.where(r_now > opt_r)
r_now[wr] = np.nan
m2_now[wr] = np.nan

wr = np.where(r_sig > opt_r)
r_sig[wr] = np.nan
sig_sig[wr] = np.nan


pylab.figure(figsize=(6,6))
pylab.plot(r_org, m2_org, marker='o', color='r', ls='',ms=7, mec='r', mew=2, label=r'$\rm MOM2_{obs}$')
#pylab.errorbar(r_org, m2_org, yerr=m2_err_org, color='r', ls='', capsize=0)
pylab.plot(r_now, m2_now, marker='^', color='none', mec='b', mew=2, ms=7, ls='', label=r'$\rm MOM2$')
pylab.plot(r_sig, sig_sig, marker=(5,1), color='none', mec='k', mew=1, ms=10, ls='', label=r'$\rm \sigma$')
pylab.axvline(x=opt_l, ls='--', lw=1., color='k')
#pylab.axvline(x=opt_r, ls='--', lw=1., color='k')
pylab.xlabel(r'R (kpc)')
pylab.ylabel(r'$\rm MOM2/ \sigma \ (km \thinspace s^{-1})$')
pylab.legend(loc=1, frameon=False, prop={'size':14}, ncol=1, handletextpad=0.1)
outfile = prefix + '_sigplot.pdf'
pylab.savefig(outfile, bbox_inches='tight')
pylab.show()
