#!/usr/bin/env python
'''
To modify the original sigma profile such as to correct for 
the actual sigma profile.
'''
import os
import sys
import numpy as np
import  matplotlib.pyplot as plt
from scipy import interpolate
from pylab import rc

if len(sys.argv) != 5:
	print '<in_org_sig_file> <sig_file_last_iter> <sig_file_this_itr> <prefix>'
	sys.exit()

sig_org_file = sys.argv[1]
sig_last_file = sys.argv[2]
sig_this_file = sys.argv[3]
prefix = sys.argv[4]

rc('font', family='serif', size=20)

data_org = np.load(sig_org_file)
data_last = np.load(sig_last_file)
data_this = np.load(sig_this_file)

r_org = data_org[0]
sig_org = data_org[1] # the observed MOMNT2/SG-sigma (M/S)

r_last = data_last[0]
sig_last = data_last[1] # m/s

r_this = data_this[0] 
sig_this = data_this[1] # m/s

func_sig = interpolate.interp1d(r_this, sig_this, kind='linear', bounds_error=False, fill_value='extrapolate')
sig_th_org = func_sig(r_org)

func_sig = interpolate.interp1d(r_last, sig_last, kind='linear', bounds_error=False, fill_value='extrapolate')
sig_last_org = func_sig(r_org)
"""
diff =.3*np.abs((sig_th_org - sig_org))
#diff = diff/np.max(diff)
#sig_new = (sig_last_org - diff*10000.) # m/s
sig_new = sig_last_org - diff
"""
sig_new =[]
diff=[]
for i in range(len(r_org)):
	if (r_org[i]<3000):
		sig_new.append(sig_last_org[i])
		diff.append(0)
	elif(3000<r_org[i]<10000):
		if (sig_org[i]>sig_this[i]):
			#diff.append(.35*(sig_org[i]-sig_this[i]))
			#sig_new.append(sig_last_org[i]+diff[i])
			sig_new.append(sig_last_org[i])
			diff.append(0)
		else:
			diff.append(.35*(sig_this[i]-sig_org[i]))
			sig_new.append(sig_last_org[i]-diff[i])
	else:
		sig_new.append(sig_last_org[i])
		diff.append(0)

outfile = prefix+'_sig_new'
tot = [r_org, sig_new]
np.save(outfile, tot)


sig_new_n = np.array(sig_new)

"""
plt.figure(figsize=(6,6))
plt.plot(r_org/1000., sig_org/1000., 'ro',label='$M2_{org}$')
plt.plot(r_this/1000., sig_this/1000., marker='s', color='none',mec='b', label='$M2_{now}$')
plt.plot(r_last/1000., sig_last/1000., marker='^', color='none',mec='k', label='$Sig_{last}$')
plt.plot(r_org/1000., sig_new_n/1000, 'mx', label='$Sig_{new}$')
plt.legend(frameon=False)
plt.xlabel('R (kpc)')
plt.ylabel('$\sigma_{HI} \ (km/s)$')
outfig = prefix + '_plot.pdf'
plt.title(prefix)
plt.savefig(outfig, bbox_inches='tight')
"""
#pylab.show()


#func_sig_ag = interpolate.interp1d(r_ag, disp_ag, kind='linear', bounds_error=False, fill_value='extrapolate')
