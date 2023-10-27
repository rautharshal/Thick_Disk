#!/usr/bin/env python
'''
To break the loop if the MOM2 is close to the observed one.
'''
import os
import sys
import numpy as np
from scipy import interpolate


def diff(infile1, infile2, rad1, rad2, p_err):
	data1 = np.load(infile1)
	data2 = np.load(infile2)

	r1 = data1[0]   # original sigma file
	sig1 = data1[1] # m/s

	r2 = data2[0]
	sig2 = data2[1] # m/s

	wr = np.where((r1 >= rad1) & (r1 <= rad2))
	tr1 = r1[wr]
	tsig1 = sig1[wr]

	func_sig = interpolate.interp1d(r2, sig2, kind='linear', bounds_error=False, fill_value='extrapolate')
	tsig2 = func_sig(tr1)

#	For comparison:
	tdiff = abs((tsig1-tsig2)/tsig1)
	A = str(np.max(tdiff))
	print np.max(tdiff)
	with open('error.txt','a') as file:
		file.write('\n'+A)
	if np.max(tdiff) < p_err : 
		flag = 1
	else:
		flag = 0
	
	return flag
