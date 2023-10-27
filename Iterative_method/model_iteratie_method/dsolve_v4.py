#!/usr/bin/env python
'''
Just solve the differential equation by using input parameters and give the 
output.
INPUT:

rkpc
rot_factor
sden_s
sden_g
sigs
sigg
percentage_err
'''
import numpy as np
import sys
import os
from scipy.integrate import ode
from scipy import integrate
from matplotlib import pylab
from scipy import interpolate
import pprocess


# To calculate the dark matter density (iso-thermal halo)
def darkrho_iso(z_kpc):
    rho_h = rho_0 / (1.0 + (rkpc ** 2 + z_kpc ** 2) / (rc_kpc ** 2))
    return rho_h


# ------------------------------
# To calculate the dark matter density (NFW halo)
def darkrho_nfw(z_kpc):
    r_run = np.sqrt(rkpc ** 2 + z_kpc ** 2)
    rho_h = rho_0 / ((r_run / rc_kpc) * (1. + r_run / rc_kpc) ** 2)
    return rho_h


# ------------------------------

# Single component solver to find extremas
def sderiv(t, z):
    term1 = z[1]
    #	term2 = -(z[0]/sig**2.0)*( 0.054*(z[0] + rho_var + rho_var2 + rho_var3 + darkrho_iso(t) ) - rot_factor ) + (z[1]**2)/z[0]
    term2 = -(z[0] / sig ** 2.0) * (0.054 * (z[0] + rho_var + rho_var2 + darkrho_nfw(t)) - rot_factor) + (
                z[1] ** 2) / z[0]
    return np.array([term1, term2])


# ------------------------------
# FINDING THE MAXIMUM & MINIMUM OF THE CENTRAL DENSITY

def single(sigs, sden, percentage_err, var_h, var_rho, var_h2, var_rho2,flag, tcden):
    global rho_var, sig, rho_var2

    sig = sigs
    func_rho = interpolate.interp1d(var_h, var_rho, kind='linear', bounds_error=False, fill_value=0.0)
    func_rho2 = interpolate.interp1d(var_h2, var_rho2, kind='linear', bounds_error=False, fill_value=0.0)

    itr = 500
    h_max = 10  # in pc
    t = np.linspace(0.0001, h_max, itr)
    rho = np.zeros(itr)

    #	Calculating initial cden assuming a scale height of 1kpc for star
    if flag == 's':  # star z0 = 1kpc
        #		cden = np.max([np.max(var_rho),np.max(var_rho2),sden/(np.pi*1000.)])
        cden = sden / (np.pi * 1000.)
    elif flag == 'ag':  # atomic gas z0=1kpc
        #		cden = np.max([np.max(var_rho),np.max(var_rho2),sden/(np.pi*1000.)])
        #		cden = sden/(np.pi*5000.)
        cden = sden / (np.pi * 5000.)
    elif flag == 'mg':  # molecular gas z0 = 500pc
        #		cden = np.max([np.max(var_rho),np.max(var_rho2),sden/(np.pi*500.)])
        #		cden = sden/(np.pi*2000.)
        cden = sden / (np.pi * 2000.)
    if tcden != 0.0:
        cden = tcden

    #	solver = ode(sderiv).set_integrator('dopri5')
    solver = ode(sderiv).set_integrator('dop853')

    for j in range(1000):
        init = [cden, 0.0]
        solver.set_initial_value(init)
        prho = func_rho(t)
        prho2 = func_rho2(t)

        for i in range(itr):
            rho_var = prho[i]
            rho_var2 = prho2[i]

            rho[i] = solver.integrate(t[i])[0]
        tsden = 2.0 * integrate.trapz(rho, t)

        #		if np.mod(j,20) == 0 and flag == 's':
        #			print flag

        #			pylab.figure()
        #			pylab.plot(t, rho)
        #			pylab.show()

        #			print sden, tsden
        #			raw_input('Halt\n')

        if abs((sden - tsden) / sden) < percentage_err / 100.0:
            #			print 'converged'
            #			print 'No of iteration is %d'%j
            #			print sden, tsden, abs((sden - tsden)/sden)
            #			raw_input('Halt\n')
            break
        else:
            try:
                #				wrr = np.where(rho < np.max(rho)/2.5)
                #				tcritic = wr[0][0]
                wr = np.where(rho < np.max(rho) / 2.0)
                fwhm = 2.0 * t[wr[0][0]]
                c = fwhm / 2.3548
                wr = np.where(t > 5.0 * c)
                t_sig = t[wr[0][0]]
                darea = (sden - tsden)
                da = darea / (np.sqrt(2.0 * np.pi) * c)
                cden = cden + da
            except:
                #print 'Need to change the step size.'
                h_max = h_max + 10.0
                t = np.linspace(0.0001, h_max, itr)
    #				print 'z max set to %d pc'%h_max
    #	print flag
    #	pylab.figure()
    #	pylab.plot(t, rho)
    #	pylab.show()
    #	raw_input('Halt\n')
    return t, rho


# ------------------------------
def error(st1, sr1, st2, sr2, gt1, gr1, gt2, gr2, mgt1, mgr1, mgt2, mgr2,err_lim):
    min_err = -1.
    if st1.size == st2.size and gt1.size == gt2.size and mgt1.size == mgt2.size:
        err1 = np.max(abs((sr2 - sr1) / sr1))
        err2 = np.max(abs((gr2 - gr1) / gr1))
        err3 = np.max(abs((mgr2 - mgr1) / mgr1))
        print 'Coming here'
        print err1, err2, err3, err_lim

        min_err = np.max([err1, err2, err3])
        if err1 < err_lim and err2 < err_lim and err3 < err_lim:
            escape = True
            min_err = np.max([err1, err2, err3])
        else:
            escape = False
    else:
        escape = False
    return escape, min_err


# ------------------------------
def control(rho0, rckpc, kpc, rf, sdens, sdeng, sdenmg, ssig, gsig, mgsig, p_err):
    global rho_0, rc_kpc, rot_factor, rkpc

    rho_0 = rho0
    rc_kpc = rckpc
    rkpc = kpc
    rot_factor = rf
    sden_s = sdens
    sden_g = sdeng
    sden_mg = sdenmg
    sigs = ssig
    sigg = gsig
    sigmg = mgsig
    percentage_err = p_err

    var_h = np.linspace(0., 100., 100)
    var_rho = np.repeat(0., 100)
    var_h2 = var_h * 1
    var_rho2 = var_rho * 1

    cden_s = 0.0
    cden_ag = 0.0
    cden_mg = 0.0

    for j in range(100):
        #		Solving for star
        ths1, trhos1 = single(sigs, sden_s, percentage_err, var_h, var_rho, var_h2, var_rho2,'s',
                              cden_s)
        cden_s = np.max(trhos1)
        var_h, var_rho = 1.0 * ths1, 1.0 * trhos1

        #		Solving for HI
        thg1, trhog1 = single(sigg, sden_g, percentage_err, var_h, var_rho, var_h2, var_rho2, 'ag',
                              cden_ag)
        cden_ag = np.max(trhog1)
        var_h2, var_rho2 = 1.0 * thg1, 1.0 * trhog1

        #		Solving for H_2 thin
        thmg1, trhomg1 = single(mgsig, sden_mg, percentage_err, var_h, var_rho, var_h2, var_rho2,
                                'mg', cden_mg)
        cden_mg = np.max(trhomg1)
        var_h, var_rho = 1.0 * thmg1, 1.0 * trhomg1



        #		Second round

        #		For star
        ths2, trhos2 = single(sigs, sden_s, percentage_err, var_h, var_rho, var_h2, var_rho2,'s',
                              cden_s)
        cden_s = np.max(trhos2)
        var_h2, var_rho2 = 1.0 * ths2, 1.0 * trhos2

        #		For HI
        thg2, trhog2 = single(sigg, sden_g, percentage_err, var_h, var_rho, var_h2, var_rho2,'ag',
                              cden_ag)
        cden_ag = np.max(trhog2)
        var_h, var_rho = 1.0 * thg2, 1.0 * trhog2

        #		For H_2 thin
        thmg2, trhomg2 = single(mgsig, sden_mg, percentage_err, var_h, var_rho, var_h2, var_rho2,
                                'mg', cden_mg)
        cden_mg = np.max(trhomg2)
        var_h2, var_rho2 = 1.0 * thmg2, 1.0 * trhomg2

	#         Third Round
	#		Solving for star
        ths3, trhos3 = single(sigs, sden_s, percentage_err, var_h, var_rho, var_h2, var_rho2,'s',
                              cden_s)
        cden_s = np.max(trhos3)
        var_h, var_rho = 1.0 * ths3, 1.0 * trhos3

        #		Solving for HI
        thg3, trhog3 = single(sigg, sden_g, percentage_err, var_h, var_rho, var_h2, var_rho2, 'ag',
                              cden_ag)
        cden_ag = np.max(trhog3)
        var_h2, var_rho2 = 1.0 * thg3, 1.0 * trhog3

        #		Solving for H_2 thin
        thmg3, trhomg3 = single(mgsig, sden_mg, percentage_err, var_h, var_rho, var_h2, var_rho2,
                                'mg', cden_mg)
        cden_mg = np.max(trhomg3)
        var_h, var_rho = 1.0 * thmg3, 1.0 * trhomg3



        #		Fourth round

        #		For star
        ths4, trhos4 = single(sigs, sden_s, percentage_err, var_h, var_rho, var_h2, var_rho2,'s',
                              cden_s)
        cden_s = np.max(trhos4)
        var_h2, var_rho2 = 1.0 * ths4, 1.0 * trhos4

        #		For HI
        thg4, trhog4 = single(sigg, sden_g, percentage_err, var_h, var_rho, var_h2, var_rho2,'ag',
                              cden_ag)
        cden_ag = np.max(trhog4)
        var_h, var_rho = 1.0 * thg4, 1.0 * trhog4

        #		For H_2 thin
        thmg4, trhomg4 = single(mgsig, sden_mg, percentage_err, var_h, var_rho, var_h2, var_rho2,
                                'mg', cden_mg)
        cden_mg = np.max(trhomg4)
        var_h2, var_rho2 = 1.0 * thmg4, 1.0 * trhomg4




	print j
        escape, min_err = error(ths2, trhos2, ths4, trhos4, thg2, trhog2, thg4, trhog4, thmg2, trhomg2, thmg4, trhomg4, percentage_err)

        #		pylab.figure();	pylab.plot(ths1, trhos1, 'r-');	pylab.plot(ths2, trhos2, 'b-');	pylab.title('star'); pylab.show(); pylab.figure(); pylab.plot(thg1, trhog1, 'r-'); pylab.plot(thg2, trhog2, 'b-'); pylab.title('Atomic gas'); pylab.show(); pylab.figure(); pylab.plot(thmg1, trhomg1, 'r-'); pylab.plot(thmg2, trhomg2, 'b-');	pylab.title('Molecular gas');pylab.figure(); pylab.plot(thmgt1, trhomgt1, 'r-'); pylab.plot(thmgt2, trhomgt2, 'b-');	pylab.title('Molecular gas thick'); pylab.show(); raw_input('Halt\n')

        if escape:
        	print 'Converged coupling after %d iteration'%(j+1)
            	break
    #		else:
    #			sys.stdout.write('\rIteration running %d %4.2f'%(j,min_err))
    #			sys.stdout.flush()
    print 'out of loop'
    return ths2, trhos2, thg2, trhog2, thmg2, trhomg2
