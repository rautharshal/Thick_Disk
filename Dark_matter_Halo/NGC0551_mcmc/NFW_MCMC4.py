import dyn_py
import os
import time
import emcee
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as ss
import time
import gc
import contextlib
import corner
import multiprocessing
from multiprocessing import Pool
from dyn_py.mge_vcirc import mge_vcirc
from astropy.table import Table
from readcol import readcol
from astropy.io import ascii
from astropy import units as u
import sys
from playsound import playsound

# **********************************  INTRO  *****************************************
# To work with this code one need to have mge parameters, gas velocity(HI or CO+HI)
# and total velocity already being calculated. Stellar and dark matter contributions
# are calculated during the fitting process here. This is an example, so all the input
# files are provided here.

# INPUT FILES:
# 1) mge parameters file: surf is the central surface brightness, sigma is the dispersion
# along the major x0 axis and qobs is flattening.
# 2) total velocity profile with errors. HI kinematics from 3D-BBarolo was used.
# 3) gas velocity profile with error. Combination of CO and HI was used here.

# OUTPUT FILES:
# 1) gal+_burning_chain.pdf - plot with the burning chain
# 2) gal+_final_chain.pdf   - plot with the posterior chain
# 3) gal+_corner.pdf        - corner plot with final parameter distributions
# 4) gal+_Table_Param.txt   - table with the parameter values and their percentiles
# 5) gal+_MCMC.pdf          - final plot with all the velocity components and fit

# All appearing plots during the run are visible only for few seconds, do not close it.
# ************************************************************************************


start_time = time.time()

# some of the settings for the future plots
plt.rcParams['xtick.major.width']=1.7
plt.rcParams['xtick.labelsize']=10
plt.rcParams['ytick.major.width']=1.7
plt.rcParams['ytick.labelsize']=10
plt.rcParams['axes.labelsize']=20
plt.rcParams["axes.edgecolor"] = "black"
plt.rcParams["axes.linewidth"]  = 1.5


################################# DATA AND PARAMETERS ##############################

gals,dists,idegs,xlim_all, ylim_all, logMstar_all, IFUp_all, QSs, NAs=\
readcol('master_table/GMRT_7gals.txt', twod = False, namecomment='#', fsep = ',')

# galaxy input parameters
gal=gals[0]
ideg=idegs[0]     # Inclination in degrees
dist=dists[0]     # distance in Mpc
xlim_in=xlim_all[0] # arcsec
ylim_in=ylim_all[0] # km/s
IFup=IFUp_all[0]  # using of IFU data in Vc

Mbh  = 1e6       # BH mass in solar masses



# read MGE parameters that are needed for modelling stellar velocity profiles
filename = 'data_input/mges/'+ gal +'_convs.txt'
mgedata  = readcol(filename, twod = False, skipline=1, namecomment='#')
surf     = mgedata[0]
sigma    = mgedata[1]
qobs     = mgedata[2]

#sys.exit()

# read Vdyn_obs, file that has observational velocity. Obtained with combining the HI and CO profiles from Barolo
# Rarcsec     - radius in arcseconds
# Vdyn_obs    - velocity in km/s
# errVdyn_obs - lower error for the velocity in km/s
# Vup         - upper error for the velocity in km/s
Rarcsec, Vdyn_obs, errVdyn_obs, Vup = readcol('data_input/Vdyns/'+ gal + '_Vdyns.txt',twod = False, skipline=1, namecomment='#')


# read Vgas, file with the gas profile
# Rgas_arcsec - radius in arcseconds
# Vgas_org    - velocity in km/s
#Rgas_arcsec, Vgas_org = readcol('data_input/vgas/'+ gal +'_gas.txt',twod = False, skipline=1, namecomment='#')
#err_gas = 11.422775775685405 # uniform error for the gas
Rgas_arcsec, Vgas_org, err_gas = readcol('data_input/vgas/'+ gal +'_vgas_werr.txt',twod = False, skipline=1, namecomment='#')
Vgas = np.interp(Rarcsec, Rgas_arcsec, Vgas_org) # gas velocity interpolation


N_p = 3   # number of parameters for this model, needed for reduced chi^2 calculation


################################################## MCMC ###########################################################

# fitting function
def log_halo(p, Rarcsec):

# p defines parameters used in fitting
# Rarcsec will define the final number of radiuses.
# In this example radiuses(Rarcsec) are taken from total velocity.

    ml = p[0]       # Mass to light ratio
    M200 = p[1]     # M_200 from NFW model
    c = p[2]        # c from NFW model

    # Vc stellar
    Vstar=mge_vcirc(surf*ml, sigma, qobs, ideg, Mbh, dist, Rarcsec)

    # V_dm halo
    Vhalo=func_NFW_mc(Rarcsec, M200, c, dist)

    # Vdyn model
    Vdyn_mod=np.sqrt(Vgas**2+Vstar**2+Vhalo**2)

    # Priors, that defines the limit for the parameter. Also the last two has a lognormal distributions around approximate numbers.
    priors = ss.uniform.logpdf(ml, loc=0.01, scale=100)+\
             ss.uniform.logpdf(M200, loc=0.01, scale=1e20)+\
             ss.uniform.logpdf(c, loc=0.01, scale=1000) +\
             ss.norm.logpdf(M200, loc=(10**11.225), scale=(10**11.225 * 10**0.7 -10**11.225))+\
             ss.norm.logpdf(c, loc=10, scale=5)+\
	     ss.norm.logpdf(ml, loc=0.6, scale=0.3)




    if np.isfinite(priors) == False:
        return -np.inf

    s1 = 0
    err=(errVdyn_obs+Vup)/2.    # mean value among upper and lower error
    for i in range(len(Vdyn_mod)):
          s1+=( (Vdyn_mod[i]-Vdyn_obs[i])**2/err[i]**2)
    p1 = s1

    lp = -0.5*p1  + priors
    if np.isnan(lp):
        return -np.inf

    return lp



print 'walkers', "\n"
# Set the walkers, dimension and number
ndim, nwalkers = 3, 130
p0 = np.zeros((nwalkers,ndim))
p0[:,0] = np.random.uniform(0.01, 3, nwalkers) # mass-to-light ratio
jk = np.random.uniform(6, 15, nwalkers)  # M200
p0[:,1]= [10**jk[i] for i in range(nwalkers)]
p0[:,2] = np.random.uniform(1, 45, nwalkers) # c





fun = lambda q: np.log(1.0+q) - q/(1.0+q)
def func_NFW_mc(Rarcsec, M200, c, dist):
# Function that makes dark matter part in case of the NFW model.
# Here the function is derived to depend only from M200, c and distance.

     R = Rarcsec/3600. * np.pi/180. * dist*1e6     # arcsec to pc
     m = np.power(M200, (1./3))                    # in Msun
     A = 0.014421*m
     R200 = 20.66*m
     X = np.divide(R, R200)
     B = np.sqrt( fun(X*c)/(fun(c)*X) )
     Vhalo = A*B
     return Vhalo



print   "Run MCMC analysis"

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 					                     BURN IN CHAIN
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# setting the number of the iteration for the burning stage
burn_steps = 3000
print 'burn_step=', burn_steps
print("Burn-in...","\n")


# This is done with parallelization processes in different cores
# Number of the cores
proc = 20
with contextlib.closing(Pool(processes=proc)) as pool:
    start = time.time()                                # start the timing
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_halo,
        args=[Rarcsec], pool=pool)
    pos, prob, state = sampler.run_mcmc(p0, burn_steps, progress=True)
    end = time.time()                                  # end the timing
    multi_time = (end - start)/60.0
    print("Multiprocessing took {0:.2f} minutes".format(multi_time))


# Plotting the chain figure
print 'Chain figure'
fig = plt.figure(figsize=(11,4))
plt.suptitle(gal,fontsize=20)
plt.subplot(2,2,1)
plt.plot(sampler.chain[:,:,0].T)
plt.ylabel(r'Chain for $ML$')
plt.subplot(2,2,2)
plt.plot(sampler.chain[:,:,1].T)
plt.ylabel(r'Chain for $M_{200}$')
#plt.ylim(1e12, 1e14)
plt.yscale("log")
plt.subplot(2,2,3)
plt.plot(sampler.chain[:,:,2].T)
plt.ylabel(r'Chain for $c$')
#---------------------------------------------------
plt.subplots_adjust(top=0.9)
fig.savefig('data_output/figures/'+ gal+"_burnin_chain"+".pdf", orientation = 'landscape', format = 'pdf',bbox_inches='tight', dpi=100)
#plt.pause(2)
#plt.close()



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 					                      POSTERION CHAIN
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# setting the number of the iteration for the posterior chain
steps= 10000
print 'steps=',steps,"\n\n"
print("Running MCMC...", "\n")


sampler.reset()
with contextlib.closing(Pool(processes=proc)) as pool:
    start = time.time()
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_halo,
        args=[Rarcsec], pool=pool)
    pos,prob,state = sampler.run_mcmc(pos, steps,progress=True)
    end = time.time()
    multi_time = (end - start)/60.0
    print("Multiprocessing_2 took {0:.2f} minutes".format(multi_time))


# Plotting the posterior chain
fig = plt.figure(figsize=(12,6))
plt.suptitle(gal,fontsize=20)
plt.subplot(2,2,1)
plt.plot(sampler.chain[:,:,0].T)
plt.ylabel(r'Chain for $ML$')
plt.subplot(2,2,2)
plt.plot(sampler.chain[:,:,1].T)
plt.ylabel(r'Chain for $M_{200}$')
plt.subplot(2,2,3)
plt.plot(sampler.chain[:,:,2].T)
plt.ylabel(r'Chain for $c$')
plt.subplots_adjust(top=0.9)
fig.savefig('data_output/figures/'+gal+"_final_chain"+".pdf", orientation = 'landscape', format = 'pdf',bbox_inches='tight', dpi=100 )
#plt.pause(2)
#plt.close()


# percentage range for corner plot
perc_u = 84
perc_d = 16
quant_u = perc_u/100.0
quant_d = perc_d/100.0

flatchain_main = np.array(sampler.flatchain)
flatchain_log = np.log10(flatchain_main)


# Plotting flatchain (triangle figure) CORNER
figure = corner.corner(flatchain_log, labels=[r"$log(M/L)$", r"$log(M_{200})$",r"$log(c)$",],color = "darkblue",hist_bin_factor=1, plot_contours=True,  quantiles=[0.25, 0.5, 0.75], label_kwargs=dict(fontsize=20, fontweight='bold'),title_fmt=".1e", use_math_text=True,verbose=True, title_args={"fontsize": 10} )
figure.suptitle(gal,fontsize=20)
plt.subplots_adjust(top=0.89)
plt.subplots_adjust(hspace=0.15)
figure.tight_layout(pad=0.5)
figure.savefig('data_output/figures/'+gal+"_corner"+".pdf", orientation = 'landscape', format = 'pdf',bbox_inches='tight', dpi=100 )
#plt.pause(2)
#plt.close()




# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 					                      RESULTS
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Extract parameters, median and percentile values
ml_dist = sampler.flatchain[:,0]
M200_dist = sampler.flatchain[:,1]
c_dist = sampler.flatchain[:,2]

ml_med = np.median(sampler.flatchain[:,0])
M200_med = np.median(sampler.flatchain[:,1])
c_med = np.median(sampler.flatchain[:,2])

ml_plus = np.percentile(sampler.flatchain[:,0], perc_u)- np.median(sampler.flatchain[:,0])
M200_plus = np.percentile(sampler.flatchain[:,1], perc_u)- np.median(sampler.flatchain[:,1])
c_plus = np.percentile(sampler.flatchain[:,2], perc_u)- np.median(sampler.flatchain[:,2])

ml_minus = np.median(sampler.flatchain[:,0]) - np.percentile(sampler.flatchain[:,0],  perc_d)
M200_minus = np.median(sampler.flatchain[:,1]) - np.percentile(sampler.flatchain[:,1],  perc_d)
c_minus = np.median(sampler.flatchain[:,2]) - np.percentile(sampler.flatchain[:,2],  perc_d)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 					                   Tables with all parameters
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename='data_output/'+ gal + '_Table_Param'+ '.txt'
parameters = ['ML', 'M200', 'c']

medians = [ ml_med, M200_med, c_med]

# Array for the upper percentiles
ups = [ ml_plus, M200_plus, c_plus]

# Array for the lower percentiles
lws = [ ml_minus, M200_minus, c_minus]

t = Table([parameters, medians, ups, lws], names = ('parameter','median','up','lw'))
t.write(filename, format = 'ascii', overwrite="True")
#-------------------------------------------------

# Printing the resulting values for the parameters
print '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&'
print 'Fitted parameters'
print 'NGC',gal
print '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', "\n"

print 'ml: ', np.median(sampler.flatchain[:,0]), \
    '+', np.percentile(sampler.flatchain[:,0],  perc_u) - np.median(sampler.flatchain[:,0]),\
    '-', np.median(sampler.flatchain[:,0]) - np.percentile(sampler.flatchain[:,0],  perc_d)
print('M200: %.5e + %.5e - %.5e' %( np.median(sampler.flatchain[:,1]), np.percentile(sampler.flatchain[:,1],  perc_u) - np.median(sampler.flatchain[:,1]), np.median(sampler.flatchain[:,1]) - np.percentile(sampler.flatchain[:,1],  perc_d)) )
print('c: %.5e + %.5e - %.5e ' %( np.median(sampler.flatchain[:,2]), np.percentile(sampler.flatchain[:,2],  perc_u) - np.median(sampler.flatchain[:,2]), np.median(sampler.flatchain[:,2]) - np.percentile(sampler.flatchain[:,2],  perc_d))  )
print("\n")



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 					                    VELOCITIES
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Calculating the velocities, parallelization is used for the stellar part.

# calculate the final Vc stellar
start = time.time()
step =  len(ml_dist)/proc
finish = 15*step
left = len(ml_dist) - end

def do(b, e, st):
    w = []
    for i in range(b, e):
         w.append(mge_vcirc(surf* ml_dist[i], sigma, qobs, ideg, Mbh, dist, Rarcsec))
    st.put(w)

st = multiprocessing.Queue()
processes = []
results = []

for i in range(0, finish, step):
    p = multiprocessing.Process(target=do, args=(i,i+step, st))
    if i % 1000 == 0:
       print("GO")
    processes.append(p)
    p.start()

for r in processes:
    ret = st.get()
    results.append(ret)

for r in processes:
    p.join()

print(np.shape(results))
end = time.time()
multi_time = (end - start)/60.0
print("Velocity fin for stars took {0:.2f} minutes".format(multi_time))

part_1 = np.reshape(results, (finish, len(Rarcsec)))

part = list(part_1)
for i in range(finish, len(ml_dist)):
    part.append(mge_vcirc(surf* ml_dist[i], sigma, qobs, ideg, Mbh, dist, Rarcsec))

Vstar_fin = np.reshape(part, (len(ml_dist), len(Rarcsec))).T




# calculate the final Vhalo
Vhalo_fin=np.zeros([len(Rarcsec),len(M200_dist)])

print 'Calculate Vhalo ...'
for j in range(len(M200_dist)):
  if j % 1000 == 0:
     print '%i / %i' % (j, len(M200_dist))
  Vhalo=func_NFW_mc(Rarcsec, M200_dist[j], c_dist[j], dist)
  Vhalo_fin[:,j]=Vhalo


# calculate the final dynamical velocity
Vdyn_fin = np.zeros([len(Rarcsec),len(M200_dist)])

print 'Calculate V, model ...'
for k in range(len(Rarcsec)):
  if k % 10 == 0:
     print '%i / %i' % (k, len(Rarcsec))
  Vdyn=np.sqrt(Vgas[k]**2+Vstar_fin[k,:]**2+Vhalo_fin[k,:]**2)
  print(len(Vdyn))
  Vdyn_fin[k,:]=Vdyn



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 					                      MASSES
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#po = (dist * 1.00)/206     # coefficient to convert arcsec to kpc
po = (np.pi/(18*36))*dist     # coefficient to convert arcsec to kpc

# function to calculate the mass
def mass(R, v):
    # R - last distance in arcsec
    # v - the last value of velocity in km/s
    # Equation M = ( R * v^2 )/ (G) G in pc M_sol^(-1) (km/s)^2. Return np.log10(M) in solar unity
    G = 4.30091e-6 #km^2-kpc/(sec^2-Msun)
    mass_tot = ((R*po * (v**2))/G)
    return np.log10(mass_tot)
#    enum = R * po * (v**2) * 1e6
#    den = 2 * 4.3
#    return np.log10(enum/den)


# median profiles of the velocities
V_dyn_mod = np.nanmedian(Vdyn_fin, axis = 1)
Vhalo = np.nanmedian(Vhalo_fin, axis = 1)
Vstars = np.nanmedian(Vstar_fin, axis = 1)

# masses at the last observational point
mg = mass(Rarcsec[len(Rarcsec)-1],   Vgas[len(Vgas)-1])                    # mass of the gas
ms = mass(Rarcsec[len(Rarcsec)-1],  Vstars[len(Vstars)-1] )                # stellar mass
mdyn_mod = mass(Rarcsec[len(Rarcsec)-1],  V_dyn_mod[len(V_dyn_mod)-1])     # total mass modelled
mdyn_obs = mass(Rarcsec[len(Rarcsec)-1],  Vdyn_obs[len(Vdyn_obs)-1])       # total mass observational

# dark matter mass, modelled and observational
mh_mod = np.log10((10.**mdyn_mod) - (10.**mg) - (10.**ms))
mh_obs = np.log10((10.**mdyn_obs) - (10.**mg) - (10.**ms))


# calculating the reduced chi^2
Points = len(Vdyn_obs)
err = (errVdyn_obs+Vup)/2.
s3 = 0
for i in range(len(V_dyn_mod)):
      s3+=( ((V_dyn_mod[i]-Vdyn_obs[i])**2.)/(err[i]**2.))
chi2 = s3/(Points-N_p)
print("chi^2_red ", chi2 )


# calculating dark matter fractions
frac_obs = (10.**mdyn_obs-10.**mg-10.**ms)/10.**mdyn_obs *100
frac_mod = (10.**mdyn_mod-10.**mg-10.**ms)/10.**mdyn_mod *100

end_time = time.time()
time_taken = (end_time - start_time)/60.0
print 'Time taken:%.2e'  %time_taken, 'minutes'

outfile = 'data_output/' + gal + '_mass_param.txt'
f = open(outfile, "w")
saveout = sys.stdout
sys.stdout = f
print 'gas mass = ', mg
print 'star mass = ', ms
print 'dynamical mass from model = ', mdyn_mod
print 'dynamical mass from observation = ', mdyn_obs
print 'helo mass from model = ', mh_mod
print 'helo mass from observation = ', mh_obs
print '\n'
print 'dark matter fraction from observation = ', frac_obs
print 'dark matter fraction from model = ', frac_mod

print '\n'
print '\n'
print "chi^2_red ", chi2
f.close()

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 					                     PLOT THE FINAL FIGURE
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Calculate the ERRORS
eVctop3 = np.nanpercentile(Vstar_fin, 75, axis = 1) - np.nanmedian(Vstar_fin, axis = 1)
eVcbot3 = np.nanmedian(Vstar_fin, axis = 1) - np.nanpercentile(Vstar_fin, 25, axis = 1)

eVctop2 = np.nanpercentile(Vhalo_fin, 75, axis = 1) - np.nanmedian(Vhalo_fin, axis = 1)
eVcbot2 = np.nanmedian(Vhalo_fin, axis = 1) - np.nanpercentile(Vhalo_fin, 25, axis = 1)

eVctop_p = np.nanpercentile(Vdyn_fin, perc_u, axis = 1) - np.nanmedian(Vdyn_fin, axis = 1)
eVcbot_p = np.nanmedian(Vdyn_fin, axis = 1) - np.nanpercentile(Vdyn_fin, perc_d, axis = 1)


# additional settings for figure
plt.rcParams["axes.edgecolor"] = "black"
plt.rcParams["axes.linewidth"]  = 1.5



# FINAL FIGURE
fig = plt.figure(figsize=(12,6))
ax = fig.add_subplot(111)

# CO from 3D-BBarolo
#ax.errorbar(Rarcsec[:IFup], Vdyn_obs[:IFup], marker='8',markersize=9, markerfacecolor='gray', yerr = (errVdyn_obs[:IFup], Vup[:IFup]), label='Vc, CO', color = 'black', zorder= 10, lw=2)

# HI from 3D-BBarolo
ax.errorbar(Rarcsec[IFup:], Vdyn_obs[IFup:], marker='s',markersize=9, markerfacecolor='silver', yerr = (errVdyn_obs[IFup:], Vup[IFup:]), label='V, HI', color = 'black', zorder= 10, lw=2)
plt.plot(Rarcsec, Vdyn_obs, lw=2, color='black')


# Vgas curve
#ax.scatter(Rarcsec, Vgas, marker='v', label='Vgas',s=85, edgecolors='black', linewidths=1.4, color = 'green', zorder=3 )
#ax.plot(Rarcsec, Vgas, color='black', lw=2, zorder=2 )
#plt.fill_between(Rarcsec, Vgas+err_gas, Vgas-err_gas, alpha=0.3,  facecolor='green', linewidth=2, linestyle='-', antialiased=True)
ax.errorbar(Rgas_arcsec, Vgas_org, yerr=err_gas, marker='v',label='Vgas', markersize=9,markerfacecolor='green', color = 'black', lw=2)

#  modelled Vhalo curve
ax.plot(Rarcsec,np.nanmedian(Vhalo_fin, axis = 1), '-.', color='blue', label='Vhalo',lw=2.5, zorder=0)
plt.fill_between(Rarcsec, np.nanmedian(Vhalo_fin, axis = 1)+eVctop2, np.nanmedian(Vhalo_fin, axis = 1)-eVcbot2, alpha=0.3,  facecolor='blue', linewidth=2, linestyle='-', antialiased=True)

# modelled Vstellar curve
ax.plot(Rarcsec,np.nanmedian(Vstar_fin, axis = 1), '--',  color='orange', label='Vstars',lw=2.5,zorder=4)
plt.fill_between(Rarcsec, np.nanmedian(Vstar_fin, axis = 1)+eVctop3, np.nanmedian(Vstar_fin, axis = 1)-eVcbot3, alpha=0.3,  facecolor='orange', linewidth=2, linestyle='-', antialiased=True)

#  modelled Vdyn curve
ax.plot(Rarcsec,np.nanmedian(Vdyn_fin, axis = 1), color='red', label = 'V, model', lw=2.5,zorder=7)
plt.fill_between(Rarcsec, np.nanmedian(Vdyn_fin, axis = 1)+eVctop_p, np.nanmedian(Vdyn_fin, axis = 1)-eVcbot_p, alpha=0.3,  facecolor='red', linewidth=2, linestyle='-', antialiased=True)


ax.set_xlabel(r"R [arcsec]")
ax.set_ylabel(r"V [km/s]")

ax.tick_params(axis = 'x', direction='in', length=5, width=2, colors='black', labelsize=20)
ax.tick_params(axis = 'y', direction='in', length=5, width=2, colors='black', labelsize=20)


leg = plt.legend(loc = "best", fontsize=15, ncol=2, bbox_to_anchor=(0.6, 0.25), fancybox=True, framealpha=1, shadow=True, borderpad=0.5)
leg.get_frame().set_linewidth(1.5)
leg.get_frame().set_edgecolor('black')



ax.set_xlim(0,xlim_in)
ax.set_ylim(0, ylim_in)
#ax.xaxis.set_ticks([ 10, 20, 30,  40, 50, 60, 70, 80, 90 ])
#ax.yaxis.set_ticks([0, 40,  80, 120,  160,  200,  240,  280, 320])

'''
textstr = '\n'.join((
    r'$\rm M_{halo}^{obs}(Rmax) = 10^{%.2f}\rm \, M_{\odot}$' % (mh_obs,  ),
    r'$\rm M_{halo}^{mod}(Rmax) = 10^{%.2f}\rm \, M_{\odot}$' % (mh_mod,  ),
    r'$\rm M_{stars}(Rmax) = 10^{%.2f}\rm \, M_{\odot}$' % (ms,  ),
    r'$\rm M_{gas}(Rmax) = 10^{%.2f}\rm \, M_{\odot}$' % (mg, ),
    r'$\rm M_{dyn}^{obs}(Rmax) = 10^{ %.2f}\rm \, M_{\odot}$' % (mdyn_obs, ),
    r'$\rm M_{dyn}^{mod}(Rmax) = 10^{ %.2f}\rm \, M_{\odot}$' % (mdyn_mod, )))

# Box with the masses
props = dict(boxstyle='square', facecolor='gold', pad=0.4, alpha=0.6)
ax.text(15, 165, textstr,  fontsize=10,  horizontalalignment='center', verticalalignment='center',  bbox=props)
'''
# Additional box with overall information
textstr2 = ';   '.join((r'Galaxy '+gal,  r'MGE: IRAC 3.6 micron'))
props = dict(boxstyle='square', facecolor='deeppink', edgecolor='black', alpha=0.5)
ax.text(30, 295, textstr2,  fontsize=20,  horizontalalignment='center', verticalalignment='center',  bbox=props)

'''
# Box with frations and chi^2
textstr3 = '\n'.join(( r'$\chi^2_{\nu}$ = %.2f ' % (chi2, ),
r'$ \rm f_{DM}^{obs} = %.2f \,  $'  % (frac_obs, ) +  r'$\%$',
r'$ \rm f_{DM}^{mod} = %.2f \, $' % (frac_mod, ) +  r'$\%$'))
props3 = dict(boxstyle='round4', facecolor='white', edgecolor='darkred', pad=0.5, alpha=0.9)
ax.text(75, 287, textstr3,  fontsize=10.,  horizontalalignment='center',  verticalalignment='center',  bbox=props3)
'''

'''
# Box with the parameters
textstr4 = '\n'.join((r'Parameters',
r'$ \rm ML = %.2f \rm \, \pm %.2f \, \, \, [M_{\odot}/L_{\odot}]$' % (ml_med, (ml_plus+ml_minus)/2., ),
r'$\rm c = %.2f \rm \, \pm %.2f $'  % (c_med, (c_plus+c_minus)/2.,),
r'$\rm M_{200} = 10^{%.2f} \rm \, \pm %.2e \, \, \, [M_{\odot}] $' % (np.log10(M200_med), (M200_plus+M200_minus)/2.,)  ))
props5 = dict(boxstyle='sawtooth', facecolor='lightskyblue',  edgecolor='darkcyan', pad=0.7, alpha=0.4)
ax.text(75, 83, textstr4,  fontsize=9,  horizontalalignment='center',color='navy', fontweight='bold', verticalalignment='center',  bbox=props5)
'''

ax2 = ax.twiny()
ax1Xs = ax.get_xticks()

ax2Xs = []
for X in ax1Xs:
    ax2Xs.append(X * po)

ax2Xs = ["%.2f"%item for item in ax2Xs]
ax2.set_xticks(ax1Xs)
ax2.set_xbound(ax.get_xbound())
ax2.set_xticklabels(ax2Xs)
ax2.tick_params( direction='in', length=5, width=2, colors='black', labelsize=20)
ax2.set_xlabel(r"R [Kpc]", color="black")
fig.savefig('data_output/figures/'+gal+'_MCMC'+'.pdf', orientation = 'landscape', format = 'pdf',bbox_inches='tight', dpi=150 )
#plt.pause(3)
#plt.close()



#############################################################################
# upper and lower limit of the velocities decided by perc_u and perc_d respectively
##############################################################################
Vstars_up = np.nanpercentile(Vstar_fin, 75, axis = 1)
Vstars_dwn = np.nanpercentile(Vstar_fin, 25, axis = 1)

V_dyn_mod_up = np.nanpercentile(Vdyn_fin, perc_u, axis = 1)
V_dyn_mod_dwn = np.nanpercentile(Vdyn_fin, perc_d, axis = 1)

Vhalo_up = np.nanpercentile(Vhalo_fin, 75, axis = 1)
Vhalo_dwn = np.nanpercentile(Vhalo_fin, 25, axis = 1)

ms_up = mass(Rarcsec[len(Rarcsec)-1],  Vstars_up[len(Vstars_up)-1] )                # stellar mass up
ms_dwn = mass(Rarcsec[len(Rarcsec)-1],  Vstars_dwn[len(Vstars_dwn)-1] )                # stellar mass down
mdyn_mod_up = mass(Rarcsec[len(Rarcsec)-1],  V_dyn_mod_up[len(V_dyn_mod_up)-1])     # upper mass modelled
mdyn_mod_dwn = mass(Rarcsec[len(Rarcsec)-1],  V_dyn_mod_dwn[len(V_dyn_mod_dwn)-1])     # lower mass modelled
mh_mod_up = np.log10((10.**mdyn_mod_up) - (10.**mg) - (10.**ms_up))
mh_mod_dwn = np.log10((10.**mdyn_mod_dwn) - (10.**mg) - (10.**ms_dwn))

mg_up = mass(Rarcsec[len(Rarcsec)-1],   (Vgas[len(Vgas)-1]+err_gas[len(Vgas)-1]))                    # mass of the gas up
mg_dwn = mass(Rarcsec[len(Rarcsec)-1],   (Vgas[len(Vgas)-1]-err_gas[len(Vgas)-1]))                    # mass of the gas down

outfile = 'data_output/' + gal + '_mass_param.txt'
f = open(outfile, "w")
saveout = sys.stdout
sys.stdout = f
print 'gas mass: median, up, down ' , mg, mg_up, mg_dwn
print 'star mass: median, up, down ', ms, ms_up, ms_dwn
print 'dynamical mass from model: mean, up, down = ', mdyn_mod, mdyn_mod_up,  mdyn_mod_dwn
print 'dynamical mass from observation = ', mdyn_obs
print 'helo mass from model: median, up, down = ', mh_mod, mh_mod_up, mh_mod_dwn
print 'helo mass from observation = ', mh_obs
print '\n'
print 'dark matter fraction from observation = ', frac_obs
print 'dark matter fraction from model = ', frac_mod

print '\n'
print '\n'
print "chi^2_red ", chi2
f.close()

outfile = 'data_output/' + gal + '_mass_param_nonlog.txt'
f = open(outfile, "w")
saveout = sys.stdout
sys.stdout = f
print 'gas mass: %.2e + %.2e - %.2e' %(10**mg, 10**mg_up-10**mg, 10**mg-10**mg_dwn)
print 'star mass:  %.2e + %.2e - %.2e' %(10**ms, 10**ms_up-10**ms, 10**ms-10**ms_dwn)
print 'dynamical mass from model:  %.2e + %.2e - %.2e' %(10**mdyn_mod, 10**mdyn_mod_up-10**mdyn_mod, 10**mdyn_mod-10**mdyn_mod_dwn)
print 'dynamical mass from observation: %.2e ' %10**mdyn_obs
print 'helo mass from model: %.2e + %.2e - %.2e' %(10**mh_mod, 10**mh_mod_up-10**mh_mod, 10**mh_mod-10**mh_mod_dwn)
print 'helo mass from observation: %.2e' %10**mh_obs
print '\n'
print 'dark matter fraction from observation = ', frac_obs
print 'dark matter fraction from model = ', frac_mod

print '\n'
print '\n'
print "chi^2_red ", chi2
f.close()

outfile = 'data_output/' + gal + '_velocity.txt'
f = open(outfile, "w")
saveout = sys.stdout
sys.stdout = f

print ('#Radius(arcsec) Vstar(km/s) eupVstar(km/s) edwnVstar(km/s) Vhalo(km/s) eupVhalo(km/s) edwnVhalo(km/s) Vtot_mod(km/s) eupVtot_mod(km/s) edwn_Vtot_mod(km/s)')
for i in range(len(Rarcsec)):
      print Rarcsec[i], Vstars[i], eVctop3[i], eVcbot3[i], Vhalo[i], eVctop2[i], eVcbot2[i], V_dyn_mod[i], eVctop_p[i], eVcbot_p[i]

f.close()

outfile = 'data_output/' + gal + '_velocity.txt'
f = open(outfile, "w")
saveout = sys.stdout
sys.stdout = f

print ('#Radius(arcsec) Vstar(km/s) eupVstar(km/s) edwnVstar(km/s) Vhalo(km/s) eupVhalo(km/s) edwnVhalo(km/s) Vtot_mod(km/s) eupVtot_mod(km/s) edwn_Vtot_mod(km/s)')
for i in range(len(Rarcsec)):
	print Rarcsec[i], Vstars[i], eVctop3[i], eVcbot3[i], Vhalo[i], eVctop2[i], eVcbot2[i], V_dyn_mod[i], eVctop_p[i], eVcbot_p[i]

f.close()

