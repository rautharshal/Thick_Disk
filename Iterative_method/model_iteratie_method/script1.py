import os
import sys
import numpy as np
import break_loop

if len(sys.argv) !=6:
    print '<base sig file> <txt file> <restart?(y/n)> <prefix> <distance to galaxy>'
    sys.exit()

obase_sig = sys.argv[1]
base_sig = sys.argv[1]
last_sig = sys.argv[1]
star_file = sys.argv[2] # 551_s_vdisp.txt
restart = sys.argv[3]
prefix = sys.argv[4] #ngc551
d = sys.argv[5]

for trm in range(20):
    itr_val = str((trm * 1 + 1))
    cmd1 = 'python cons_xco_thick.py NGC0551_hi_surf.dat 551_h2_surface_density_1.txt 551_s_surface_density.txt '+star_file+' 551_h2_vel_disp.txt '+d+' 2000 11000 100 12. p 30 '+prefix+'_sden'+itr_val+'.fits ' + prefix+'_agden'+itr_val+'.fits ' + prefix+'_mgden'+itr_val+'.fits ' + prefix+'_sig'+itr_val+'.fits'
    outfile = prefix+'_sden'+itr_val+'.fits'

    exe = os.system('ls ' + outfile)
    if exe != 0 or restart == 'y':
        print cmd1
        #        raw_input('Halt\n')
        os.system(cmd1)
    else:
        print
        print '%s file exists.' % outfile
        print 'Skipping 3comp diffsolve.'
        print
    # 	Building the cubes
    print 'Building the cubes'

    itr_val = str((trm * 1 + 1)) # pixresolution # 85 for v1200 150 for v500 2.16
    cmd_4comp = 'python sclco_4comp_momgen_cube_v3.py ' +prefix+'_sden'+itr_val+'.fits ' + base_sig + ' 72 50 700 2.16 2.16 0. 5192 10 63 30 ' + prefix + '_' + itr_val

    outfile = prefix + '_final_cube_' + itr_val + '.fits'
    exe = os.system('ls ' + outfile)
    if exe != 0 or restart == 'y':
        print cmd_4comp

        f = open('command_4comp', 'w')
        saveout = sys.stdout
        sys.stdout = f
        print cmd_4comp
        sys.stdout = saveout
        f.close()
        cmd_chmod = 'chmod a+x command_4comp'
        os.system(cmd_chmod)
        cmd4 = 'mpirun -np 20 bash command_4comp'
        # raw_input('Halt\n')
        os.system(cmd4)
    else:
        print
        print '%s file exists.' % outfile
        print 'Skipping momgen_cube.'
        print
#Collecting the cubelets and Convolving and reshampling.

    #raw_input('Halt\n')

    itr_val = str((trm*1+1))
    cmd_collect = 'python collect_sclco_v2.py ' + prefix+'_'+itr_val+'_cube0.fits ' + prefix+'_'+itr_val+'_kernel.fits ' + prefix+'_convl_'+itr_val

    outfile = prefix+'_final_cube_'+itr_val+'.fits'
    exe = os.system('ls '+outfile)
    if exe != 0 or restart == 'y':
        print cmd_collect
        f = open('command_collect','w')
        saveout = sys.stdout
        sys.stdout = f
        print cmd_collect
        f.close()
        sys.stdout = saveout
        cmd_chmod = 'chmod a+x command_collect'
        os.system(cmd_chmod)
        cmd5 = 'mpirun -np 10 bash command_collect'
        #raw_input('Halt\n')
        os.system(cmd5)
    else:
        print
        print '%s file exists.'%outfile
        print 'Skipping convolution_cube.'
        print
    # 	Collecting the convolved cubelets.

    itr_val = str((trm * 1 + 1))
    cmd6 = 'python collect_convl.py ' + prefix + '_convl_' + itr_val + '_ccube0.fits ' + prefix + '_final_cube_' + itr_val + '.fits'

    outfile = prefix + '_final_cube_' + itr_val + '.fits'
    exe = os.system('ls ' + outfile)
    if exe != 0 or restart == 'y':
        print cmd6
        os.system(cmd6)
        #raw_input('Halt\n')
    else:
        print
        print '%s file exists.' % outfile
        print 'Skipping collecting_cubes.'
        print
    # 	Removing the intermediate files.
    #raw_input('Halt\n')

    itr_val = str((trm * 1 + 1))
    cmd_rm = 'rm ' + prefix + '_convl_' + itr_val + '_ccube*.fits'
    os.system(cmd_rm)
    cmd_rm = 'rm ' + prefix + '_' + itr_val + '_cube*.fits'
    os.system(cmd_rm)
    cmd_rm = 'rm ' + prefix + '_' + itr_val + '_kernel.fits'
    os.system(cmd_rm)
    cmd_rm = 'rm command_4comp command_collect'
    os.system(cmd_rm)

    #	Making MOMNT maps instead of decomposing Gaussian-Hermite

    itr_val = str((trm * 1 + 1))
    cmd_momnt = 'python momgen.py ' + prefix + '_final_cube_' + itr_val + '.fits ' + prefix + '_' + itr_val

    outfile = prefix + '_' + itr_val + '_mom2.fits'
    exe = os.system('ls ' + outfile)
    if exe != 0 or restart == 'y':
        print cmd_momnt
        os.system(cmd_momnt)
        # raw_input('Halt\n')
    else:
        print
        print '%s file exists.' % outfile
        print 'Skipping momgen.'
        print
    itr_val = str((trm * 1 + 1))
    cmd_momex = 'python mom2_extract.py ' + prefix + '_' + itr_val + '_mom2.fits ' + obase_sig + ' 500 stack_' + itr_val

    outfile = 'stack_' + itr_val + '_sig.npy'
    exe = os.system('ls ' + outfile)
    if exe != 0 or restart == 'y':
        print cmd_momex
        os.system(cmd_momex)
        # raw_input('Halt\n')
    else:
        print
        print '%s file exists.' % outfile
        print 'Skipping extract mom2 profile.'
        print


    sys.exit()
# 	Checking the breaking loop
    itr_val = str((trm * 1 + 1))
    flag = break_loop.diff(obase_sig, 'stack_' + itr_val + '_sig.npy', 3000, 10000, 0.05)
    print 'obase_sig %s' % obase_sig, 'stack_ %s' % itr_val, '_sig.npy'
    print
    print 'The value of the flag is %d' % flag
    print
    print

    # 	Creating the new sigma file
    if flag == 1:
        itr_val = str((trm * 1 + 1))
        print 'Converged after %d iteration' % (trm + 1)
        cmd_plot = 'python plot_sig.py ' + 'stack_' + itr_val + '_sig.npy ' + obase_sig + ' ' + base_sig + ' ' + '3000 10000 ' + prefix + '_' + itr_val
        print cmd_plot
        os.system(cmd_plot)

        break;
    else:
        itr_val = str((trm * 1 + 1))
        cmd9 = 'python sig_modify.py ' + obase_sig + ' ' + last_sig + ' stack_' + itr_val + '_sig.npy ' + ' nstack_' + itr_val
        print cmd9
        os.system(cmd9)
        # raw_input('Halt\n')

        base_sig = 'nstack_' + itr_val + '_sig_new.npy'
        last_sig = 'nstack_' + itr_val + '_sig_new.npy'

	"""
	arcsec = np.loadtxt(star_file,dtype=float,usecols=0)
	vel = np.loadtxt(star_file,dtype=float,usecols=1)
	dist = 1e6*d*np.pi*(arcsec/3600)/180.
	vel_m = vel*1000
	arr = [dist,vel_m]
	np.save('551_star_sig',arr)
	"""
	a = np.load(base_sig)
	d1 = float(d)
	data = np.column_stack(((a[0]*180*3600)/(1e6*d1*np.pi),a[1]/1000.))
	star_file ='551_s_vdisp'+itr_val+'.txt' 
	np.savetxt(star_file,data,delimiter='\t')
	"""
	data = np.loadtxt(star_file,delimiter='\t',dtype='float',usecols=range(2))
	nan_rows=np.isnan(data.astype('float')).any(axis=1)
	data = data[~nan_rows]
	np.savetxt(star_file,data,delimiter='\t')
	"""
#    raw_input('Halt\n')
