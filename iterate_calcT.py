#!/usr/bin/env python2
import sys
import os
import math
import linecache
import struct
import time

import numpy as np

import matplotlib as mpl
mpl.use('Agg')   # to avoid the window in plotting

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogLocator

#sys.path.append('/turquoise/users/sjin/python/')


######### ######### ######### ######### ######### ######### 
###### the function that read one temperature file (binary)

def readtemperature_1d(filename):
    debug = False
 
    grid = linecache.getline('amr_grid.inp', 6)
    ngrids = [int(s) for s in grid.split() if s.isdigit()]
    nrad  = ngrids[0]
    nthet = ngrids[1]
    nphi  = ngrids[2]
    
    #temp_3d = np.array(temp_3d).reshape(nphi*nthet*nrad)

    head = np.fromfile(filename, count=4, dtype=np.int64)
    if debug:
	print "head[0]", head[0]
        print "head[1]", head[1]
        print "head[2]", head[2]
        print "head[3]", head[3]

    if head[2]!=(nrad*nthet*nphi):
        print ' ERROR'
        print ' Number of grid points in '+filename+' is different from that in amr_grid.inp'
        quit()

    if head[1] == 8:
        f = open(filename, "rb")  # reopen the file
        f.seek(32, os.SEEK_SET)  # seek
        temp_1d = np.fromfile(f, count=-1, dtype=np.float64)
    elif head[1]==4:
        f = open(filename, "rb")  # reopen the file
        f.seek(32, os.SEEK_SET)  # seek
        temp_1d = np.fromfile(f, count=-1, dtype=np.float32)
    else:
        print 'ERROR'
        print 'Unknown datatype in '+filename
    f.close()

    print " dust_temp reading: done."

    return temp_1d

########################################################
############ the function do one run

def onerun(result_dir):

    gas_radmc_path = "./gas_line_radmc"
    
    debug = False        # info printing
    
    # creat amr_grid.inp and initial dust_temperature.bdat
    os.system('./manu_gas_model_create.py')
    
    # creat gas_...surf_density.ascii using 5 parameters
    os.system('./gas.manu_density_based_amr_hpc.py')
    
    # calculate the dust_density.binp 
    # NOTE make sure gas freezing out and photo dissosciation are off.
    os.system('./gas_numberdensity_create_sph.py')
    
    # run radmc3d for the first time
    os.system('cp dustkappa_dust.inp.temp dustkappa_dust.inp -f')
    os.system('./radmc3d mctherm')
    
    
    #########################################
    # read amr_grid to get the grid info
    
    grid = linecache.getline('amr_grid.inp', 6)
    ngrids = [int(s) for s in grid.split() if s.isdigit()]
    nrad  = ngrids[0]
    nthet = ngrids[1]
    nphi  = ngrids[2]
    
    
    for i in range(10):
        os.system('./gas_numberdensity_create_sph.py')
        os.system('mv dust_temperature.bdat dust_temperature.old.bdat')
        os.system('./radmc3d mctherm')
        filename_old = "%s%s" % ("./","dust_temperature.old.bdat")
        filename_new = "%s%s" % ("./","dust_temperature.bdat")
        temp_1d_old = readtemperature_1d(filename_old)
        temp_1d_new = readtemperature_1d(filename_new)
        temp_diff = (temp_1d_new - temp_1d_old)/temp_1d_new
        text = "%s%s%s%s" % (np.average(temp_diff), "  ", temp_diff.max(), "\n")
        with open("nI_difMax_diffAve.dat", "a") as myfile:
            myfile.write(text)
        text = ""
        print i, "  ", np.average(temp_diff), "  ", temp_diff.max()
        time.sleep(10)
        if (np.average(temp_diff) < 0.03) and (temp_diff.max() < 0.03) and  (i > 1):
            break
    
    os.system('cp amr_grid.inp ./gas_line_radmc')
    os.system('./dust.manu_density_based_amr.py')
    os.system('cp dustkappa_dust.inp.imag dustkappa_dust.inp -f')
    
    os.system('./radmc3d image lambda 866.9627371918420522  npix 400 incl 51 posang 151 sizeau 1400 noline')
    os.system('./image_conv_azimth_extract.py')
    cmd = "mv AU_fluxJy ./results/%s/Au_fluxJy_dust"%(result_dir)
    os.system(cmd)
    os.system('mv image.out image.out.12co')
    
    #os.system('./radmc3d image lambda 906.8455757791763290 npix 400 incl 51 posang 151 sizeau 1400 noline')
    #os.system('mv image.out image.out.13co')
    
    #os.system('./radmc3d image lambda 910.3083452186419890 npix 400 incl 51 posang 151 sizeau 1400 noline')
    #os.system('mv image.out image.out.c18o')
    
    os.system('cp dust_density.binp ./gas_line_radmc')
    os.system('cp dust_temperature.bdat ./gas_line_radmc')
    
    parentdir = "../"
    
    os.system('./gas_numberdensity_create_sph_12CO.py')
    os.chdir( gas_radmc_path )
    os.system('cp ../camera_wavelength_micron.inp_12CO ./camera_wavelength_micron.inp')
    os.system('cp ../image.out.12co ./image.out_lambda_0')
    os.system('./image_line.sh')
    os.system('./subtract_convolveGuassian_and_plot_12CO.py')
    cmd = "mv AU_fluxJy ../results/%s/Au_fluxJy_12CO"%(result_dir)
    os.system(cmd)
    os.chdir( parentdir )
    
    #os.system('./gas_numberdensity_create_sph_13CO.py')
    #os.chdir( gas_radmc_path )
    #os.system('cp ../camera_wavelength_micron.inp_13CO ./camera_wavelength_micron.inp')
    #os.system('cp ../image.out.13co ./image.out_lambda_0')
    #os.system('./image_line.sh')
    #os.system('./subtract_convolveGuassian_and_plot_13CO.py')
    #cmd = "mv AU_fluxJy ../results/%s/Au_fluxJy_13CO"%(result_dir)
    #os.system(cmd)
    #os.chdir( parentdir )
    
    #os.system('./gas_numberdensity_create_sph_C18O.py')
    #os.chdir( gas_radmc_path )
    #os.system('cp ../camera_wavelength_micron.inp_C18O ./camera_wavelength_micron.inp')
    #os.system('cp ../image.out.c18o ./image.out_lambda_0')
    #os.system('./image_line.sh')
    #os.system('./subtract_convolveGuassian_and_plot_C18O.py')
    #cmd = "mv AU_fluxJy ../results/%s/Au_fluxJy_C18O"%(result_dir)
    #os.system(cmd)
    #os.chdir( parentdir )
    
#  main program
run_grid = np.genfromtxt('run_grid.txt', dtype=float)
print run_grid.shape
for i in range(run_grid.shape[0]):
    result_dir = ('%04d'%(run_grid[i,0]))
    print result_dir
    cmd = "sed 's/mdiskSS/%s/;s/rcSS/%s/;s/gammaSS/%s/;s/rc_atanSS/%s/;s/gamma_atanSS/%s/' gas.manu_density_based_amr_hpc.py.ss > gas.manu_density_based_amr_hpc.py"%(run_grid[i,1], run_grid[i,2], run_grid[i,3], run_grid[i,4], run_grid[i,5])
    os.system(cmd)
    os.system('chmod +x gas.manu_density_based_amr_hpc.py')
    #cmd = "sed 's/mdiskSS/%s/;s/rcSS/%s/;s/gammaSS/%s/;s/rc_atanSS/%s/;s/gamma_atanSS/%s/' dust.manu_density_based_amr.py.ss > dust.manu_density_based_amr.py"%(run_grid[i,1], run_grid[i,2], run_grid[i,3], run_grid[i,4], run_grid[i,5])
    #os.system(cmd)
    os.system('chmod +x dust.manu_density_based_amr.py')
    result_dir_path = "./results/%s"%(result_dir)
    if not os.path.exists(result_dir_path):
        os.mkdir(result_dir_path, 0755)
        print result_dir_path, " is created", "\n"
    onerun(result_dir)
