#!/usr/bin/env python2
import sys
import os
import linecache
import numpy as np
import math
import struct

#sys.path.append('/turquoise/users/sjin/python/')

import matplotlib as mpl
mpl.use('Agg')   # to avoid the window

##################################
#### creat gas density ascii

fit_gas  = True 
#fit_gas  = False 
#write_gas3d_density  = True 
write_gas3d_density  = False

inner_hole = False
hole_size = 7.5  # cavity size in AU

write_temp = True
tan_crit = 0.1  # two-layer model temperature boundary

## Sigma_ = SigmaC * (r/RC)**-gamma * exp( -(r/RC)**(2-gamma)) * atan((r/RC)**gamma_atan)

if fit_gas:
    Mdisk = mdiskSS         # gas disk mass in solar mass
    gamma = gammaSS          #  power-law index
    r_c_tan   = rc_atanSS     #  the radial size of the gas disk in AU
    r_c   = rcSS        #  the radial size of the gas disk in AU
    gamma_atan = gamma_atanSS

    minimum  = 0.10
    #
    gas2dust = 100
    min_number = True
    minimum_v  = 0.000
    innerdisk_r = 5     # boundary of the inner disk
    innerdisk_hr_scale = 5  # hr_innerdisk/hr
    innerdisk_D = 0.00  # density of the inner disk

#filename1 = "%s%s" % ("./","dust_temperature.bdat")

# Define some constants
AU    = 1.4959787e13   # Astronomical unit [cm]
GG    = 6.67259e-8     # Gravitational constant
Msun  = 1.9891e33      # Solar mass [g]
Mstar = 2.30           # Stellar mass in Msun



###################
######## Start Function Declaration

######### ######### ######### ######### ######### ######### 
###### the function that read one temperature file (binary)

def readtemperature(filename):
    debug = False
 
    global temp_3d, nrad, nthet, nphi
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
        temp_3dread = np.fromfile(f, count=-1, dtype=np.float64)
    elif head[1]==4:
        f = open(filename, "rb")  # reopen the file
        f.seek(32, os.SEEK_SET)  # seek
        temp_3dread = np.fromfile(f, count=-1, dtype=np.float32)
    else:
        print 'ERROR'
        print 'Unknown datatype in '+filename
    f.close()

    print " dust_temp reading: done."

    temp_3dread = np.array(temp_3dread).reshape(nphi,nthet,nrad)
    temp_3d = temp_3dread
    print " temp_3d transformed to nphi,nthet,nrad"

    return temp_3d

######## End Function Declaration
###################



################### ################### 
################### the main program
################### ################### 


################################# 
############### read the grids

grid = linecache.getline('amr_grid.inp', 6)
print grid
ngrids = [int(s) for s in grid.split() if s.isdigit()]
nrad  = ngrids[0]
nthet = ngrids[1]
nphi  = ngrids[2]

nphi_grid = nphi + 1
nrad_grid = nrad + 1

print "nrad", nrad
print "nthet", nthet
print "nphi", nphi

grids = np.genfromtxt('amr_grid.inp', skip_header=6, dtype=float)
print grids.shape

nrad_b  = nrad+1
nthet_b = nthet+1
nphi_b  = nphi+1
print "ngrids[0]", ngrids[0]
print "ngrids[1]", ngrids[1]
print "ngrids[2]", ngrids[2]

nrad_b1 = grids[0:nrad:1]
rad_grid = grids[0:nrad+1:1]
nrad_b2 = grids[1:nrad+1:1]
rad = (nrad_b1+nrad_b2)/2.0

nthet_b1 = grids[nrad+1:nrad+1+nthet:1]
thet_grid = grids[nrad+1:nrad+2+nthet:1]
nthet_b2 = grids[nrad+2:nrad+2+nthet:1]
thet = (nthet_b1+nthet_b2)/2.0

nphi_b1 = grids[nrad+2+nthet:nrad+2+nthet+nphi:1]
phi_grid = grids[nrad+2+nthet:nrad+3+nthet+nphi:1]
nphi_b2 = grids[nrad+3+nthet:nrad+3+nthet+nphi:1]
phi = (nphi_b1+nphi_b2)/2.0
 

################################# 
######  Declare the variables and read density data 

surf1g = np.zeros(nrad_grid, dtype=float)


###################################################### 
######## manually set the density profile

print "write the gas surface density ascii", "\n"
if 1.0-gamma == -1:
    SigmaC = Mdisk*Msun/(2.0*math.pi)/((r_c*AU)**2.0)/(math.log(600.0/r_c)-math.log(5.0/r_c))
if 1.0-gamma != -1:
    SigmaC = Mdisk*Msun/(2.0*math.pi)/((r_c*AU)**2.0)/((600.0/r_c)**(2.0-gamma)/(2.0-gamma)-(5.0/r_c)**(2.0-gamma)/(2.0-gamma))



if fit_gas:
    dens2d_name = "%s%s" % ("gas_surface_density", ".ascii")
    print "write the gas surface density ascii", "\n"
    print "SigmaC:", SigmaC
    rhogas_2d = np.zeros((nphi_grid, nrad_grid), dtype=float)
    for j in range(nrad_grid):
        rr = rad_grid[j]
        isurf_1 = SigmaC * (rad_grid[j]/(AU*r_c))**(-gamma) * math.exp(-(rad_grid[j]/(AU*r_c))**(2.0-gamma))
        if min_number:
            surf_min = minimum_v 
        else:
            surf_min = isurf_1 * minimum
        isurf = isurf_1 * math.atan((rr/(AU*r_c_tan))**gamma_atan)  / 1.57079632679
        if (isurf < surf_min) and (rr/AU < r_c) and (rr/AU > innerdisk_r) :
            isurf = surf_min
        if inner_hole:
            if rr/AU < hole_size:
                isurf = 0
                print isurf
        if (rr/AU < innerdisk_r) :
            isurf = innerdisk_D
        rhogas_2d[:,j] = isurf
        surf1g[j] = isurf
    # calculate total disk mass
    np.savetxt('surf1g.dat', np.transpose([rad_grid/AU,surf1g]))
    disk_mass = 0.0
    for i in range(nrad_grid-1):
        disk_mass = disk_mass + math.pi*(rad_grid[i+1]**2.0-rad_grid[i]**2.0)*(rhogas_2d[0,i]+rhogas_2d[0,i+1])/2.0
        #disk_mass = disk_mass + 2.0*math.pi*((rad_grid[i]+rad_grid[i+1])/2.0)*(rad_grid[i+1]-rad_grid[i])*(rhogas_2d[0,i]+rhogas_2d[0,i+1])/2.0
    print "Disk mass in Mdisk, disk_mass:", Mdisk, disk_mass/Msun, "\n"
    # write the asccii file
    wfile = open("gas_surface_density.ascii", 'w')
    for i in range(nphi_grid):
        for j in range(nrad_grid):
            wfile.write('%.7e %.7e %.7e \n'%(phi_grid[i],rad_grid[j],rhogas_2d[i,j]))
    wfile.close()


#### write temp file

if write_temp:
    print "Write temperature:"
    temp_3d = np.zeros((nphi, nthet, nrad), dtype=float)

    for j in range(nrad):
        rr = (rad_grid[j]+rad_grid[j+1])/2.0
        for i in range(nthet):
            if math.tan(math.pi/2.0-(thet_grid[i]+thet_grid[i+1])/2.0) < tan_crit:
                #temp_3d[:,i,j] = 15.0*(rr*math.sin((thet_grid[i]+thet_grid[i+1])/2.0)/1.0/AU)**-0.5
                temp_3d[:,i,j] = 140.0*(rr*math.sin((thet_grid[i]+thet_grid[i+1])/2.0)/16.0/AU)**-0.5
            else:
                temp_3d[:,i,j] = 140.0*(rr*math.sin((thet_grid[i]+thet_grid[i+1])/2.0)/16.0/AU)**-0.5

    print "tranformed to 1 row array for write", "\n"
    temp_3d = np.array(temp_3d).reshape(nphi*nthet*nrad)

    with open('dust_temperature.bdat', 'wb') as f:
        f.write(struct.pack(4*'q',1,4,nphi*nthet*nrad,1))
        for num1 in temp_3d:
            f.write(struct.pack('f', num1))

