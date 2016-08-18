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

#fit_dust = False
fit_dust = True

#fit_gas  = True 
fit_gas  = False 
#write_gas3d_density  = True 
write_gas3d_density  = False

inner_hole = False
hole_size = 7.5  # cavity size in AU

write_temp = False
tan_crit = 0.1  # two-layer model temperature boundary

## Sigma_ = SigmaC * (r/RC)**-gamma * exp( -(r/RC)**(2-gamma)) * atan((r/RC)**gamma_atan)

if fit_gas:
    Mdisk = 1.2e-4         # gas disk mass in solar mass
    gamma = 0.2          #  power-law index
    r_c_tan   = 50.0        #  the radial size of the gas disk in AU
    r_c   = 150.0        #  the radial size of the gas disk in AU
    gamma_atan = 3
    minimum  = 0.10
    #
    gas2dust = 100
    min_number = True
    minimum_v  = 0.0005
    innerdisk_r = 5     # boundary of the inner disk
    innerdisk_hr_scale = 5  # hr_innerdisk/hr
    innerdisk_D = 0.03  # density of the inner disk

# fitted parameter for dust 
if fit_dust:
    Mdisk = 0.0076         # gas disk mass in solar mass
    # 0.0076  34.23       # gas disk mass in solar mass
    # 0.0075 chi 42.36
    # high 0.01
    # low 0.005
    gamma = -0.15 #  power-law index
    r_c   = 68.0        #  the radial size of the gas disk in AU
    r_c_tan   = 69.0        #  the radial size of the gas disk in AU
    gamma_atan = 4.35

    minimum  = 0.000
    #
    inner_cavity = 0
    #
    outerB = 1000 #AU
    outerD = 0.0 # g/cm^2

sett        = 0.1
h_2_r       = 0.05


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
if fit_dust:
    density_3d_1 = np.zeros((nphi,nthet,nrad), dtype=float)

surf1d = np.zeros(nrad, dtype=float)
surf1g = np.zeros(nrad_grid, dtype=float)


###################################################### 
######## manually set the density profile

print "write the gas surface density ascii", "\n"
if 1.0-gamma == -1:
    SigmaC = Mdisk*Msun/(2.0*math.pi)/((r_c*AU)**2.0)/(math.log(600.0/r_c)-math.log(5.0/r_c))
if 1.0-gamma != -1:
    SigmaC = Mdisk*Msun/(2.0*math.pi)/((r_c*AU)**2.0)/((600.0/r_c)**(2.0-gamma)/(2.0-gamma)-(5.0/r_c)**(2.0-gamma)/(2.0-gamma))


if fit_dust:
    for j in range(nthet):
        for k in range(nrad):
            rr = rad[k]*math.sin(thet[j])
            rz = rad[k]*math.cos(thet[j])
            hr = h_2_r*20*AU* (rr/20.0/AU)**1.25
            isurf_1 = SigmaC*(rad[k]/(AU*r_c))**(-gamma) * math.exp(-(rad[k]/(AU*r_c))**(2.0-gamma)) 
            surf_min = isurf_1 * minimum
            isurf = isurf_1 * math.atan((rr/(AU*r_c_tan))**gamma_atan) / 1.57079632679
            if isurf < surf_min:
                isurf = surf_min
            if (rr/AU < inner_cavity):
                isurf = 0
            if (rr/AU > outerB) :
                isurf = outerD
            #
            surf1d[k] = isurf
            rho0 = isurf/math.sqrt(math.pi*2.0)/(hr*sett)
            density_3d_1[:,j,k] = rho0*math.exp(-rz*rz/2.0/(hr*sett)/(hr*sett))

    np.savetxt('surf1d.dat', np.transpose([rad/AU,surf1d]))
    ########### write the new dust_density.binp
    print "Write the new density data"
    density_first = np.array(density_3d_1).reshape(nphi*nthet*nrad)
    print "tranformed to 1 row array for write"
    with open('dust_density.binp', 'wb') as f:
        f.write(struct.pack(4*'q',1,4,nphi*nthet*nrad,1))
        for num1 in density_first:
            f.write(struct.pack('f', num1))

