#!/usr/bin/env python2
import numpy as np
import math
import os
import struct

from matplotlib.colors import LogNorm
import matplotlib as mpl
mpl.use('Agg')   # to avoid the window

#then
import matplotlib.pyplot as plt


## number of grids
nrad=150
#nrad=200
nthet=30
nphi=1


rad_min = 5.0
rad_max = 600.0
scaling_radius = 1.0
lograd = False #True

logthet = True
thet_min = 0.6981317008
#thet_min = 1.0472 #0.6981317008
thet_max = np.pi/2.0
# dthet  = (thet_max-thet_min)/pow(logthetweight,nthet)
logthetweight = 1.15

phi_min = 0.0
phi_max = np.pi*2.0

dens2d_name = "%s%s" % ("gas_surface_density", ".ascii")

write_amr = True
write_temp = False
write_vgas = True
write_rhogas = False
write_turbvgas = False

gas_radmc_path = "./gas_line_radmc"

turbvgas  = 1000.0   # turbulent velocity of gas in cm/s

## Sigma_ = SigmaC * (r/RC)**-gamma * exp( -(r/RC)**(2-gamma))
Mdisk = 0.12         # gas disk mass in solar mass
r_c   = 115.0        #  the radial size of the gas disk in AU
gamma = 0.8          #  power-law index

abundance = 5.0e-5   # molecule abundance.  rho_mole = rho_H2 * abundance
                     # abundance(12CO) = 3.0e-4
mu_molecule = 28     # for CO12, 12+16

moleName  = "12c16o" # the numberDens file will be numberdens_moleName.binp
dbase_sty = "leiden" # the datebase format of gas radiative properties 

gas_freeze  = True   # include the freeze out of gas or not
freeze_T    = 20     # the freezing temperature of this molecule      

photodissociation = True   # include photodissociation or not
sigma_photodisso  = 0.003  # critical sigma of photodissociation, in g cm^-2

# Define some constants
AU    = 1.4959787e13   # Astronomical unit [cm]
GG    = 6.67259e-8     # Gravitational constant
Msun  = 1.9891e33      # Solar mass [g]
Mstar = 0.97           # Stellar mass in Msun
kB    = 1.380658e-16   # Bolzmann's constant [erg/K]
mh    = 1.67262158e-24 # Mass of proton [g]

Y = 0.24
X = 1.0-Y
mu1   = 1.0/(X/2.006+Y/4.008)
mu2   = 2.3000e0   # Mean molecular weight (H2 + He + metals)
mu    = mu2


############################

## number of grids
nrad_grid=nrad+1
nthet_grid=nthet+1
nphi_grid=nphi+1

rin  = rad_min  * scaling_radius * AU
rout = rad_max  * scaling_radius * AU
#  // define the radial grid

rad_grid = np.zeros(nrad_grid)
if lograd:
    dr = (rout/rin)**(1./(nrad))
    for i in range(nrad_grid):
        rad_grid[i] = rin*dr**i
else:
    dr = (rout-rin)/(nrad) 
    for i in range(nrad_grid):
        rad_grid[i] = rin  + i*dr
      
      
thet_grid = np.zeros(nthet_grid)
if logthet:
    dthet = (thet_max-thet_min)/(logthetweight**nthet)
    for i in range(nthet_grid):
        thet_grid[i] = thet_max - dthet*(logthetweight**(nthet-i))
    thet_grid[nthet] = thet_max 
else:
    dthet  = (thet_max-thet_min)/(nthet)
    for i in range(nthet_grid):
        thet_grid[i] = thet_min+i*dthet


dphi  = (phi_max-phi_min)/(nphi)
phi_grid = np.zeros(nphi_grid)
for i in range(nphi_grid):
    phi_grid[i] = phi_min+i*dphi


print ""
print "nrad", nrad
print "lograd", lograd
print "nthet", nthet
print "logthet", logthet
print "nphi", nphi
print "nrad_grid", nrad_grid
print "nthet_grid", nthet_grid
print "nphi_grid", nphi_grid, "\n"


##################################
#### write amr_grid file
if write_amr:
    fname = 'amr_grid.inp'
    print 'Writing '+fname
    wfile = open(fname, 'w')
    wfile.write('%d\n'%1)     # Format number
    wfile.write('%d\n'%0)     # AMR self.style (0=regular self. NO AMR)
    wfile.write('%d\n'%100)   # Coordinate system (0-99 cartesian, 100-199 spherical)
    wfile.write('%d\n'%0)     # Gridinfo
    wfile.write('%d %d %d \n'%(1, 1, 0))             # Which dimension is active
    wfile.write('%d %d %d \n'%(nrad, nthet, nphi))  # Grid size (x,y,z or r,phi,thet)
    for i in range(nrad_grid): wfile.write('%.9e\n'%rad_grid[i])
    for i in range(nthet_grid): wfile.write('%.9e\n'%thet_grid[i])
    for i in range(nphi_grid): wfile.write('%.9e\n'%phi_grid[i])
    wfile.close()


##################################
#### write temp file

if write_temp:
    print "Write temperature:"
    temp_3d = np.zeros((nphi, nthet, nrad), dtype=float)

    for j in range(nrad):
        rr = (rad_grid[j]+rad_grid[j+1])/2.0
        for i in range(nthet):
            temp_3d[:,i,j] = 65.0*(rr*math.sin((thet_grid[i]+thet_grid[i+1])/2.0)/100.0/AU)**-0.5

    print "tranformed to 1 row array for write", "\n"
    temp_3d = np.array(temp_3d).reshape(nphi*nthet*nrad)

    with open('dust_temperature.bdat', 'wb') as f:
        f.write(struct.pack(4*'q',1,4,nphi*nthet*nrad,1))
        for num1 in temp_3d:
            f.write(struct.pack('f', num1))


##################################
#### write gas velocity

if write_vgas:
    ### create a new dir for line transfer
    if not os.path.exists(gas_radmc_path):
        os.mkdir(gas_radmc_path, 0755)
        print gas_radmc_path, " is created", "\n"
    print "write the gas_velocity.inp"
    vgas_3d = np.zeros((nphi, nthet, nrad, 3), dtype=float)
    for j in range(nrad):
        rr = (rad_grid[j]+rad_grid[j+1])/2.0
        for i in range(nthet):
            thet = (thet_grid[i]+thet_grid[i+1])/2.0
            rrr = rr * math.sin(thet)
            vgas_3d[:,i,j,2] = (GG*Mstar*Msun/rrr)**0.5
    gasV_name = "%s%s" % (gas_radmc_path, "/gas_velocity.inp")
    wfile = open(gasV_name, 'w')
    wfile.write('%d\n'%1)     # Format number
    wfile.write('%d \n'%(nrad*nthet*nphi))  # Grid size (x,y,z or r,phi,thet)
    for i in range(nphi):
        for j in range(nthet):
            for k in range(nrad):
                wfile.write('%.7e %.7e %.7e \n'%(vgas_3d[i,j,k,0],vgas_3d[i,j,k,1],vgas_3d[i,j,k,2]))
    wfile.close()
if write_turbvgas:
    turbV_name = "%s%s" % (gas_radmc_path, "/microturbulence.inp")
    wfile = open(turbV_name, 'w')
    wfile.write('%d\n'%1)     # Format number
    wfile.write('%d \n'%(nrad*nthet*nphi))  # Grid size (x,y,z or r,phi,thet)
    for i in range(nphi):
        for j in range(nthet):
            for k in range(nrad):
                wfile.write('%.7e \n'%(turbvgas))
    wfile.close()
    ### create a new dir for line transfer

##################################
#### creat gas density ascii

if write_rhogas:
    print "write the gas surface density ascii", "\n"
    SigmaC = Mdisk*Msun*(2.0-gamma)/(2.0*math.pi*(r_c*AU)**2.0)
    print "SigmaC:", SigmaC
    rhogas_2d = np.zeros((nphi_grid, nrad_grid), dtype=float)
    for i in range(nphi_grid):
        for j in range(nrad_grid):
            rhogas_2d[i,j] = SigmaC * (rad_grid[j]/(AU*r_c))**(-gamma) * math.exp(-(rad_grid[j]/(AU*r_c))**(2.0-gamma))
    # calculate total disk mass
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

