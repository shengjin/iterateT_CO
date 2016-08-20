#!/usr/bin/env python2
import sys
import os
import math
import linecache
import struct

import numpy as np

import matplotlib as mpl
mpl.use('Agg')   # to avoid the window in plotting

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogLocator

#sys.path.append('/turquoise/users/sjin/python/')

########################################################
############ Options and Constants

##### Ngrids in the z-dirction used in the cylindrical coordinate 
###     in the calculation of gas volume density
nr_cyl=250
nz_cyl=200
# scale: LIN or LOG
r_scale  = "LIN"
z_scale  = "LOG"
zmax2rmax = 0.5 # z_max = r_max * zmax2rmax

# Path to be created for gas line radmc
gas_radmc_path = "./gas_line_radmc"


# scale the gas surface density
scale_dens = 1   # scale the gas surface density


# create a spherical dustDens for iteratiive T_calculation?
interp_sphere_dustDens = False  
microdust2gas = 0.001  # dust2gas ratio for micron size dust for T calculation.

#moleName  = "12c16o" # the numberDens file will be numberdens_moleName.binp
moleName  = "12c18o" # the numberDens file will be numberdens_moleName.binp
dbase_sty = "leiden" # the datebase format of gas radiative properties 

#abundance = 1.4e-4   # molecule abundance in mass.  rho_mole = rho_H2 * abundance
#abundance = 3.0e-4   # molecule abundance in mass.  rho_mole = rho_H2 * abundance
                     # abundance(12CO) = 3.0e-4
#mu_molecule = 28     # for CO12, 12+16

abundance = 5.5e-7   # 13co/c180 = 5.5 (solar system)
    #abundance = 4.2e-5   # 13co/c180 = 5.5 (solar system)
mu_molecule = 30     # for C18O, 12+18

#abundance = 2.3e-6   # molecule abundance in mass.  rho_mole = rho_H2 * abundance
    #abundance = 2.3e-5   # molecule abundance in mass.  rho_mole = rho_H2 * abundance
                     # abundance(13CO) = 3.0e-4
#mu_molecule = 29     # for 13CO, 13+16


#abundance = 5e-5   # mass ratio of moleculae 

#gas_freeze  = False #True   # include the freeze out of gas or not
gas_freeze  = True #True   # include the freeze out of gas or not
freeze_T    = 20     # the freezing temperature of this molecule      

#photodissociation = False    # include photodissociation or not
photodissociation = True #True    # include photodissociation or not
#sigma_photodisso  = 0.003  # critical sigma of photodissociation, in g cm^-2

sigma_photodisso  = 3e-3  # critical sigma of photodissociation, in g cm^-2
sigma_photodisso_H = False  # surface density of H2

use_number = True         # use number_photodisso instead of sigma_photodisso
number_photodisso = 5e20  # number surf_dens of photodissociation, H used
fnumb_H    = 0.706        # number ratio of H nuclei
print sigma_photodisso
print number_photodisso 
print use_number

# Define some constants
AU2CM = 1.4959787e13   # Astronomical unit [cm]
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

temp3d_name = "%s%s" % ("dust_temperature", ".bdat")


# dens2d = dens2d.ascii *scale_dens
dens2d_name = "%s%s" % ("gas_surface_density", ".ascii")


##### make some plots, data_output for debug?
debug = False        # info printing
ddebug = False       # CAUTION: info printing inside loops
debug_plot = False    # make some plots
debug_data = False   # data output
debug_tmp = False      # tmp DEBUG

print ""
if debug:
    print "Using the gas temperatrue from: ", temp3d_name
    print "Using the gas surface density from: ", dens2d_name, "\n"

########### Function Definition

################################
# define the Func to transform cylindrical coordinates to spherical 
def cyl2sph(c_r,c_z):
    # input: r,z in cyl, output: r,thet in sph
    s_r = (c_r**2.0+c_z**2.0)**0.5
    s_thet = math.atan(c_r/c_z)
    return s_r,s_thet

################################
# define the Func to transform spherical coordinates to cylindrical
def sph2cyl(s_r,s_thet):
    # input: r,thet in sph, output: r,z in cyl
    c_r = s_r*math.sin(s_thet)
    c_z = s_r*math.cos(s_thet)
    return c_r,c_z

###############
###############
# define the interpolate function
def LinIntP(x,x1,y1,x2,y2):
    """
    makes linear interpolation for f(x) at the point x, if (x1,f(x1)=y1) and (x2,f(x2)=y2)
    using Lagrange's formula
    """
    return ((x-x2)/(x1-x2))*y1+((x-x1)/(x2-x1))*y2


#################################
############### read the grids

grid = linecache.getline('amr_grid.inp', 6)
ngrids = [int(s) for s in grid.split() if s.isdigit()]
nrad  = ngrids[0]
nthet = ngrids[1]
nphi  = ngrids[2]
if debug:
    print "amr_grid.inp grids: rad, thet, phi"
    print "nrad", nrad
    print "nthet", nthet
    print "nphi", nphi, "\n"

grids = np.genfromtxt('amr_grid.inp', skip_header=6, dtype=float)
if debug:
    print grids.shape
    print "Temp3d grids: rad, thet, phi"
    print "ngrids[0]", ngrids[0]
    print "ngrids[1]", ngrids[1]
    print "ngrids[2]", ngrids[2], "\n"

nrad_b1 = grids[0:nrad:1]
nrad_b2 = grids[1:nrad+1:1]
rad = (nrad_b1+nrad_b2)/2.0
rad0 = nrad_b1[0]
radE = nrad_b2[-1]
rad_grid = grids[0:nrad+1:1]
if debug:
    print "rad0, radE", rad0, radE
    print "shape of rad_grid", rad_grid.shape, "\n"

nthet_b1 = grids[nrad+1:nrad+1+nthet:1]
nthet_b2 = grids[nrad+2:nrad+2+nthet:1]
thet = (nthet_b1+nthet_b2)/2.0
thet0 = nthet_b1[0]
thetE = nthet_b2[-1]
thet_grid  = grids[nrad+1:nrad+2+nthet:1]
if debug:
    print "thet0, thetE", thet0, thetE
    print "shape of thet_grid", thet_grid.shape, "\n"

nphi_b1 = grids[nrad+2+nthet:nrad+2+nthet+nphi:1]
nphi_b2 = grids[nrad+3+nthet:nrad+3+nthet+nphi:1]
phi = (nphi_b1+nphi_b2)/2.0
phi0 = nphi_b1[0]
phiE = nphi_b2[-1]
phi_grid  = grids[nrad+2+nthet:nrad+3+nthet+nphi:1]
if debug:
    print "phi0, phiE", phi0, phiE, "\n"
    print "shape of phi_grid", phi_grid.shape, "\n"



#############################################
############### read the temperature data

head = np.fromfile(temp3d_name, count=4, dtype=int)
if debug:
    print "head[0]", head[0]
    print "head[1]", head[1]
    print "head[2]", head[2]
    print "head[3]", head[3]

if head[2]!=(nrad*nthet*nphi):
    print ' ERROR'
    print ' Number of grid points in '+temp3d_name+' is different from that in amr_grid.inp'
    quit()

print "Reading the 3d dust_temp file"

if head[1] == 8:
    f = open(temp3d_name, "rb")  # reopen the file
    f.seek(32, os.SEEK_SET)  # seek
    temp_3d = np.fromfile(f, count=-1, dtype=np.float64)
elif head[1]==4:
    f = open(temp3d_name, "rb")  # reopen the file
    f.seek(32, os.SEEK_SET)  # seek
    temp_3d = np.fromfile(f, count=-1, dtype=np.float32)
else:
    print 'ERROR'
    print 'Unknown datatype in '+temp3d_name
f.close()

print "dust_temp3d reading: done."

temp_3d = np.array(temp_3d).reshape(nphi,nthet,nrad)
print "temp_3d transformed to nphi,nthet,nrad", "\n"

if debug_plot:
    plt.pcolormesh(rad/AU2CM, thet, temp_3d[0,:,:])#, norm=LogNorm())
    plt.ylim([thet[nthet-1],thet[0]])
    #plt.pcolormesh(rad, phi, temp_3d[:,0,:])
    plt.xlabel('Radial Distance [AU]')
    plt.ylabel('Thet')
    ticks_1=[0,20,40,60,80,100,300]
    cbar1=plt.colorbar(ticks=ticks_1)
    cbar1.set_label('Temperature ( K ) ')
    plt.savefig('verticaltemp.png')
    plt.clf()


################################
##### read the 2d gas density 

print "Read the gas2d_dens file: ", dens2d_name, "\n"
phi_r_rho_gas = np.genfromtxt(dens2d_name, dtype=float)
phi_r_rho_gas[:,2] = phi_r_rho_gas[:,2]*scale_dens

if debug_plot:
    print "plot the 2d gas density"
    x_gas=phi_r_rho_gas[:,1]*np.cos(phi_r_rho_gas[:,0])/AU2CM
    y_gas=phi_r_rho_gas[:,1]*np.sin(phi_r_rho_gas[:,0])/AU2CM
    rho_gas=phi_r_rho_gas[:,2]
    #stokesnumber = (1.0/rho)*0.3*1.57079637051
    axisadd=-50
    plt.scatter(x_gas, y_gas, c=rho_gas, edgecolors='none', marker="o", s=(phi_r_rho_gas[:,1]/AU2CM/10), cmap="jet")
    plt.axis([x_gas.min()-axisadd, x_gas.max()+axisadd, y_gas.min()-axisadd, y_gas.max()+axisadd])
    plt.xlabel('Radial Distance [AU]')
    plt.ylabel('Radial Distance [AU]')
    ticks_at1 = [0,10,20,50,100,200]
    cbar2=plt.colorbar(ticks=ticks_at1)
    cbar2.set_label(r'$\mathrm{Surface \hspace{0.4} Density \hspace{0.4} [ \hspace{0.3}g/cm^2 \hspace{0.3}]}$', rotation=90, fontsize=14,  labelpad=0)
    plt.tight_layout()
    plt.savefig('gas2d_from_ascii.png')
    plt.clf()


#############################################
##### Create a Cylindrical coordinate system
##### to calulate the gas density ...

print "Create a cylindrical temperature 3d data-grid"

nr_cyl_grid = nr_cyl + 1
r_cyl_min = rad0
r_cyl_max = radE
if debug:
    print "r_cyl_min, r_cyl_max", r_cyl_min, r_cyl_max
    print "r_scale: ", r_scale, "\n"
#  // define the radial grid
r_cyl_grid = np.zeros(nr_cyl_grid)
if (r_scale == "LIN"):
    dr = (r_cyl_max-r_cyl_min)/(nr_cyl)
    for i in range(nr_cyl_grid):
        r_cyl_grid[i] = r_cyl_min  + i*dr
if (r_scale == "LOG"):
    dr = (r_cyl_max-r_cyl_min)**(1./(nr_cyl))
    for i in range(nr_cyl_grid):
        r_cyl_grid[i] = r_cyl_min*dr**i
if debug_data:
    np.savetxt('r_cyl_grid.out', np.transpose(r_cyl_grid))


nz_cyl_grid = nz_cyl + 1
z_cyl_min = 0.0
z_cyl_max = zmax2rmax*radE
#z_cyl_max = radE*np.cos(thet0)
if debug:
    print "z_cyl_min, z_cyl_max", z_cyl_min, z_cyl_max
    print "z_scale: ", z_scale, "\n"
z_cyl_grid = np.zeros(nz_cyl_grid)
if (z_scale == "LIN"):
    dz = (z_cyl_max-z_cyl_min)/(nz_cyl)
    for i in range(nz_cyl_grid):
        z_cyl_grid[i] = z_cyl_min + i*dz
if (z_scale == "LOG"):
    dz = (z_cyl_max-z_cyl_min)**(1./(nz_cyl))
    for i in range(nz_cyl_grid):
        z_cyl_grid[nz_cyl-i] = z_cyl_max/dz**i
    z_cyl_grid[0] = z_cyl_min
if debug_data:
    np.savetxt('z_cyl_grid.out', np.transpose(z_cyl_grid))


########################
### create a new dir for line transfer

if not os.path.exists(gas_radmc_path):
    os.mkdir(gas_radmc_path, 0755)
    print gas_radmc_path, " is created", "\n"


################################
# create the new temp3d in cylindrial coordinates
print "Create the new temp3d in cylindrial coordinates"
temp_cyl_3d = np.zeros((nr_cyl,nphi,nz_cyl), dtype=float)

# create the r_cyl and z_cyl
r_cyl = np.zeros(nr_cyl)
z_cyl = np.zeros(nz_cyl)
for i in range(nr_cyl):
    r_cyl[i] = (r_cyl_grid[i]+r_cyl_grid[i+1])/2.0
for i in range(nz_cyl):
    z_cyl[i] = (z_cyl_grid[i]+z_cyl_grid[i+1])/2.0
if debug_data:
    np.savetxt('r_cyl.out', np.transpose(r_cyl))
    np.savetxt('z_cyl.out', np.transpose(z_cyl))

# Calculate temp_cyl_3d by shooting.
print "Calculate temp_cyl_3d by shooting:"
for j in range(nr_cyl):
    for k in range(nz_cyl):
        s_r, s_thet = cyl2sph(r_cyl[j],z_cyl[k])
        nshoot_rad=0
        nshoot_thet=0
        for i in range(nrad):
            if s_r > rad_grid[nrad]:
                nshoot_rad = 9999
                if ddebug:
                    print nshoot_rad
                break
            elif (s_r > rad_grid[i]) and (s_r < rad_grid[i+1]):
                nshoot_rad = i
                if ddebug:
                    print nshoot_rad
                break
        for i in range(nthet):
            if s_thet < thet_grid[0]:
                nshoot_thet = 9999
                if ddebug:
                    print nshoot_rad
                break
            elif (s_thet > thet_grid[nthet]) and (s_thet < thet_grid[nthet]*1.001):
                nshoot_thet = nthet-1
                if ddebug:
                    print nshoot_rad
                break
            elif (s_thet > thet_grid[i]) and (s_thet < thet_grid[i+1]):
                nshoot_thet = i
                if ddebug:
                    print nshoot_rad
                break
        if (nshoot_rad == 9999) or (nshoot_thet == 9999):
            temp_cyl_3d[j,:,k] = 0.0
        else:
            temp_cyl_3d[j,:,k] = temp_3d[:,nshoot_thet,nshoot_rad]
print "Calculate temp_cyl_3d by shooting. DONE", "\n"

if debug_plot:
    #plt.pcolormesh(r_cyl/AU2CM, z_cyl/AU2CM, temp_cyl_3d[:,0,:].T, vmin=10, vmax=200, norm=LogNorm())
    plt.pcolormesh(r_cyl/AU2CM, z_cyl/AU2CM, temp_cyl_3d[:,0,:].T, vmin=10, vmax=temp_cyl_3d[:,0,:].max(), norm=LogNorm())
    #plt.ylim([thet[nthet-1],thet[0]])
    #plt.pcolormesh(rad, phi, temp_3d[:,0,:])
    plt.xlabel('Radial Distance [AU]')
    plt.ylabel('Z [AU]')
    # axes[i,j].set_xticklabels(['-600', '-400', '-200', '0', '200', '400', '600'],rotation=90)
    ticks_1=[10,20,100,200,500]
    cbar3=plt.colorbar(ticks=ticks_1)
    #cbar3=plt.colorbar(ticks = LogLocator(subs=range(10)))  
    cbar3.ax.set_yticklabels(["10", "20", "100", "200", "500"])
    cbar3.set_label('Temperature ( K ) ')
    plt.savefig('verticaltemp_cyl.png')
    plt.clf()


################################
# 2D (nrad+1,nphi+1)  gas_surf from gas_surface_density.ascii

gas_surf = np.zeros((nrad+1,nphi+1), dtype=float)
gas_surf = np.array(phi_r_rho_gas[:,2]).reshape(nphi+1,nrad+1).T

### output the 2D gas_surf to check the trasform is right
if debug_plot:
    plt.pcolormesh(rad_grid/AU2CM, phi_grid, gas_surf.T)
    plt.colorbar()
    plt.savefig('gas_2d_from_ascii_tranform.png')
    plt.clf()

gas_surf_cyl = np.zeros((nr_cyl,nphi), dtype=float)

print "Calculate gas_surf_cyl by interpolation:"
for i in range(nr_cyl):
    rangeout=False
    ilow=0
    iup=0
    #
    for ii in range(nrad):
        if r_cyl[i] < rad_grid[0]:
            rangeout=True
            break
        elif r_cyl[i] > rad_grid[nrad]:
            rangeout=True
            break
        elif (r_cyl[i] < rad_grid[ii+1]):
        #elif (r_cyl[i] > rad_grid[ii]) and (r_cyl[i] < rad_grid[ii+1]):
            iup=ii+1
            if (r_cyl[i] > rad_grid[ii]):
                ilow=ii
                break
            else:
                if (r_cyl[i] > rad_grid[ii-1]):
                    ilow=ii
                    break
                else:
                    print "WRONG"
                quit()
    #
    for j in range(nphi):
        if rangeout == True:
            gas_surf_cyl[i,j] = 0
        else:
            rho_up = (gas_surf[iup,j]+gas_surf[iup,j+1])/2.0
            rho_low = (gas_surf[ilow,j]+gas_surf[ilow,j+1])/2.0
            gas_surf_cyl[i,j] = LinIntP(r_cyl[i],rad_grid[ilow],rho_low,rad_grid[iup],rho_up)
print "Calculate gas_surf_cyl by interpolation. DONE", "\n"

if debug_plot:
    phi_cyl = np.zeros(nphi)
    for i in range(nphi):
        phi_cyl[i] = (phi_grid[i]+phi_grid[i+1])/2.0
    plt.pcolormesh(r_cyl/AU2CM, phi_cyl, gas_surf_cyl.T)
    plt.colorbar()
    plt.savefig('gas_2d_from_ascii_tranform_interpl.png')
    plt.clf()



################################################ 
##### Solve the volume density of gas using 
####   Eq. 12 in Rosenfeld et al. 2013

# create the gas3d in cylindrial coordinates
print "Create the gas_rho_3d in cylindrial coordinates"
gas_rho_3d = np.zeros((nr_cyl,nphi,nz_cyl), dtype=float)


# physical quantities needed for solving the differential equation
rho_z = np.zeros(nz_cyl)
T_z   = np.zeros(nz_cyl)
#lnT   = np.zeros(nz_cyl)
cs_z  = np.zeros(nz_cyl)

z_len = np.zeros(nz_cyl)
for i in range(nz_cyl):
    z_len[i] = z_cyl_grid[i+1] - z_cyl_grid[i]

dz    = np.zeros(nz_cyl-1)
for i in range(nz_cyl-1):
    dz[i] = (z_cyl_grid[i+2] - z_cyl_grid[i])/2

print "Solve the differential equation for hydro-equilibrium:"
GM = GG*Mstar*Msun
for i in range(nr_cyl):
    # print a progress bar
    if(i%max(1,int(0.01*nr_cyl))==0):
        print "%s%s" % (100*i/nr_cyl, "% finished ...")
    #
    r_square = r_cyl[i]**2.0
    #
    rho_z = rho_z*0.0
    for j in range(nphi):
        gas_surf_once = gas_surf_cyl[i,j]
        T_z  = temp_cyl_3d[i,j,:]
        cs_z = (T_z*kB/mu/mh)**0.5
        #
        rho_z[0] = 100.0
        if T_z[0] == 0:
            print i, j, "(Tz0 = 0)"
            break
        for k in range(nz_cyl-1):
            if T_z[k+1] == 0:
                rho_z[k+1] = 0
                gas_rho_3d[i,j,k] = rho_z[k+1]
            else:
                del_lnT = np.log(T_z[k+1]) - np.log(T_z[k])
                right_1   = del_lnT/dz[k]
                right_2   = GM*z_cyl[k+1]/(r_square+z_cyl[k+1]**2.0)**1.5*1.0/(cs_z[k+1]**2.0)
                ln_rho    = -(right_1+right_2)*dz[k] + np.log(rho_z[k])
                rho_z[k+1] =  np.exp(ln_rho)
                gas_rho_3d[i,j,k] = rho_z[k+1] 
print "Solve the differential equation for hydro-equilibrium. DONE", "\n"
        
# write a vertical slice of gas density
if debug_data:
    np.savetxt('rho_z_vertical_1d.out', np.transpose(gas_rho_3d[1,0,:]))
# plot a vertical slice of the 3d gas density.
if debug_plot:
    plt.pcolormesh(r_cyl/AU2CM, z_cyl/AU2CM, gas_rho_3d[:,0,:].T)
    #plt.ylim([thet[nthet-1],thet[0]])
    #plt.pcolormesh(rad, phi, temp_3d[:,0,:])
    plt.xlabel('Radial Distance [AU]')
    plt.colorbar()
    plt.ylabel('Z [AU]')
    plt.savefig('vertical_gas_rho_cyl.png')
    plt.clf()


############################################################
####### Normalization  of gas density based on gas_surf_cyl

## the scaling factor 
scaling_dens = np.zeros((nr_cyl,nphi), dtype=float)

print "Scaling the gas density based the hydro surface density:"
#########
# calculate the scaling factor
for i in range(nr_cyl):
    for j in range(nphi):
        rho_surf = 0.0
        for k in range(nz_cyl):
            rho_surf = rho_surf + gas_rho_3d[i,j,k]*z_len[k]
        scaling_dens[i,j] = gas_surf_cyl[i,j]/rho_surf

# scaling the gas density using scaling factor
for i in range(nr_cyl):
    for j in range(nphi):
        gas_rho_3d[i,j,:] = gas_rho_3d[i,j,:]*scaling_dens[i,j]*0.5
print "Scaling the gas density based the hydro surface density. DONE", "\n"

### plot a vertical slice of the new gas_rho_3d
if debug_plot:
    #plt.plot(z_cyl/AU2CM, gas_rho_3d[15,0,:])
    plt.loglog(z_cyl/AU2CM, gas_rho_3d[15,0,:])
    #plt.pcolormesh(r_cyl/AU2CM, z_cyl/AU2CM, gas_rho_3d[:,0,:].T, vmin=1e-5, vmax=gas_rho_3d[:,0,:].max(), norm=LogNorm())
    plt.ylabel('Gas Density ( g/cm^3 ) ')
    plt.xlabel('Z [AU]')
    plt.savefig('vertical_gas_rho_cyl_scaled.png')
    plt.clf()

### Re-Calculate the sur_face density and make a plot
if debug_plot:
    print "Re-Calculate the sur_face density and make a plot.", "\n"
    for i in range(nr_cyl):
        for j in range(nphi):
            rho_surf = 0.0
            for k in range(nz_cyl):
                rho_surf = rho_surf + gas_rho_3d[i,j,k]*z_len[k]
            gas_surf_cyl[i,j] = rho_surf
    # make a figure
    phi_cyl = np.zeros(nphi)
    for i in range(nphi):
        phi_cyl[i] = (phi_grid[i]+phi_grid[i+1])/2.0
    plt.pcolormesh(r_cyl/AU2CM, phi_cyl, gas_surf_cyl.T)
    plt.colorbar()
    plt.savefig('gas_2d_from_ascii_tranform_regenerated_integrated3D.png')
    plt.clf()


###################################################
####### Calculate the number density by gas_rho_3d(nr_cyl,nphi,nz_cyl)
#######      + temp_cyl_3d(nr_cyl,nphi,nz_cyl)

print "Convert the volume density of gas to number density:"
gas_number_3d = np.zeros((nr_cyl,nphi,nz_cyl), dtype=float)

## tweak the number density to account for the freeze out of gas
if gas_freeze:
    print "gas were freezed out at region with T < ", freeze_T, " K. ..."
    for i in range(nr_cyl):
        for j in range(nphi):
            for k in range(nz_cyl):
                if temp_cyl_3d[i,j,k] < freeze_T:
                    gas_rho_3d[i,j,k] = 0.0
                    if ddebug:
                        print i,j,k,temp_cyl_3d[i,j,k],"Freeze"
    if debug_plot:
        #plt.plot(z_cyl/AU2CM, gas_rho_3d[15,0,:])
        plt.loglog(z_cyl/AU2CM, gas_rho_3d[15,0,:])
        plt.ylabel('Gas Density ( g/cm^3 ) ')
        plt.xlabel('Z [AU]')
        plt.savefig('vertical_gas_rho_cyl_scaled_freezed.png')
        plt.clf()

## tweak the number density to account for photodissociation
if photodissociation:
    print "including photodissociation... "
    ## find the layer that above which gas is dissociated
    nlayer_disso = np.zeros((nr_cyl,nphi), dtype=float)
    nlayer_disso_ratio = np.zeros((nr_cyl,nphi), dtype=float)
    if use_number:
        for i in range(nr_cyl):
            for j in range(nphi):
                nlayer = 0.0
                number_tmp = 0.0
                for k in range(nz_cyl):
                    k_rv = nz_cyl-1-k
                    surfdens_once = gas_rho_3d[i,j,k_rv]*z_len[k_rv]
                    number_once = (surfdens_once/mu/mh)*fnumb_H
                    number_tmp = number_tmp + number_once
                    if number_tmp > number_photodisso:
                        nlayer = k_rv
                        break
                nlayer_disso[i,j] = nlayer
                nlayer_disso_ratio[i,j] = float(nz_cyl-nlayer)/float(nz_cyl)
    else:
        for i in range(nr_cyl):
            for j in range(nphi):
                nlayer = 0.0
                sigma_tmp = 0.0
                for k in range(nz_cyl):
                    k_rv = nz_cyl-1-k
                    if sigma_photodisso_H:
                        sigma_tmp = sigma_tmp + gas_rho_3d[i,j,k_rv]*z_len[k_rv]
                    else:
                        sigma_tmp = sigma_tmp + gas_rho_3d[i,j,k_rv]*z_len[k_rv]*abundance
                    if sigma_tmp > sigma_photodisso:
                        nlayer = k_rv
                        break
                nlayer_disso[i,j] = nlayer
                nlayer_disso_ratio[i,j] = float(nz_cyl-nlayer)/float(nz_cyl)
    # plot the nlayer
    if debug_plot:
        plt.pcolormesh(r_cyl/AU2CM, phi_cyl, nlayer_disso_ratio.T) 
        plt.xlabel('Radial Distance [AU]')
        plt.ylabel('phi')
        cbar5=plt.colorbar()
        cbar5.set_label('percent of top layers photodiso... ed')
        plt.savefig('nlayer_photodisso.png')
        plt.clf()
    #
    # dissociate the gas above nlayer_disso
    for i in range(nr_cyl):
        for j in range(nphi):
            for k in range(nz_cyl):
                k_rv = nz_cyl-1-k
                if k_rv > nlayer_disso[i,j]:
                    gas_rho_3d[i,j,k_rv] = 0.0
                    if ddebug:
                        print gas_rho_3d[i,j,k_rv], i,j,k 
                else:
                    break
    ### plot a vertical slice of the new gas_rho_3d
    if debug_plot:
        #plt.plot(z_cyl/AU2CM, gas_rho_3d[15,0,:])
        plt.loglog(z_cyl/AU2CM, gas_rho_3d[15,0,:])
        plt.ylabel('Gas Density ( g/cm^3 ) ')
        plt.xlabel('Z [AU]')
        plt.savefig('vertical_gas_rho_cyl_scaled_photodissoed.png')
        plt.clf()


### Convert the volume densit to number densty
gas_number_3d = gas_rho_3d*abundance/mu_molecule/mh

if debug_plot:
    #plt.plot(z_cyl/AU2CM, gas_number_3d[15,0,:])
    plt.loglog(z_cyl/AU2CM, gas_number_3d[15,0,:])
    plt.ylabel('Number Density')
    plt.xlabel('Z [AU]')
    plt.savefig('vertical_gas_rho_cyl_final_nuberDens.png')
    plt.clf()
    #plt.pcolormesh(z_cyl/AU2CM, r_cyl/AU2CM, gas_number_3d[:,1,:].T, norm=LogNorm())
    plt.pcolormesh(r_cyl/AU2CM, z_cyl/AU2CM, gas_number_3d[:,0,:].T,  norm=LogNorm())
    plt.xlabel('AU')
    plt.ylabel('Z [AU]')
    ticks_1=[10,100,1000,10000]
    cbar1=plt.colorbar(ticks=ticks_1)
    cbar1.set_label('Number density ( 1 / cm^-3 ) ')
    plt.savefig('2d_vertical_gas_rho_cyl_final_nuberDens.png')
    plt.clf()

print "Convert the volume density of gas to number density. DONE", "\n"


################################
#### convert the numberdens_... to spherical co-ordt

gas_number_sph_3d = np.zeros((nphi,nthet,nrad), dtype=float)

# Calculate gas_number_3d_sph by shooting.
print "Calculate gas_number_sph_3d by shooting:"
for j in range(nthet):
    for k in range(nrad):
        c_r, c_z = sph2cyl(rad[k],thet[j])
        nshoot_r=0
        nshoot_z=0
        for i in range(nr_cyl):
            if c_r < r_cyl_grid[0]:
                nshoot_r = 9999
                if ddebug:
                    print nshoot_r
                break
            elif (c_r > r_cyl_grid[i]) and (c_r < r_cyl_grid[i+1]):
                nshoot_r = i
                if ddebug:
                    print nshoot_r
                break
        for i in range(nz_cyl):
            if c_z > z_cyl_grid[nz_cyl]:
                nshoot_z = 9999
                if ddebug:
                    print nshoot_z
                break
            elif (c_z < z_cyl_grid[0]) and (c_z > -z_cyl_grid[1]):
                nshoot_z = 0
                if ddebug:
                    print nshoot_z
                break
            elif (c_z > z_cyl_grid[i]) and (c_z < z_cyl_grid[i+1]):
                nshoot_z = i
                if ddebug:
                    print nshoot_z
                break
        if (nshoot_r == 9999) or (nshoot_z == 9999):
            gas_number_sph_3d[:,j,k] = 0.0
        else:
            gas_number_sph_3d[:,j,k] = gas_number_3d[nshoot_r,:,nshoot_z]
if debug_plot:
    plt.pcolormesh(rad/AU2CM, thet, gas_number_sph_3d[0,:,:]) #, norm=LogNorm())
    #plt.pcolormesh(rad/AU2CM, thet, gas_number_sph_3d[1,:,:], norm=LogNorm())
    plt.xlabel('AU')
    plt.ylim(thetE,thet0)
    plt.ylabel('polar Deg')
    plt.colorbar()
    plt.savefig('2d_vertical_gas_rho_sph_final_nuberDens_shooted.png')
    plt.clf()

print "Calculate gas_number_shp_3d by shooting. DONE", "\n"


################################
#### write the numberdens_moleName.binp file 

numberdens_name = "%s%s%s" % ("numberdens_", moleName, ".binp")

print "Write the numberdens file: ",  numberdens_name
numberdens_1Darray = np.array(gas_number_sph_3d).reshape(nphi*nthet*nrad)
numberdens_file = "%s%s%s" % (gas_radmc_path, "/", numberdens_name)
with open(numberdens_file, 'wb') as f:
    f.write(struct.pack(3*'q',1,4,nphi*nthet*nrad))
    for num1 in numberdens_1Darray:
        f.write(struct.pack('f', num1))
print "Write the numberdens file. DONE", "\n"


#############################
##### create a spherical dustDens before molecule freeze 
##    and photodissociation for iteratiive T_calculation.
### write a new density file for micro-size density

# NOTE: turn off freeze and dissosication
if interp_sphere_dustDens:
    mu_dust_density_from_numberdens = "%s%s" % ("./dust_density", ".binp")
    convert_factor = abundance/mu_molecule/mh
    micro_dust_density = gas_number_sph_3d/convert_factor*microdust2gas
    print "Write the dust_density_micron-size file: ", mu_dust_density_from_numberdens 
    newdens_1Darray = np.array(micro_dust_density).reshape(nphi*nthet*nrad)
    newdens_file = "%s%s" % ("./", mu_dust_density_from_numberdens)
    with open(newdens_file, 'wb') as f:
        f.write(struct.pack(4*'q',1,4,nphi*nthet*nrad,1))
        for num1 in newdens_1Darray:
            f.write(struct.pack('f', num1))
    print "Write the numberdens file. DONE", "\n"



################################################
####### write the lines.inp

print "write the lines.inp for line-transfer", "\n"
lines_inp_name = "%s%s" % (gas_radmc_path, "/lines.inp")
wfile = open(lines_inp_name, 'w')
wfile.write("%d\n"%2)        # File format
wfile.write("%d\n"%1)       # Nr of gas species
wfile.write("%s %s %d %d %d\n"%(moleName, dbase_sty, 0, 0, 0))
wfile.close()

print "Done"
