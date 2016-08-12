#!/usr/bin/env python2

import os
import numpy as np

import matplotlib as mpl
mpl.use('Agg')  # to avoid window

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
#from matplotlib.colors import SymLogNorm

from readimage import *
from natconst import *
from matplotlib import cm
import math
from math import sin, cos


remove_dust = True

# set it to True if an new image.out just has been generated
#newfile = False
newfile = True

RST_freq = 3.457959900000E+11
freq_delt = 2.422247656860E+05

#RST_freq = 3.305912230320E+11
#freq_delt = 2.315717818604E+05

#RST_freq = 3.293335270120E+11
#freq_delt = 6.591170815430E+05

cc = 299792.458  #km.s^-1
dv = freq_delt/RST_freq*cc

# inclination of the disk
incli = 51.0
#incli = 46.72

# position angle
# use it to scale the center Radius VS T_B
posang = 151.0
#posang = 138.02

# rotation angle used in the center ellipse to
rotang = posang
#rotang = 90-posang

# radius_wanted
n_rad = 60
d_rad = 10
inneravoid = 0
gap_width = 10.0 # AU

dpc = 140.0  # distance in pc

remove_tiny = False
tiny = 1e-30

# LkCa 15
# dust
#gaus_ma = 6.323864890469E-05    # maj axis of guassian in deg
#gaus_mi = 4.604819748137E-05    # min
#gaus_pa = 2.383443450928E+01    # PA of maj axis of Guass measure East from North

# 12CO
gaus_ma = 9.943533274863E-05                                                  
gaus_mi = 6.375448157390E-05
gaus_pa =  3.046934509277E+01                    

# 13CO
#gaus_ma  =   7.848847243521E-05                                                  
#gaus_mi  =   5.940155022674E-05                                                  
#gaus_pa  =   2.329190826416E+01 

# C18O
#gaus_ma  =   8.189328842693E-05                                                  
#gaus_mi  =   6.247880558173E-05                                                  
#gaus_pa  =   2.386120986938E+01 

# HL Tau
#gaus_ma = 8.367259158856E-06*3.5
#gaus_mi = 5.335457002123E-06*3.5
#gaus_pa = -1.758292236328E+02 

gaus_ma_arcsec = gaus_ma*3600
gaus_mi_arcsec = gaus_mi*3600
if gaus_pa > 0:
    gaus_pa_ctclk  = 180.0-gaus_pa  # PA in imConv
else:
    gaus_pa_ctclk  = -gaus_pa  # PA in imConv

fwhm = [gaus_ma_arcsec, gaus_mi_arcsec]

#debug = True
debug = False

AU = 1.496e13

# the n_th image for azimuthal extraction
n_image = 0

#########################
#### print the gaussian
print "gaus_pa_ctclk"
print gaus_pa_ctclk
print "gaus_ma_arcsec, gaus_mi_arcsec"
print gaus_ma_arcsec, gaus_mi_arcsec
gaus_ma_rad = gaus_ma*math.pi/180.0
gaus_mi_rad = gaus_mi*math.pi/180.0
gaus_ma_AU = math.sin(gaus_ma_rad)*dpc*pc/AU
gaus_mi_AU = math.sin(gaus_mi_rad)*dpc*pc/AU
print "gaus_ma_AU, gaus_mi_AU"
print gaus_ma_AU, gaus_mi_AU


#debug = True
debug = False

colorlog = False

###################### ###################### 
########### Read image and scale the flux
###################### ###################### 

if newfile:
    os.system('mv image.out image_gas.out')
    #os.system('mv image_gas.out image_gas.out.old')

#image = readImage('./image.out')
image = readImage('./image_gas.out')
dpc = 140.0  # distance in pc

#################################
# NOTE: in this case, nx=ny=m=n
# assuming the same number of grids in x,y axes
nx = image.nx
NN = nx
pix_x = np.linspace(-NN/2,NN/2,nx, endpoint=False)
ny = image.ny
pix_y = np.linspace(-NN/2,NN/2,ny, endpoint=False)

im_x_au = (pix_x)*image.sizepix_x/au
im_y_au = (pix_y)*image.sizepix_y/au


######################
print "INFO:","\n"
print "image_dust.nwav"
print image.nwav
print "image.wav"
print image.wav


#######################################################
#####  Remove the dust continuum at each wav for gas image
#########################################

if remove_dust:
    nwav = image.nwav
    for i in range(image.nwav):
        imagedust = "%s%s" % ("image.out_lambda_", 0)
        image_dust = readImage(imagedust)
        image.image[:,:,i] = image.image[:,:,i] - image_dust.image[:,:,0]
        image.imageJyppix[:,:,i] = image.imageJyppix[:,:,i] - image_dust.imageJyppix[:,:,0]



#######################################################
#####  Create average image
#########################################

image_ave = readImage('./image.out_lambda_0')

image_ave.image[:,:,0] = image_ave.image[:,:,0]*0.0
image_ave.imageJyppix[:,:,0] = image_ave.imageJyppix[:,:,0]*0.0
for i in range(image.nwav):
    image_ave.image[:,:,0] = image_ave.image[:,:,0] + image.image[:,:,i]*dv
    image_ave.imageJyppix[:,:,0] = image_ave.imageJyppix[:,:,0] + image.imageJyppix[:,:,i]*dv

if remove_tiny:
    for i in range(nx):
        for j in range(ny):
            if image_ave.image[i,j,0] < tiny:
                image_ave.image[i,j,0] = 0.0
            if image_ave.imageJyppix[i,j,0] < tiny:
                image_ave.imageJyppix[i,j,0] = 0.0




#######################################################
#####  Plotting

print " Plotting... ", "\n"

#######################################################
#####  Make some Plot 
#########################################
n_row = 5
n_col = 6
jump = 0
axisadd = -0


#######################################################
#####  Convolve image with gaussian
#########################################

print "Convolve the image wiht Guassian"
image_conv = image.imConv(fwhm, gaus_pa_ctclk, dpc)
image_ave_conv = image_ave.imConv(fwhm, gaus_pa_ctclk, dpc)



##########################################
# define azimuthal extract function
#  could be ellipse or circle
##########################################

def azimuthal_Jy_avg(image, gap_au, gap_width, e_h, e_k, incli, rotang, sizepix):
    # input: image, gap_au, gap_width, e_h, e_k, incli, rotang, sizepix, nx
    # output: fin_vis
    ##################################
    # screen out the point in the ring at [i_in, i_out]

    au = 1.496e13

    fin0_x = []
    fin0_y = []
    fin0_vis = []
    
    # parameters for the center ellipse
    gap_min = gap_au - gap_width*0.5
    gap_max = gap_au + gap_width*0.5
    if gap_min < 0.0:
        gap_min =0.0
    
    # assuming sizepix_x = sizepix_y
    e_a_grid_max = gap_max*au/sizepix  # long semi-axis
    e_a_grid_min = gap_min*au/sizepix  # long semi-axis
    
    inclination = math.radians(incli)
    e_b_grid_max = e_a_grid_max * cos(inclination)     # short semi-axis
    e_b_grid_min = e_a_grid_min * cos(inclination)     # short semi-axis
    rotang = math.radians(rotang)
    m = image.shape[0]
    # convert integer to float in order to make sure we find
    #    the center of the image.
    # image.out: 1) 20 grids
    #               python array: 0, 1, ..., 19
    #               center is at point 10
    #            2) 19 grids
    #               python array: 0, 1, ..., 18
    #               center is at point 9.5
    i_2_au = sizepix/au
    for ii in range(m):
        i=float(ii)
        for jj in range(m):
            j=float(jj)
            if (e_a_grid_min == 0.0) and (e_b_grid_min == 0.0):
                if ( ((i-e_h)*cos(rotang)+(j-e_k)*sin(rotang))**2/e_a_grid_max**2 + ((i-e_h)*sin(rotang)-( j-e_k)*cos(rotang))**2/e_b_grid_max**2 <= 1.0**2):
                    fin0_x.append(i)
                    fin0_y.append(j)
                    fin0_vis.append(image[ii,jj])
            else:
                if ( ((i-e_h)*cos(rotang)+(j-e_k)*sin(rotang))**2/e_a_grid_max**2 + ((i-e_h)*sin(rotang)-( j-e_k)*cos(rotang))**2/e_b_grid_max**2 <= 1.0**2) and ( ((i-e_h)*cos(rotang)+(j-e_k)*sin(rotang))**2/e_a_grid_min**2 + ((i-e_h)*sin(rotang)-(j-e_k)*cos(rotang))**2/e_b_grid_min**2 > 1.0**2) :
                    fin0_x.append(i)
                    fin0_y.append(j)
                    fin0_vis.append(image[ii,jj])
    fin_x = np.asarray(fin0_x)
    fin_y = np.asarray(fin0_y)
    fin_vis = np.asarray(fin0_vis)
    n_fin = fin_x.shape[0]
    if n_fin > 0: 
        np.savetxt('x_y_vis', np.transpose([fin_x,fin_y,fin_vis]))
    #total=np.sum(fin_vis), "\n"
    #avg = total/n_fin
    return fin_vis


##########################################
##########################################

#### define the ImagePlot nx*nx array to extract and plot azimuthal ring
ImagePlot = image_ave_conv.imageJyppix[:,:,n_image]


sizepix = image_ave_conv.sizepix_x
print "image_ave_conv.sizepix_x/au,image.sizepix_x/au"
print image_ave_conv.sizepix_x/au,image.sizepix_x/au
print nx, image_ave_conv.nx

nfloat = float(nx)
nhalfpix = nfloat/2
e_h = nhalfpix
e_k = nhalfpix


r_Jy = np.zeros([n_rad,2], dtype=float64)
for i in range(n_rad):
    gap_au = float(i)*d_rad+inneravoid
    r_Jy[i,0] = gap_au
    avg = azimuthal_Jy_avg(ImagePlot, gap_au, gap_width, e_h, e_k, incli, rotang, sizepix)
    n_num = avg.shape[0]
    if n_num == 0:
        r_Jy[i,1] = 0
        print "WARNNING: no points found around ", r_Jy[i,0], " AU!!"
    else:
        f_n_num = float(n_num)
        #print avg
        total = np.sum(avg)
        #print total
        r_Jy[i,1] = total/f_n_num
        print r_Jy[i,0], r_Jy[i,1]
np.savetxt('AU_fluxJy', np.transpose([r_Jy[:,0], r_Jy[:,1]]))
    


quit()
#######################################################
#####  Make some Plot 
#########################################

#######################################################
#########################################
#plt, axes = plt.subplots(nrows=1, ncols=3, figsize=(13, 6), dpi=80, sharex=True, sharey=True)
plt, axes = plt.subplots(n_row, n_col, figsize=(15, 12), dpi=80, sharex=True, sharey=True)


v_min   = 0.0
v_max   = image_conv.imageJyppix[:,:,:].max()*1000

for i in range(n_row):
    for j in range(n_col):
        if colorlog:
            im=axes[i,j].pcolormesh(im_x_au, im_y_au, image_conv.imageJyppix[:,:,n_col*i+j+jump]*1000.0, cmap='jet', norm=LogNorm(linthresh=0.03, linscale=0.03, vmin=v_min, vmax=v_max))
        else:
            im=axes[i,j].pcolormesh(im_x_au, im_y_au, image_conv.imageJyppix[:,:,n_col*i+j+jump]*1000, cmap='jet', vmin=v_min, vmax=v_max)
        axes[i,j].tick_params(axis='both', which='major', labelsize=8)
        #tick.label.set_fontsize(14) 
        axes[i,j].axis([im_x_au.min()-axisadd, im_x_au.max()+axisadd, im_y_au.min()-axisadd, im_y_au.max()+axisadd])
        axes[i,j].set_yticks( [-600, -400, -200, 0, 200, 400, 600] )
        axes[i,j].set_xticks( [-600, -400, -200, 0, 200, 400, 600] )
        axes[i,j].set_xticklabels(['-600', '-400', '-200', '0', '200', '400', '600'],rotation=90)
        if j==0:
            axes[i,j].set_ylabel('AU')
        axes[i,j].set_xlabel('AU')
        #for tick in axes[i,j].get_xticklabels():
        #    tick.set_rotation(45)

#plt.tight_layout()
plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9, hspace=0.0, wspace=0.0)
cax = plt.add_axes([0.82, 0.1, 0.02, 0.8])
plt.colorbar(im, cax=cax, label='m Janskey/beam')
#cax = plt.axes([0.85, 0.1, 0.075, 0.8])
#plt.colorbar(cax=cax)

#plt.subplots_adjust(hspace=0.0, wspace=0.0)
plt.savefig('gas_lines_radmc_convolved.png')

plt.clf()


