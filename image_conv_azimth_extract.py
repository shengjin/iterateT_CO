#!/usr/bin/env python2

import numpy as np

import matplotlib as mpl
mpl.use('Agg')  # to avoid window

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from readimage import *
from natconst import *
from matplotlib import cm
import math
from math import sin, cos


remove_tiny = True
tiny  = 1e-30

# inclination of the disk
incli = 51.0
#incli = 46.72

# position angle
# use it to scale the center Radius VS T_B
#posang = 0.0
posang = 151.0
#posang = 138.02

# rotation angle used in the center ellipse to
rotang = posang
#rotang = 90-posang

# radius_wanted
n_rad = 23
d_rad = 10
inneravoid = 0
gap_width = 10.0 # AU

n_image = 0

# center of the ellipse in pixel

dpc = 140.0  # distance in pc
n_image = 0  # the n_th image in image.out

# LkCa 15
gaus_ma = 6.323864890469E-05    # maj axis of guassian in deg
gaus_mi = 4.604819748137E-05    # min
gaus_pa = 2.383443450928E+01    # PA of maj axis of Guass measure East from North

# HL Tau
#gaus_ma = 8.367259158856E-06                                                  
#gaus_mi = 5.335457002123E-06                                                  
#gaus_pa = -1.758292236328E+02 

gaus_ma_arcsec = gaus_ma*3600
gaus_mi_arcsec = gaus_mi*3600
gaus_pa_ctclk  = gaus_pa  # PA in imConv (counts from North counterclockwise)

fwhm = [gaus_ma_arcsec, gaus_mi_arcsec]

#debug = True
debug = False

AU = 1.496e13


#########################
#### print the gaussian
print "gaus_ma_arcsec, gaus_mi_arcsec"
print gaus_ma_arcsec, gaus_mi_arcsec
gaus_ma_rad = gaus_ma*math.pi/180.0
gaus_mi_rad = gaus_mi*math.pi/180.0
gaus_ma_AU = math.sin(gaus_ma_rad)*dpc*pc/AU
gaus_mi_AU = math.sin(gaus_mi_rad)*dpc*pc/AU
print "gaus_ma_AU, gaus_mi_AU"
print gaus_ma_AU, gaus_mi_AU

###################### ###################### 
########### Read image and scale the flux
###################### ###################### 

image = readImage()

nx = image.nx
ny = nx

if remove_tiny:
    for i in range(nx):
        for j in range(ny):
            if image.image[i,j,0] < tiny:
                image.image[i,j,0] = 0.0
            if image.imageJyppix[i,j,0] < tiny:
                image.imageJyppix[i,j,0] = 0.0



image_conv = image.imConv(fwhm, gaus_pa_ctclk, dpc)

ImagePlot = image_conv.imageJyppix[:,:,n_image]


##########################################
# set the grids in pixel
########################################## 
# NOTE: in this case, nx=ny=m=n
# assuming the same number of grids in x,y axes
nx = image.nx
NN = nx
pix_x = np.linspace(-NN/2,NN/2,nx, endpoint=False)
ny = image.ny
pix_y = np.linspace(-NN/2,NN/2,ny, endpoint=False)

# set the grids in AU
im_x_au = (pix_x)*image.sizepix_x/au
im_y_au = (pix_y)*image.sizepix_y/au



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
        gap_min = 0.0
    
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
    #print i_2_au,"i_2_au"
    #print e_h, e_k,"e_h,e_k"
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
        #np.savetxt('x_y_vis', np.transpose([fin_x,fin_y,fin_vis]))
        print n_fin 
    else:
        print fin_x,fin_y,n_fin
    #total=np.sum(fin_vis), "\n"
    #avg = total/n_fin
    return fin_vis

##########################################
##########################################

sizepix = image_conv.sizepix_x
print image_conv.sizepix_x/au,image.sizepix_x/au
print nx, image_conv.nx

nfloat = float(nx)
nhalfpix = nfloat/2
e_h = nhalfpix
e_k = nhalfpix


r_Jy = np.zeros([n_rad,2], dtype=float64)
for i in range(n_rad):
    gap_au = float(i)*d_rad+inneravoid
    print "gap_au", gap_au
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
        #print r_Jy[i,0], r_Jy[i,1]
np.savetxt('AU_fluxJy', np.transpose([r_Jy[:,0], r_Jy[:,1]]))
    
quit()
    

plt.xlabel('AU')
plt.ylabel("m Janskey / beam")
#plt.ylim(-180,180)
#plt.xlim(0,n_rad+inneravoid)
plt.plot(r_Jy[:,0], r_Jy[:,1]*1000)
#plt.pcolormesh(im_x_au, im_y_au, ImageJyppix_scaled*1000, vmin=1e-6, vmax=0.002,norm=LogNorm(),cmap='RdBu')
plt.savefig('azimuthalavg_fluxJy.png')
plt.clf()



#######################################################
#####  Make some Plot 
#########################################

plt.xlabel('AU')
plt.ylabel('AU')
plt.ylim(-100,100)
plt.xlim(-100,100)
plt.pcolormesh(im_x_au, im_y_au, ImagePlot*1000, vmin=0, vmax=18) #, norm=LogNorm())
#plt.pcolormesh(im_x_au, im_y_au, ImagePlot*1000,cmap='RdBu')
#plt.pcolormesh(im_x_au, im_y_au, ImageJyppix_scaled*1000, vmin=1e-6, vmax=0.002,norm=LogNorm(),cmap='RdBu')
cbar1=plt.colorbar()
cbar1.set_label("m Janskey / beam")
plt.title("Flux density")
plt.savefig('flux_density.png')
plt.clf()






