#!/usr/bin/env python2
import numpy as np

gamma = [0.5, 1, 1.5, 2]
rc  = [12500]
rc_atan  = [15, 25, 35, 45]
gamma_atan = [1, 2, 3, 4]
Sigma_0 = [7.5e-4, 3.0e-3, 7.5e-3, 2.0e-2]

color= [5, 3, 7, 9]
# yellow, blue, red, purple


wfile = open("plot.txt", 'w')
h = 0
for i in range(len(Sigma_0)):
    for j in range(len(rc)):
	for k in range(len(gamma)):
	    for m in range(len(rc_atan)):
		for n in range(len(gamma_atan)):
                    wfile.write('%s%04d%s %01d %s\n'%("'./", h, "/Au_fluxJy_12CO' u 1:($2*1000) w l lt ", color[k], ", \\"))
		    h = h + 1
wfile.close()
