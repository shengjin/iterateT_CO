#!/usr/bin/env python2
import numpy as np

gamma = [0.5, 1, 1.5, 2]
rc  = [12500]
rc_atan  = [15, 25, 35, 45]
gamma_atan = [1, 2, 3, 4]
Sigma_0 = [7.5e-4, 3.0e-3, 7.5e-3, 2.0e-2]

#nrun = len(gamma)*len(rc)*len(gamma_atan)*len(rc_atan)*len(Sigma_0)

wfile = open("run_grid.txt", 'w')
h = 0
for i in range(len(Sigma_0)):
    for j in range(len(rc)):
	for k in range(len(gamma)):
	    for m in range(len(rc_atan)):
		for n in range(len(gamma_atan)):
		    wfile.write('%04d %.2e %.2e %.2e %.2e %.2e \n'%(h, Sigma_0[i], rc[j], gamma[k], rc_atan[m], gamma_atan[n]))
		    #wfile.write('%4i %.4f %.1f %.1f %.1f %.1f \n'%(h, Sigma_0[i], rc[j], gamma[k], rc_atan[m], gamma_atan[n]))
		    h = h + 1
wfile.close()

result_dir = "0001"
#onerun(result_dir)
    
