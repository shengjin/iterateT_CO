#!/usr/bin/env python2
import numpy as np

Sigma_0 = [0.6e-3, 3.0e-3, 6.0e-3]
rc  = [12500]
gamma = [0.5, 0.8, 1.1, 1.4]
rc_atan  = [5, 15, 25, 35]
gamma_atan = [1, 2, 3, 4]

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

#result_dir = "0001"
#onerun(result_dir)
    
