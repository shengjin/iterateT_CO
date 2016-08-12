#########################################
# read amr_grid to get the grid info

grid = linecache.getline('amr_grid.inp', 6)
print grid
ngrids = [int(s) for s in grid.split() if s.isdigit()]
nrad  = ngrids[0]
nthet = ngrids[1]
nphi  = ngrids[2]
print "nrad", nrad
print "nthet", nthet
print "nphi", nphi


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


######### ######### ######### 
###### the function that store one temperature date (nmodify layers)

def storetemperature(nfile):
    global nrad, nphi, nmodify, temp_nmodifytotal, temp_3d
    global temp_crit
    debug = True

    if debug:
	print "nfile", nfile

    for i in range(nmodify):
        for j in range(nphi):
            for k in range(nrad):
		if (temp_3d[j,nthet-1-i,k] > temp_crit):
		    temp_nmodifytotal[i,j,k,0] = temp_nmodifytotal[i,j,k,0] + temp_3d[j,nthet-1-i,k]
		    temp_nmodifytotal[i,j,k,1] = temp_nmodifytotal[i,j,k,1] + 1.0

    return temp_nmodifytotal


######## End Function Declaration
###################



################### ################### 
################### the main program
################### ################### 

######  Declare the array used to hold the temperature data 
temp_3d = np.zeros((nphi,nthet,nrad), dtype=float)
# 0: T Sum;   1: n_number 
temp_nmodifytotal = np.zeros((nmodify,nphi,nrad,2), dtype=float)

# remove the cell if T < temp_crit
#temp_crit = -10.01
temp_crit = 0.01

###################



#### Read all the temperature data for the nmodify layers

print "start reading mult files"

for i in range(ntotal):
    nfile = i+1 
    ext = nfile 
    ext = "%03d" % ext
    filename = "%s%s%s" % ("./", ext, "/dust_temperature.bdat")
    if os.path.isfile(filename):
        print "i, nfile", i, nfile
	print "Reading:", filename
	readtemperature(filename)
	storetemperature(nfile)

print "sum up mult files"

for i in range(nmodify):
## Note: we use the upper nthet-nmodify layers of the last temp_3d read
    for j in range(nphi):
        for k in range(nrad):
	    if ( temp_nmodifytotal[i,j,k,1] >= 1.0):
	       	temp_3d[j,nthet-1-i,k] = temp_nmodifytotal[i,j,k,0] / temp_nmodifytotal[i,j,k,1]
	    else:
	       	temp_3d[j,nthet-1-i,k] = 0.0


print "Modify the last %s layers of dust_temperature data" % nmodify



########### write the new temp_3d to dust_temperature.bdat
 
print "Write the new temperature data"


temp_3d = np.array(temp_3d).reshape(nphi*nthet*nrad)
#np.savetxt('temp.dat', temp_3d) for ascii file
print "tranformed to 1 row array for write"

with open('dust_temperature.bdat', 'wb') as f:
    f.write(struct.pack(4*'q',1,4,nphi*nthet*nrad,1))
    for num1 in temp_3d:
        f.write(struct.pack('f', num1))

