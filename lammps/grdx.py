#!/usr/bin/python3
# calculate g(r) grom g(dx,dy,dz) 

import sys
import argparse
import numpy
import math
import os.path
from scipy import stats
from scipy import signal
import matplotlib.pyplot as plt

#-----------------------------------------------------------
def readdxfile(ifile,np,bins):
    f = open(ifile,"r") #open file to read
    line = f.readline()
    while(line[0]=="#"):  # loop over comments
        print(line)
        line = f.readline()

    words = line.split()
    # "object 1 class gridpositions counts "+str(np[0])+" "+str(np[1])+" "+str(np[2])+"\n"
    np[0] = int(words[5])
    np[1] = int(words[6])
    np[2] = int(words[7])
    
    #"origin " + str(bins[0]) + " " + str(bins[0]) + " " + str(bins[0]) +"\n"
    words = f.readline().split()
    # print(words) 
    bins[0] = float(words[1])

    # "delta "+str(bins[2])+" 0 0\ndelta 0 "+str(bins[2])+" 0\ndelta 0 0 "+str(bins[2])+"\n"
    words = f.readline().split()
    # print(words) 
    bins[2] = float(words[1])
    bins[1] = (np[0]-1)*bins[2]+bins[0]

    f.readline() # next deltas
    f.readline() # next delta
    
    #hdr = "object 2 class gridconnections counts "+str(np[0])+" "+str(np[1])+" "+str(np[2])+"\n"
    line = f.readline()
    # hdr = "object 3 class array type double rank 0 items "+str(np[0]*np[1]*np[2])+" data follows\n"
    line = f.readline()
    print(line)
    grdx = numpy.zeros((np[0],np[1],np[2]))
    print(grdx.shape)
    #read data
    for i in range(np[0]):
        for j in range(np[1]):
            for k in range(np[2]): 
                words = f.readline().split()
                # print(words)
                grdx[i][j][k] = float(words[0])

    f.close
    print("Parsed dxfile",ifile)
    return(grdx)
    
#-------------------------------------------------------------------------

# read command line arg using argparser
parser = argparse.ArgumentParser(description="Read in OpenDX file and create g(r)")
parser.add_argument('input', help='input dx file')
#parser.add_argument('output', help='output file name')

# read in arguments, now translate to variables
args = parser.parse_args()
print (args)
infile = args.input # 'test.lammpstrj'
#outfile = args.output

#Grids
np = numpy.zeros(3,dtype=int)
bins = numpy.zeros(3)

grid = readdxfile(infile,np,bins)

print(np,bins,grid.shape)

pt = numpy.zeros(np[0]*np[1]*np[2])
wt = numpy.zeros(np[0]*np[1]*np[2])

n=0
npt = 0
fact = bins[2]**3/(4.*math.pi*bins[2])
for i in range(np[0]):
    x = i*bins[2]+bins[0]
    x2 = x*x
    for j in range(np[1]):
        y = j*bins[2]+bins[0]
        y2 = y*y
        for k in range(np[2]):
            z = k*bins[2]+bins[0]
            #print(x,y,z,grid[i][j][k])
            r = math.sqrt(x2+y2+z*z)
            # print(r,grid[i][j][k])
            if(r<=bins[1]):
                npt += 1
            pt[n] = r
            wt[n] = grid[i][j][k] 
            n += 1

print(n,npt)

gr,hbins = numpy.histogram(pt,bins=np[0]*2,range=(0.,bins[1]),weights=wt)
dh = hbins[1]-hbins[0]
xdata = hbins[:-1]-dh
for i in range(len(gr)):
    gr[i] *= bins[2]**3 * 3./(4.0*math.pi*(hbins[i+1]**3-hbins[i]**3)) # normalize

plt.plot(xdata,gr)

#kernel = stats.gaussian_kde(pt)
#plt.plot(xdata,kernel(xdata))


gr2 = numpy.zeros(len(gr))
idx = numpy.where(pt<=hbins[-1])
ptn = pt[idx]
wt = wt[idx]
ptn = ptn/dh # get bin numbers (min dist = 0)
nr = ptn.astype(int) # get integer bin numbers
ptn -= nr # get fraction of bin

hlen = len(gr)-1
for i in range(len(nr)):
    ii = nr[i]
    #gr2[ii] +=  wt[i]
    if(ii<hlen):
        gr2[ii] += (1.-ptn[i])*wt[i]
        gr2[ii+1] +=  ptn[i]*wt[i]
    else:
        gr2[ii] += (1-ptn[i])*wt[i]
        

for i in range(len(gr)): # normalize
    gr2[i] *= bins[2]**3 * 3./(4.0*math.pi*(hbins[i+1]**3-hbins[i]**3))  

plt.plot(xdata,gr2)

y = signal.savgol_filter(gr2,5,2)  # smooth data
plt.plot(xdata,y)

plt.show()
