#!/usr/bin/python3
# parse lammps log files

import sys
import numpy
from scipy.stats import skew,kurtosis

fi = open(sys.argv[1],"r")
lines = fi.readlines()
fi.close()

ls = 0
ns = 0
for l in lines:
    words = l.split()
    ls += 1
    if(len(words)>0 and words[0]=="Step"):
        if(ns!=0):
            print("Warning more then one \"Step\" found (using latest)")
        ns = ls
        #print("found step")
        #break

if(ns==0):
    print("Didn't find \"Step\" in lammps log file... exiting")
    exit(1)
    
header = lines[ns-1].split()
nc = len(header)
n = len(lines[ns:])
data = numpy.empty((nc,n),dtype=float) #maybe more then we need :-)
nd = 0
for l in range(n):
    words = lines[ns+l].split()
    if(words[0].isdigit()==False):
        break
    for i in range(nc):
        data[i][l] = float(words[i])
    nd += 1
#print(data[0][:nd])
print("#",sys.argv,len(header),"columns of data, with",nd,"lines of data")
print("# column header average stdev skew kurtosis min max range [slop intercept]")
for i in range(nc):
    da = data[i][:nd] # only data that was parsed
    dmin = numpy.min(da)
    dmax = numpy.max(da)
    model = numpy.polyfit(range(nd),da,1)
    print(i,header[i],numpy.average(da),numpy.std(da),skew(da),kurtosis(da),dmin,dmax,dmax-dmin,model)
