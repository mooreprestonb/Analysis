#!/usr/bin/python3
# parse lammps log files

import sys
import numpy
import matplotlib.pyplot as plt
from scipy.stats import skew,kurtosis
from scipy import stats

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

ndata = data[:,:nd] # slice out only data that we parsed

print("#",sys.argv,len(header),"columns of data, with",nd,"lines of data")
print("# column header   average         stdev        skew    kurtosis       min       max       range   [slop intercept]")
for i in range(nc):
    da = ndata[i] # only data that was parsed
    dmin = numpy.min(da)
    dmax = numpy.max(da)
    model = numpy.polyfit(range(nd),da,1)
    print(f'{i:3} {header[i]:8} {numpy.average(da):13.4f} {numpy.std(da):13.4f} {skew(da):10.4f} {kurtosis(da):10.4f} {dmin:10.4g} {dmax:10.4g} {dmax-dmin:10.4g} [{model[0]:10.4g} {model[1]:10.4g}]')
    #print(i,header[i],numpy.average(da),numpy.std(da),skew(da),kurtosis(da),dmin,dmax,dmax-dmin,model)

    fn = numpy.arange(nd)
    plt.figure()
    plt.subplot(121)
    plt.plot(fn,ndata[i])
    plt.ylabel(header[i])
    plt.xlabel("Frames")
    plt.subplot(122)
    hst, bins = numpy.histogram(ndata[i],bins="auto",density=True)
    xdata = bins[:-1]+(bins[1]-bins[0])/2.0
    plt.plot(xdata,hst)
    kernel = stats.gaussian_kde(ndata[i])
    plt.plot(xdata,kernel(xdata))
    plt.xlabel(header[i])
    plt.show()
