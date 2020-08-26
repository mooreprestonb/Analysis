#!/usr/bin/python3

import sys
import numpy
import math
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def histogram(r,hist,hmin,dx,nbins):
    pt = (r-hmin)/dx
    i = int(pt)
    #        if (i < (nbins-1)) and (i >= 0):
    #                hist[i] += 1
    fr = pt-i # fraction of bin width
    print(r,i,fr)
    hist[i] += 1.0-fr
    hist[i+1] += fr

def model2g(x,a1,c1,w1,a2,c2,w2):
    return a1*numpy.exp(-(x-c1)**2/w1)+a2*numpy.exp(-(x-c2)**2/w2)

print (sys.argv)
data = numpy.loadtxt(sys.argv[1])
dens = data[:,1]

nbins = 200
hist = numpy.zeros(nbins)
xhist = numpy.zeros(nbins)
rmin = numpy.amin(dens)
rmax = numpy.amax(dens)
dx = (rmax-rmin)/float(nbins-3)
rmin -= dx
print(rmin,rmax,dx,nbins)

for i in range(len(dens)):
    histogram(dens[i],hist,rmin,dx,nbins)
    
for i in range(len(xhist)):
    xhist[i] = i*dx+rmin

numpy.savetxt("hist.dat",numpy.column_stack((xhist,hist)))

ivals = [10,rmin+dx,.001,10,rmax-dx,.001] # inital guess
yi = model2g(xhist,*ivals)  # splat or unpacking argument list

fvals, covar = curve_fit(model2g,xhist,hist,p0=ivals)
print(fvals)
yf = model2g(xhist,*fvals)

plt.plot(xhist, hist, 'bo')
plt.plot(xhist, yi, 'g-', label='init')
plt.plot(xhist, yf, 'r-', label='fit')
plt.legend(loc='best')
plt.show()


