#!/usr/bin/python3

import sys
import numpy
import argparse
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

#------------------------------------------------------------
# read command line arg using argparser
parser = argparse.ArgumentParser(description="Read in density file of slab geometry and get gas and liquid density")
parser.add_argument('input', help='input density file name')
parser.add_argument('-s', default=3.0, help="stdev to away from phase transition",dest='stdev',type=float)
parser.add_argument('-p', default=False, action="store_true", dest="p",
		    help="Show graphs (Default off)")

# read in arguments, now translate to variables
args = parser.parse_args()
data = numpy.loadtxt(args.input)

h = data[1][0] - data[0][0]
x = data[:,0]
dens = data[:,1]
densp = (dens[2:]-dens[:-2])/(2*h)

avg = numpy.mean(densp)
std = numpy.std(densp)

def model2g(x,a1,c1,w1,a2,c2,w2): # two gaussians
    return a1*numpy.exp(-0.5*(x-c1)**2/w1**2)+a2*numpy.exp(-0.5*(x-c2)**2/w2**2)

ximin = numpy.argmin(densp)
ymin = densp[ximin]
ximax = numpy.argmax(densp)
ymax = densp[ximax]

ivals = [ymax,x[ximax],1.0,ymin,x[ximin],1.0] # inital guess
yi = model2g(x[1:-1],*ivals)  # splat or unpacking argument list
fvals, covar = curve_fit(model2g,x[1:-1],densp,p0=ivals)
yf = model2g(x[1:-1],*fvals)

stdint = args.stdev

print("fit",fvals,"h",h,"s",stdint)
intface = numpy.zeros(4,dtype=int) # Take interface stdint stdev.
intface[0] = max(int((fvals[1]-x[0])/h-stdint*fvals[2]/h),0)
intface[1] = int((fvals[1]-x[0])/h+stdint*fvals[2]/h)+1
intface[2] = int((fvals[4]-x[0])/h-stdint*fvals[5]/h)
intface[3] = int((fvals[4]-x[0])/h+stdint*fvals[5]/h)+1
print("Interface:",intface)

if(args.p):
    print(fvals)
    plt.plot(x, dens, 'bo-', label= 'Density')
    plt.plot(x[1:-1], densp, 'ko', label="d Density/dx")
    plt.plot(x[1:-1], yi, 'g-', label='init')
    plt.plot(x[1:-1], yf, 'r-', label='fit')
    plt.plot(x[1:intface[0]], dens[1:intface[0]], 'c-', label= 'Gas1')
    plt.plot(x[intface[1]:intface[2]], dens[intface[1]:intface[2]], 'c-', label= 'Liquid')
    plt.plot(x[intface[3]:-2], dens[intface[3]:-2], 'c-', label= 'Gas2')
    plt.legend(loc='best')
    plt.show()

# don't use first of last points
avg1 = numpy.mean(dens[1:intface[0]])
std1 = numpy.std(dens[1:intface[0]])
avg2 = numpy.mean(dens[intface[1]:intface[2]])
std2 = numpy.std(dens[intface[1]:intface[2]])
avg3 = numpy.mean(dens[intface[3]:-2])
std3 = numpy.std(dens[intface[3]:-2])

print(dens[1:intface[0]])
print(dens[intface[1]:intface[2]])
print(dens[intface[3]:-2])
print(sys.argv[1])
print("Phase 1",avg1,std1)
print("Phase 2",avg2,std2)
print("Phase 3",avg3,std3)


