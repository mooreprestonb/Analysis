#!/usr/bin/python3

import sys
import numpy
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

#print (sys.argv)
data = numpy.loadtxt(sys.argv[1])
h = data[1][0] - data[0][0]
x = data[:,0]
dens = data[:,1]
densp = (dens[2:]-dens[:-2])/(2*h)

avg = numpy.mean(densp)
std = numpy.std(densp)

def model2g(x,a1,c1,w1,a2,c2,w2):
    return a1*numpy.exp(-(x-c1)**2/w1)+a2*numpy.exp(-(x-c2)**2/w2)

ximin = numpy.argmin(densp)
ymin = densp[ximin]
ximax = numpy.argmax(densp)
ymax = densp[ximax]

ivals = [ymax,x[ximax],1.0,ymin,x[ximin],1.0] # inital guess
yi = model2g(x[1:-1],*ivals)  # splat or unpacking argument list
fvals, covar = curve_fit(model2g,x[1:-1],densp,p0=ivals)
yf = model2g(x[1:-1],*fvals)

print(fvals)
plt.plot(x, dens, 'bo', label= 'Density')
plt.plot(x[1:-1], densp, 'ko', label="d Density/dx")
plt.plot(x[1:-1], yi, 'g-', label='init')
plt.plot(x[1:-1], yf, 'r-', label='fit')
plt.legend(loc='best')
plt.show()

#print(fvals,h,x[0])
intface = numpy.zeros(4,dtype=int) # liquid side is 2*std.
intface[0] = int((fvals[1]-x[0])/h-fvals[2]/h)
intface[1] = int((fvals[1]-x[0])/h+2*fvals[2]/h)
intface[2] = int((fvals[4]-x[0])/h-2*fvals[5]/h)
intface[3] = int((fvals[4]-x[0])/h+fvals[5]/h)
print(intface)

avg1 = numpy.mean(dens[1:intface[0]])
std1 = numpy.std(dens[1:intface[0]])
avg2 = numpy.mean(dens[intface[1]:intface[2]])
std2 = numpy.std(dens[intface[1]:intface[2]])
avg3 = numpy.mean(dens[intface[3]:-2])
std3 = numpy.std(dens[intface[3]:-2])

print("phase 1 density: ",dens[1],dens[interface1[0]-2])
print("phase 2 density: ",dens[interface1[-1]+3],dens[interface2[0]-1])
print("phase 3 density: ",dens[interface2[-1]+5],dens[-2])
print(sys.argv[1],"AvgGas:",avg1,avg3," StdGas:",std1,std3,"AvgLiquid:",avg2,std2)

numpy.savetxt("hist.dat",numpy.column_stack((xhist,hist)))

