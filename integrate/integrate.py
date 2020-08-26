#!/usr/bin/python3

import numpy
from scipy.integrate import simps

filename = "0.5mgr6all.dat"
data = numpy.loadtxt(filename)

#print(data)

xdata = data[:,0] # get first column
ydata = data[:,1] # get second column

IS = simps(ydata,xdata)
IT = numpy.trapz(ydata,xdata)
print(IT,IS)

sum = 0
for i in range(len(xdata)-1):
    sum += ((ydata[i]+ydata[i+1])/2)*(xdata[i+1]-xdata[i])

print(sum)
