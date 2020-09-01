#!/usr/bin/python3

import sys
import numpy
from scipy.integrate import simps

def traptest(ydat,xdat):
    IT = numpy.trapz(ydat,xdat) # trapazod rule
    print(IT)

    # trapazod rule by hand (slow)
    sum = 0
    for i in range(len(xdata1)-1):
        sum += ((ydata1[i]+ydata1[i+1])/2)*(xdata1[i+1]-xdata1[i])
        
    print(IT,sum,abs(IT-sum))


filename = sys.argv[1]
#filename = "0.5mgr6all.dat"
maxx = 15
data = numpy.loadtxt(filename)

#print(data)

xdata = data[:,0] # get first column
ydata = data[:,1] # get second column

# find maxx index and linear extrapolate
idx = numpy.searchsorted(xdata,15)
#print(idx,xdata[idx-1:idx+1],ydata[idx-1:idx+1])

xdata1 = xdata[:idx+1]
ydata1 = ydata[:idx+1]
dx = xdata1[-1]-xdata1[-2]
dy = ydata1[-1]-ydata1[-2]
#print (xdata1[-2:],ydata1[-2:],dy,dx,(dy/dx),(maxx-xdata1[-2]))
xdata1[-1] = maxx
ydata1[-1] = ydata1[-2] + (maxx-xdata1[-2])*dy/dx
#print (xdata1[-2:],ydata1[-2:])

IS = simps(ydata1,xdata1) # simpsons rule
print(IS)
#traptest(ydata1,xdata1)

