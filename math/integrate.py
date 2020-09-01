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


maxx = 15
if(len(sys.argv)!=2):
        print(sys.argv[0]," <filename>")
        print("Performs numerical integration (using simpsons rule) on xy data")
        print("with an x cutoff of ",maxx)
        exit(1)
        
filename = sys.argv[1]
#filename = "0.5mgr6all.dat"
data = numpy.loadtxt(filename)

#print(data)

xdata = data[:,0] # get first column
ydata = data[:,1] # get second column

# find maxx index and linear extrapolate
idx = numpy.searchsorted(xdata,maxx)
#print(idx,xdata[idx-1:idx+1],ydata[idx-1:idx+1])
xdata1 = xdata[:idx+1]
ydata1 = ydata[:idx+1]
dx = xdata1[-1]-xdata1[-2]
dy = ydata1[-1]-ydata1[-2]
#print (xdata1[-2:],ydata1[-2:],dy,dx,(dy/dx),(maxx-xdata1[-2]))
# set last point
xdata1[-1] = maxx
ydata1[-1] = ydata1[-2] + (maxx-xdata1[-2])*dy/dx
#print (xdata1[-2:],ydata1[-2:])

#perform integration
IS = simps(ydata1,xdata1) # simpsons rule
print(IS)
#traptest(ydata1,xdata1) # test using trap rule?

