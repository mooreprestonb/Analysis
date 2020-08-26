#!/usr/bin/python3
# curve fitting (guassian in this example)
# for polynomials could use ... model = numpy.polyfit(x,y,degree,...)

from scipy.optimize import curve_fit
import numpy
import sys

def func(x,a,b,c):
    return a* numpy.exp((-(x-b)**2)/(2.*c*c))

if(len(sys.argv)!=3):
    print("Usage:",sys.argv[0]," <infile> <outfile>")
    exit(1)

data = numpy.loadtxt(sys.argv[1])
xdata = data[:,0]
ydata = data[:,1]
print(xdata[-1],xdata[0])
parms =[numpy.max(ydata),(xdata[-1]-xdata[0])/2.+xdata[0],((xdata[-1]-xdata[0])/4.)]
print("Initial guess at parms:",parms)
ytest = func(xdata,parms[0],parms[1],parms[2])
#print(xdata,ydata,ytest)

popt, pcov = curve_fit(func,data[:,0],data[:,1],p0=parms)
print("parms:",popt)

fdata = func(xdata,popt[0],popt[1],popt[2])
hdr = ""
for i in range(len(sys.argv)):
    hdr += " " + sys.argv[i]
hdr = "Parms:"
for i in range(len(popt)):
    hdr += " " + str(popt[i])

numpy.savetxt(sys.argv[2],numpy.c_[xdata,ydata,fdata],header=hdr,fmt='%g')
