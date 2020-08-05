#!/usr/bin/python3
# kernel estimation using guassian kernal of x-y data

from scipy import stats
import numpy
import sys

if(len(sys.argv)<2):
    print("Usage",sys.argv[0]," <infile> <outfile> [column] [npoints]")
    exit(1)

npt = 'auto'
column = 2
if(len(sys.argv)>3):
    column = int(sys.argv[3])

if(len(sys.argv)>4):
    npt = int(sys.argv[4])
    
vals = numpy.loadtxt(sys.argv[1])
data = vals[:,column-1]
#print(data)
kern = stats.gaussian_kde(data)

amin = numpy.min(data)
amax = numpy.max(data)

hist, ebins = numpy.histogram(data,bins=npt,density=True)
#print(hist,ebins)
xval = ebins[:-1]+(ebins[1]-ebins[0])/2.
print (len(xval),ebins[0],ebins[-1])
#print(xval)
yval = kern(xval)

hdr = ""
for i in range(2):
    hdr += " " + sys.argv[i]
hdr += " " + str(column) + " " + str(npt)
numpy.savetxt(sys.argv[2],numpy.c_[xval,yval,hist],fmt='%g',header=hdr)

