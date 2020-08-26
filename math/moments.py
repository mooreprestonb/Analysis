#!/usr/bin/python3

import scipy.stats
import numpy
import sys

if(len(sys.argv) != 2):
    print("Error need an data file to parse!")
    print(sys.argv)
    exit(1)
data = numpy.loadtxt(sys.argv[1])

#print(data.shape)
print("#Col average stdev variance kurtosis skew")
for i in range(data.shape[1]):
    d = data[:,i] # parse out i'th column
    print(i+1,numpy.mean(d),numpy.std(d),numpy.var(d),scipy.stats.skew(d),scipy.stats.kurtosis(d,fisher=True))
