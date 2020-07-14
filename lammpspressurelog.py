#!/usr/bin/python3

import sys
import numpy
import math
import time

fi = open(sys.argv[1],"r")
lines = fi.readlines()
fi.close()

ls = 0
for l in lines:
    w1 = l.split()[0]
    ls += 1
    if(w1=="Step"):
        #print("found step")
        break

n = len(lines[ls:])
print(n," lines of data")
press = numpy.empty((6,n),dtype=float)
for l in range(n):
    words = lines[ls+l].split()
    for i in range(6):
        press[i][l] = float(words[10+i])

for i in range(6):
    print(i,numpy.average(press[i]),numpy.std(press[i]))

