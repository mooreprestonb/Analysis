#!/usr/bin/python3

# parse lammps log files

import sys
import numpy

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

header = lines[ls-1].split()
nc = len(header)
n = len(lines[ls:])
print(len(header),"Columns of data")
print(n,"lines of data")
data = numpy.empty((nc,n),dtype=float)
for l in range(n):
    words = lines[ls+l].split()
    if(words[0].isdigit()==False):
        break
    for i in range(nc):
        data[i][l] = float(words[i])

for i in range(nc):
    print(i,header[i],numpy.average(data[i]),numpy.std(data[i]))

