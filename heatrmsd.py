#!/usr/bin/python3
import numpy
import sys

print("hello")

#readin data
filein = open("avg-receptor-5vai-heatmap-100-frames.hm",'r')
lines = filein.readlines()
filein.close()

fileout = open("betavals-5vai.dat",'w')
fileout.write("# "+str(sys.argv)+"\n")
data = numpy.empty((394,100)) # should get from first line
#print(lines[:10])
for i in range(394): # loop over heatmap lines got get numbers
    vals = lines[10+i].split(";") # get numbers
    vals[0] = vals[0].split(":")[2] # fix fisrt value (remove residue, bracket 2 picks up at 3rd value)
    vals.pop() # remove end of line characture
    val = numpy.array([float(f) for f in vals])
    data[i] = val
    fileout.write(str(i+1)+" "+str(numpy.average(val))+" "+str(numpy.std(val))+"\n")
#print(data)

fileout.close()
print("goodbye")
