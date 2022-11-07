#!/usr/bin/python3
# calculate g(r) grom g(dx,dy,dz) 

import sys
import argparse
import numpy
import math
import os.path
from scipy import stats
import matplotlib.pyplot as plt

#-----------------------------------------------------------
def readdxfile(ifile,np,bins):
    f = open(ifile,"r") #open file to read
    line = f.readline()
    while(line[0]=="#"):  # loop over comments
        print(line)
        line = f.readline()
    print(line)
    words = line.split()
    # "object 1 class gridpositions counts "+str(np[0])+" "+str(np[1])+" "+str(np[2])+"\n"
    np[0] = int(words[5])
    np[1] = int(words[6])
    np[2] = int(words[7])
    
    #"origin " + str(bins[0]) + " " + str(bins[0]) + " " + str(bins[0]) +"\n"
    words = f.readline().split()
    # print(words) 
    bins[0] = float(words[1])

    # "delta "+str(bins[2])+" 0 0\ndelta 0 "+str(bins[2])+" 0\ndelta 0 0 "+str(bins[2])+"\n"
    words = f.readline().split()
    # print(words) 
    bins[2] = float(words[1])
    bins[1] = (np[0]-1)*bins[2]+bins[0]
    print(np,bins)
    exit()
    f.readline() # next deltas
    f.readline() # next delta
    
    #hdr = "object 2 class gridconnections counts "+str(np[0])+" "+str(np[1])+" "+str(np[2])+"\n"
    line = f.readline()
    # hdr = "object 3 class array type double rank 0 items "+str(np[0]*np[1]*np[2])+" data follows\n"
    line = f.readline()
    print(line)
    grdx = numpy.zeros((np[0],np[1],np[2]))
    print(grdx.shape)
    #read data
    for i in range(np[0]):
        for j in range(np[1]):
            for k in range(np[2]): 
                words = f.readline().split()
                # print(words)
                grdx[i][j][k] = float(words[0])

    f.close
    print("Parsed dxfile",ifile)
    return(grdx)
    
#-------------------------------------------------------------------------

# read command line arg using argparser
parser = argparse.ArgumentParser(description="Read in OpenDX file and create g(r)")
parser.add_argument('input', help='input dx file')
parser.add_argument('output', help='output file name')

# read in arguments, now translate to variables
args = parser.parse_args()
print (args)
infile = args.input # 'test.lammpstrj'
outfile = args.output

#Grids
np = numpy.zeros(3,dtype=int)
bins = numpy.zeros(3)

grid = readdxfile(infile,np,bins)

print(np,bins,grid.shape)

