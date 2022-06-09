#!/usr/bin/python3

import sys              # system arguments
import argparse         # command line arguments
import numpy            # numerical routines
import re               # regular expressions

#------------------------------------------------------------
# read command line arg using argparser
parser = argparse.ArgumentParser(description="Read in number density file and create density in g/cc along z")
parser.add_argument('input', help='input density file name')
parser.add_argument('output', help='output density file name')

# read in arguments, now translate to variables
args = parser.parse_args()
#nbins = args.nbins
infile = args.input 
outfile = args.output

print(sys.argv)
data = numpy.loadtxt(infile) # import data

#read header 
f = open(infile,'r')
h1 = f.readline()
h2 = f.readline()
f.close()

#print(data)
#print(infile,outfile)
#print(h1)
#print(h2)

# parse header
words = h2.split()
#print (words)
dz = float(words[3])
box = [float(words[5]),float(words[6]),float(words[7])]
nconf = int(words[-1])

ntypes = int(words[8][6:])
artypes = words[9:-2]
print (nconf,dz,box)

if(ntypes != len(artypes)):
    print("Error! # artypes != ntypes")
    exit(1)
    
types = {}
for ar in artypes:
    m = re.match(r"(\d+):(\d+)",ar)
    # print(m.group(1),m.group(2))
    types[int(m.group(1))] = int(m.group(2))
print (types)

if (data.shape[1] == 3):
    mass = numpy.array((16.00,1.00)) # masses of each type (just water)
elif (data.shape[1] == 5):
    mass = numpy.array((16.00,1.00, 22.99,35.45)) # masses of each type H20 + NaCl
else:
    print("# of data types doesn't match mass!",data.shape)
    exit(1)

print(mass)
# .602214 = 6.02214x10^23/(1*x10^8)^3  atoms/bin to mol/cc
fact = 1./(.602214*box[0]*box[1]*dz) # convert from n/bin to g/cc

ws = numpy.sum(data[:,1:]*mass*fact,axis=1) # get weighted sum
numpy.savetxt(outfile,numpy.column_stack([data[:,0],ws]),fmt='%g') # save output
