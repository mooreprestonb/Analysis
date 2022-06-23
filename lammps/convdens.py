#!/usr/bin/python3

import sys              # system arguments
import argparse         # command line arguments
import numpy            # numerical routines
import re               # regular expressions

#------------------------------------------------------------
# read command line arg using argparser
parser = argparse.ArgumentParser(description="Read in number density file and create density in g/cc along z")
parser.add_argument('input', help='input density file name')
parser.add_argument('mass', help='mass file')
parser.add_argument('output', help='output density file name')
parser.add_argument('-u','--units',help='units (default = real)',choices=['real','lj'],default='real',nargs='?')

# read in arguments, now translate to variables
args = parser.parse_args()
#nbins = args.nbins
infile = args.input 
massfile = args.mass
outfile = args.output

print(sys.argv)
data = numpy.loadtxt(infile) # import data

#read header 
f = open(infile,'r')
h1 = f.readline()
h2 = f.readline()
f.close()

# parse header
words = h2.split()
#print (words)
dz = float(words[3])
box = [float(words[5]),float(words[6]),float(words[7])]
nconf = int(words[-1])

ntypes = int(words[8][6:])
artypes = words[9:-2]
print("Nconf=",nconf,"dz = ",dz,"Box = ",box)

if(ntypes != len(artypes)):
    print("Error! # artypes != ntypes")
    exit(1)
    
types = {}
for ar in artypes:
    m = re.match(r"(\d+):(\d+)",ar)
    # print(m.group(1),m.group(2))
    types[int(m.group(1))] = int(m.group(2))
ntypes = len(types)
print("Types {type:#}:",ntypes,types)

masses = numpy.loadtxt(massfile) # import mass data
if(ntypes==1): # only 1 mass!
    print(ntypes,masses.size,masses)
    if(masses.size!=2):
        print("Only 1 mass type and more than 1 type listed?")
        print(masses)
        exit(1)
    else:
        mass = numpy.array([masses[1]])
else:
    if(masses.shape[0] != ntypes):
       print("types and number of masses don't match!",ntypes,masses.shape[0])
       print(masses)
       exit(1)
    else:
       mass = numpy.zeros(ntypes)
       for i in range(ntypes):
           if(int(masses[i][0]) != i+1):
               print("Type out of order?",i,masses[i][0])
               print(masses)
               exit(1)
           mass[i] = masses[i][1]
       
print("Mass = ",mass)
if(args.units=='lj'):
    print("Using 'lj' units")
    fact = 1.0
else: # convert from n/bin to g/cc
    # .602214 = 6.02214x10^23/(1*x10^8)^3  atoms/bin to mol/cc
    print("Converting to g/cc 'units'")
    fact = 1./0.602214 
fact /= box[0]*box[1]*dz

ws = numpy.sum(data[:,1:]*mass*fact,axis=1) # get weighted sum
numpy.savetxt(outfile,numpy.column_stack([data[:,0],ws]),fmt='%g') # save output
