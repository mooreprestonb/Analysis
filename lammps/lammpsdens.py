#!/usr/bin/python3
# calculate density of types from lammps traj

import sys
import argparse
import numpy
import math
import time
import subprocess
#------------------------------------------------------------
def getlammpsbox(configname): # box in lammps traj
    f = open(configname,'r')

    box = numpy.zeros(3)
    line = f.readline() # read header or EOF
    if (line==""): # end of file! EOF!
        return 0
    if (line.rstrip() != "ITEM: TIMESTEP"):
        print("ERROR! first line not \"ITEM: TIMESTEP\"")
        exit(1)
        
    line = f.readline() # read timestep
    ts = int(line)

    line = f.readline() # read "ITEM: NUMBER OF ATOMS"
    if (line.rstrip() != "ITEM: NUMBER OF ATOMS"):
        print("ERROR! 3rd line not \"ITEM: NUMBER OF ATOMS\"")
        exit(1)
    line = f.readline() # natoms
    natms = int(line)
            
    line = f.readline() # read "ITEM: BOX BOUNDS pp pp pp"
    if (line.rstrip()[:16] != "ITEM: BOX BOUNDS"):
        print("ERROR! 3rd line not \"ITEM: BOX BOUNDS pp pp pp\"")
        exit(1)
    line = f.readline() # xlo xhi
    data = line.split()
    box[0] = float(data[1])-float(data[0])
    line = f.readline() # ylo yhi
    data = line.split()
    box[1] = float(data[1])-float(data[0])
    line = f.readline() # zlo zhi
    data = line.split()
    box[2] = float(data[1])-float(data[0])

    f.close()
    return box

#-----find number of configs (quicker with "wc" but not as portable) ------
def getlammpsnconf(configname,natoms): 
    lines = 0
    print("Counting configures of ",configname)
    blines = subprocess.check_output(["wc","-l",configname])
    # lines = int((str(blines,'utf-8').split())[0])
    lines = int((str(blines.decode()).split())[0])
#    for line in open(configname): # above is 5* faster, but this always works
#        lines += 1
    nconf = int(lines/(natoms+9))
    print ("# of configs = ",nconf)
    return (nconf)
#------------------------------------------------------------
def getlammpstypes(configname,natoms): # find # types in lammps traj
    types={}
    f = open(configname,'r')
    for i in range(9): # readover header
        lines = f.readline()
    for i in range(natoms): # readover header
        ntype = int(f.readline().split()[1]) # get atom type
        # print(ntype)
        if ntype in types: # add types 
            types[ntype] += 1
        else :
            types[ntype] = 1
    f.close()
    return types

#------------------------------------------------------------ 
def getlammpsatoms(configname): # find # atoms in lammps header
    f = open(configname,'r')

    line = f.readline() # read header or EOF
    if (line==""): # end of file! EOF!
        print("EOF while reading header?!?")
        return 0
    if (line.rstrip() != "ITEM: TIMESTEP"):
        print("ERROR! first line of config is not \"ITEM: TIMESTEP\"")
        exit(1)
    line = f.readline() # read timestep
    ts = int(line)

    line = f.readline() # read "ITEM: NUMBER OF ATOMS"
    if (line.rstrip() != "ITEM: NUMBER OF ATOMS"):
        print("ERROR! 3rd line of config is not \"ITEM: NUMBER OF ATOMS\"")
        exit(1)
    line = f.readline() # natoms
    natms = int(line)
    f.close()
    return(natms)

#------------------------------------------------------------
def readconfig(f,natoms,box,pos,atypes,config): # read lammps configuration

    boxt = numpy.zeros(3)
    
    line = f.readline() # read header or EOF
    if (line==""): # end of file! EOF!
        return 0
    if (line.rstrip() != "ITEM: TIMESTEP"):
        print("ERROR! first line not \"ITEM: TIMESTEP\"")
        exit(1)
        
    line = f.readline() # read timestep
    ts = int(line)

    line = f.readline() # read "ITEM: NUMBER OF ATOMS"
    if (line.rstrip() != "ITEM: NUMBER OF ATOMS"):
        print("ERROR! 3rd line not \"ITEM: NUMBER OF ATOMS\"")
        exit(1)
    line = f.readline() # natoms
    natms = int(line)
    if(natoms != natms):
        print("ERROR! natoms can not change!",natoms,natms)
        exit(1)
            
    line = f.readline() # read "ITEM: BOX BOUNDS pp pp pp"
    if (line.rstrip()[:16] != "ITEM: BOX BOUNDS"):
        print("ERROR! 3rd line not \"ITEM: BOX BOUNDS pp pp pp\"")
        exit(1)
    line = f.readline() # xlo xhi
    data = line.split()
    boxt[0] = float(data[1])-float(data[0])
    line = f.readline() # ylo yhi
    data = line.split()
    boxt[1] = float(data[1])-float(data[0])
    line = f.readline() # zlo zhi
    data = line.split()
    boxt[2] = float(data[1])-float(data[0])
    for i in range(3):
        if(boxt[i] != box[i]):
            print("ERROR! Box has changed!",boxt,"!=",box)
            exit(1)
    
    line = f.readline() # read "ITEM: ATOMS id type xs ys zs" (scaled coord)
    if (line.rstrip() != "ITEM: ATOMS id type xs ys zs"):
        print("ERROR! 3rd line not \"ITEM: ATOMS id type xs ys zs\"")
        exit(1)
        
    # loop over atoms
    for i in range(natoms):
        line = f.readline()
        data = line.split()
        num = int(data[0])-1 # lammps id's start with 1
        pos[num][0] = float(data[2])*box[0]
        pos[num][1] = float(data[3])*box[1]
        pos[num][2] = float(data[4])*box[2]
        if (config==0):
            atypes[num] = int(data[1])
        else:
            if(atypes[num] != int(data[1])):
                print("ERROR! atom type change on",num+1)
    return 1 # read in configuration
# end readlammpsconfig

#------------------------------------------------------------
def savehist(outfile,hist,nconfig,hdr):
    data = numpy.array(hist)
    data[:,1:] /= nconfig # normalize all but first column
    hdr1 = hdr + " nconfigs: " + str(nconfig)
    numpy.savetxt(outfile,data,fmt='%g',header=hdr1)

    # hist[:,1:] /= box[0]*box[1]*dx
# average and normalize
#for i in range(nbins):
#    hist[i][0] = dx*i + hmin
#    for j in range(1,ntypes+1):
#        hist[i][j] /= nconfig
# hist[i][j] /= types[j]
# might want div by box[0]*box[1]*dx to norm density

#------------------------------------------------------------
# read command line arg using argparser
parser = argparse.ArgumentParser(description="Read in Lammps trajectory file and calculate histogram of types along z")
parser.add_argument('input', help='input lammpstrj file')
parser.add_argument('output', help='output data file name')
parser.add_argument('-nbins',type=int,dest='nbins',default=100,help="Number of bins to use in histogram")
parser.add_argument('-time',type=int,dest='time',default=10,help="Report and write progress every so many seconds")

# read in arguments, now translate to variables
args = parser.parse_args()
nbins = args.nbins
configname = args.input # 'test.lammpstrj'
outfile = args.output
itime = args.time

print("Processing",configname,"to",outfile)
natoms = getlammpsatoms(configname)
print("Number of atoms:",natoms)

types = getlammpstypes(configname,natoms)
ntypes =len(types.keys())
print("Number of atom types:",ntypes,types)

box = getlammpsbox(configname)
print("Box dimensions = ",box)

nconf = getlammpsnconf(configname,natoms)

pos = numpy.zeros((natoms,3))
atypes = numpy.zeros(natoms,dtype=int)
hist = numpy.zeros((nbins,ntypes+1)) # types go from 1-n, 0 will be x
hmin = 0 # min in the z (goes from 0 to z)
hmax = box[2]
dx = box[2]/(nbins)
hist[:,0] = numpy.arange(hmin,hmax,dx) # create bin values

# create header for file
hdr = ""
for i in range(len(sys.argv)):
    hdr += " " + sys.argv[i]
hdr += "\n dz = " + str(dx) + " box: " + str(box[0]) + " " + str(box[1]) + " " + str(box[2]) + " Types:" + str(ntypes)
for k in sorted(types): # iterate over key, value pairs 
    hdr += " " + str(k) + ":" + str(types[k])

f = open(configname,'r')
nconfig = 0  # read initial configuration
tnow = time.time()
ttime = tnow
print("Processing ",configname)
while(readconfig(f,natoms,box,pos,atypes,nconfig)):
    #print(natoms,box,pos,atypes)
    nconfig += 1
    #abin = pos[:,2]/dx.astype(int)
    for i in range(natoms):
        #ibin = abin[i]
        pt = pos[i][2]/dx
        ibin = int(pt)
        #print(i,ibin,pos[i][2]/dx,pos[i][2])
        if((ibin<0) or (ibin>(nbins-2))):
            print("Warning: atom position out of range",i,ibin,pos[i][2]/dx,pos[i],box[2])
            ibin = max(ibin,0)
            ibin = min(ibin,nbins-2)
        fr = pt-ibin
        hist[ibin  ][atypes[i]] += 1. - fr
        hist[ibin+1][atypes[i]] += fr
        #hist[ibin][atypes[i]] += 1

    tnowi = time.time()
    if(itime < tnowi-tnow):
        runtime = tnowi-ttime
        etime = runtime*nconf/nconfig
        print('configs = {}/{} ~ {:2.1%}, time={:g}/{:g}'.format(nconfig,nconf,nconfig/nconf,runtime,etime))
        savehist(outfile,hist,nconfig,hdr)
        tnow = time.time()

print("# configs read in:",nconfig)
#print(hist)
savehist(outfile,hist,nconfig,hdr)
