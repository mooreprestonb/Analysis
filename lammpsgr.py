#!/usr/bin/python3
# calculate g(r) of all types from lammps traj

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
        print("ERROR! natoms can not change!",natoms,natms," config",config)
        exit(1)
            
    line = f.readline() # read "ITEM: BOX BOUNDS pp pp pp"
    if (line.rstrip()[:16] != "ITEM: BOX BOUNDS"):
        print("ERROR! 5th line not \"ITEM: BOX BOUNDS pp pp pp\"")
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
    
    line = f.readline() # read "ITEM: ATOMS id type xs ys zs" (scaled coord)
    if (line.rstrip() == "ITEM: ATOMS id type xs ys zs"):
        ibox = 1
    elif (line.rstrip() == "ITEM: ATOMS id type x y z ix iy iz"):
        ibox = 2
    else:
        print("ERROR! 3th line not \n\"ITEM: ATOMS id type xs ys zs\" or ")
        print("\"ITEM: ATOMS id type x y z ix iy iz\"")
        exit(1)
        
    # loop over atoms
    for i in range(natoms):
        line = f.readline()
        data = line.split()
        num = int(data[0])-1 # lammps id's start with 1
        pos[num][0] = float(data[2])
        pos[num][1] = float(data[3])
        pos[num][2] = float(data[4])
        if(ibox ==1):
            pos[num] *= box
        if (config==0):
            atypes[num] = int(data[1])
        else:
            if(atypes[num] != int(data[1])):
                print("ERROR! atom type change on",num+1)
    return 1 # read in configuration
# end readlammpsconfig

#------------------------------------------------------------
def savehist(outfile,hist,nconfig,hdr,hnorm):
    data = numpy.array(hist*hnorm)
    hdr1 = hdr + " nconfigs: " + str(nconfig)
    numpy.savetxt(outfile,data,fmt='%g',header=hdr1)

#------------------------------------------------------------
# read command line arg using argparser
parser = argparse.ArgumentParser(description="Read in Lammps trajectory file and calculate g(r) for types")
parser.add_argument('input', help='input lammpstrj file')
parser.add_argument('output', help='output data file name')
parser.add_argument('-nbins',type=int,dest='nbins',default=100,help="Number of bins to use in histogram")
parser.add_argument('-time',type=int,dest='time',default=10,help="Report and write progress every so many seconds")
parser.add_argument('-max',type=float,dest='rmax',default=10,help="Max distancs in g(r)")
parser.add_argument('-min',type=float,dest='rmin',default=0,help="Min Distance in g(r)")

# read in arguments, now translate to variables
args = parser.parse_args()
nbins = args.nbins
configname = args.input # 'test.lammpstrj'
outfile = args.output
itime = args.time
rmax = args.rmax
rmin = args.rmin

print("Processing",configname,"to",outfile)
natoms = getlammpsatoms(configname)
print("Number of atoms:",natoms)

types = getlammpstypes(configname,natoms)
ntypes =len(types.keys())
print("Number of atom types:",ntypes,types)

box = getlammpsbox(configname)
print("Initial box dimensions = ",box)

nconf = getlammpsnconf(configname,natoms) # count configurations

# allocate arrays
pos = numpy.zeros((natoms,3))
atypes = numpy.zeros(natoms,dtype=int)
ninter = int(ntypes*(ntypes-1)/2+ntypes)
hist = numpy.zeros((nbins,ninter+1)) # types 1-1, 1-2 ... 1-n, 2-2.... 2-n, ... n-n
histnorm = numpy.zeros((nbins,ninter+1)) 
dx = (rmax-rmin)/(nbins-1)
hist[:,0] = numpy.arange(rmin,dx*nbins+rmin,dx) # create bin values

ioff = numpy.zeros(ntypes,dtype=int)
ioff[0] = 1
for i in range(1,ntypes):
    ioff[i] = ioff[i-1]+ntypes-i
#print(ioff)
    
# create norm array for historgram
fact = 1./(4.*math.pi*dx)
for i in range(nbins):
    histnorm[i][0] = 1
    xpt = (i+.5)*dx+rmin
    for j in range(ntypes):
        for k in range(j,ntypes):
            joff = ioff[j] + k
            histnorm[i][joff] = fact/(types[j+1]*types[k+1]*xpt*xpt)
            # print(i,j,k,joff,fact,types[j+1],types[k+1],xpt)
            if(j==k) :
                histnorm[i][joff] *= 2.  # diagonal terms need a factor of 2
#print(fact,histnorm[:2,:])

# create header for file
hdr = ""
for i in range(len(sys.argv)):
    hdr += " " + sys.argv[i]
hdr += "\n dz = " + str(dx) + " Types:" + str(ntypes)
for k in sorted(types): # iterate over key, value pairs 
    hdr += " " + str(k) + ":" + str(types[k])
hdr += "\n time "
for j in range(ntypes):
    for k in range(j,ntypes):
        hdr += str(j+1) + "-" + str(k+1) + " " 
f = open(configname,'r')
nconfig = 0  # read initial configuration
tnow = time.time()
ttime = tnow
svol = 0
svol2 = 0
print("Processing ",configname)
while(readconfig(f,natoms,box,pos,atypes,nconfig)):
    #print(natoms,box,pos,atypes)
    nconfig += 1
    print("Processing config: ",nconfig,box)
    svol += box[0]*box[1]*box[2]
    svol2 += box[0]*box[1]*box[2]*box[0]*box[1]*box[2]

    for i in range(natoms-1):
        dist = pos[i+1:]-pos[i]
        dist = dist-numpy.floor(dist/box+.5)*box # get periodic boundaries
        #for j in range(len(dist)):
        #    dist[k] -= box[k]*(math.floor(dist[k]/box[k]+.5))
        rdist = numpy.linalg.norm(dist,axis=1)
        ibin = ((rdist-rmin)/dx).astype(int)
        for j in range(len(rdist)):
            if(ibin[j]>=0 and ibin[j]<nbins-1):
                jn = atypes[i]-1
                jm = atypes[i+1+j]-1
                if(jn>jm): # make sure we are in correct values
                    temp = jn
                    jn = jm
                    jm = temp
                jbin = ioff[jn] + jm
                hist[ibin[j]][jbin] += 1
                #print(i,i+1+j,pos[i],pos[i+1+j],rdist[j],ibin[j],jbin,ntypes,jn,jm)
                #frac = (rdist[j]-rmin)/dx - ibin[j]
                #hist[ibin[j]  ][jbin] += 1.-frac
                #hist[ibin[j]+1][jbin] += frac                
            #else:
                #print("Warning: atom position out of range",i,j,pos[i],pos[j],rdist)
                #ibin = max(ibin,0)
                #ibin = min(ibin,nbins-1)

    tnowi = time.time()
    if(itime < tnowi-tnow):
        runtime = tnowi-ttime
        etime = runtime*nconf/nconfig
        avol = svol/nconfig
        print('configs = {}/{} ~ {:2.1%}, time={:g}/{:g}  <vol> = {}'.format(nconfig,nconf,nconfig/nconf,runtime,etime,avol))
        fact = numpy.array(histnorm)
        fact[:,1:] *= avol/nconfig
        savehist(outfile,hist,nconfig,hdr,fact)
        tnow = time.time()

avol = svol/nconfig
stdvol = math.sqrt((svol2/nconfig - avol*avol))
print("# configs read in:",nconfig," <vol> =",avol, "stdvol = ",stdvol)
#print(hist)
histnorm[:,1:] *= avol/nconfig
savehist(outfile,hist,nconfig,hdr,histnorm)
