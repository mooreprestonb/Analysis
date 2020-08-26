#!/usr/bin/python3
# readin lammps configuration

import sys
import argparse
import numpy

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
        print("ERROR! natoms can not change!",natoms,natms)
        exit(1)
            
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

# read command line arg using argparser
parser = argparse.ArgumentParser(description="Read in Lammps trajectory file")
parser.add_argument('input', help='input lammpstrj file')
parser.add_argument('output', help='output data file name')

# readin arguments, now translate to variables
args = parser.parse_args()

configname = args.input # 'test.lammpstrj'
natoms = getlammpsatoms(configname)
print("Number of atoms:",natoms)

pos = numpy.zeros((natoms,3))
atypes = numpy.zeros(natoms,dtype=int)
box = numpy.zeros(3)

f = open(configname,'r')
nconfig = 0  # read initial configuration
while(readconfig(f,natoms,box,pos,atypes,nconfig)):
    nconfig += 1
    #print(natoms,box,pos,atypes)
print("# configs read in:",nconfig)
