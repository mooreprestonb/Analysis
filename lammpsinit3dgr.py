#!/usr/bin/python3
# read lammps data and create 3dgr for water

import sys
import re
import argparse
import numpy
import os.path

def indexline(al,value): # get line of value
    n = 0
    for x in al:
        n += 1
        #print x.rstrip()
        w = x.rstrip().split()
        if w:  # not empty )
            if w[0] == value: 
                break
    if n == len(al):
        print ("Warning: can't fine", value)
        n += 1
#        return -1
#        exit(1)    
    return n
#----------------------------
def getnatoms(al): # get natoms 
    l = 0
    for x in al:
        l += 1
        #print (x.rstrip())
        w = x.rstrip().split()
        if len(w) > 1:  # list has at least 2 elements (and not empty ;-) )
            if w[1] == "atoms": # found keyword "atoms"
                natoms = int(w[0])
                if (l != 3):
                    print ("Atoms not on line 3! ",l)
                    exit(1)
                break
    return(natoms)
#----------------------------
def getatmtypes(lines): #Find number of atoms types
    atmtypes = -1
    for x in lines:
        #print(x.rstrip())
        if(re.search("atom types",x)):
            w = x.rstrip().split()
            atmtypes = int(w[0])
            break
    return(atmtypes)
#----------------------------
def getbox(lines):
# find box dimensions
    l=0
    lbox = -1
    box = []
    for x in lines:
        l += 1
        #print(x.rstrip())
        w = x.rstrip().split()
        if (lbox > -1):  # box line
            lbox += 1
            box.extend([float(w[1])-float(w[0])])
            if (lbox > 2): # read in three lines!
                break
        if (not w and l>3):  # first empty line after line 3
            lbox += 1
    return box
#----------------------------
def getmass(lines):
    nmass = -1
    mass = []
    l=0
    for x in lines:
        l += 1
        if (nmass > -1): # reading in masses
            w = x.rstrip().split()
            if w: # not a blank line
                # print(w,natm)
                if not w[0].isdigit():
                    # print(w[0],"not digit!")
                    break
                nmass += 1
                mass.append([int(w[0]),float(w[1])])
                  
        if(re.search("Masses",x)):
            print("Found Masses: line ",l)
            nmass = 0
    return mass
#----------------------------
def writevmd(i,pos,v1,v2,v3):
    if os.path.isfile("check.vmd"):
        f = open("check.vmd","a")
    else:
        f = open("check.vmd","w")
        f.write("mol new\n")
        
    radius = .2
    f.write("draw color red\n")
    f.write("draw sphere {%f %f %f} radius %f\n" % (pos[i][0], pos[i][1], pos[i][2], radius))
    f.write("draw color white\n")
    f.write("draw sphere {%f %f %f} radius %f\n" % (pos[i+1][0], pos[i+1][1], pos[i+1][2], .5*radius))
    f.write("draw sphere {%f %f %f} radius %f\n" % (pos[i+2][0], pos[i+2][1], pos[i+2][2], .5*radius))
    f.write("draw color %d\n"%(1))
    f.write("draw line {%f %f %f} {%f %f %f}\n" % (pos[i][0], pos[i][1], pos[i][2], pos[i][0]+v1[0], pos[i][1]+v1[1], pos[i][2]+v1[2]))
    f.write("draw color %d\n"%(2))
    f.write("draw line {%f %f %f} {%f %f %f}\n" % (pos[i][0], pos[i][1], pos[i][2], pos[i][0]+v2[0], pos[i][1]+v2[1], pos[i][2]+v2[2]))
    f.write("draw color %d\n"%(3))
    f.write("draw line {%f %f %f} {%f %f %f}\n" % (pos[i][0], pos[i][1], pos[i][2], pos[i][0]+v3[0], pos[i][1]+v3[1], pos[i][2]+v3[2]))
    f.close()

#----------------------------

parser = argparse.ArgumentParser(description="Read lammps init files")
parser.add_argument("infile",help="Lammps init file name to process")
parser.add_argument("outfile",help="Lammps output 3dgr file to create")
args = parser.parse_args()

#lammpsfile = "lammps.test.init"
lammpsfile = args.infile
ofile = args.outfile

print("Reading in Lammps init file",lammpsfile)
f = open(lammpsfile,"r")
#lines = list(f)
lines = f.readlines()
f.close()

Otype = 3
Htype = 4

#Grids
nx = 50
ny = 50
nz = 50

rmin = -10.0
rmax =  10.0

natoms = getnatoms(lines)
print("natoms in header = ",natoms)

atmtypes = getatmtypes(lines)
print("Atoms types = ",atmtypes)

box = getbox(lines)
print("box read in", box)

#find masses
mass = getmass(lines)
nmass = len(mass)
print("Masses found = ",nmass)

if (nmass != atmtypes):  # error check
    print("Something wrong, atom types do not match", atmtypes, " masses found !=", nmass)
    exit()

delta = [(rmin-rmax)/nx,(rmin-rmax)/ny,(rmin-rmax)/nz]

pos = numpy.zeros((natoms,3))
atype = numpy.zeros((natoms))
#read natoms
natm = -1
hl = 0
nl = 0
for x in lines:
    w = x.rstrip().split()
    nl += 1
    if w: # not a blank line
        # print(w,natm)
        if (natm > -1):
            if not w[0].isdigit():
                # print(w[0],"not digit!")
                break
            natm += 1
            indx = int(w[0])  # atom index
            gindx = int(w[1]) # group index
            pos[indx-1] = [float(w[4]),float(w[5]),float(w[6])]
            atype[indx-1] = int(w[2])
        if w[0] == "Atoms": # start counting Atoms!
            natm = 0
            hl = nl
print("natoms found = ",natm)
if (natm != natoms):  # error check
    print("Something wrong, natoms do not match", natm, " found !=", natoms)
    exit()

gr3d = numpy.zeros((nx,ny,nz))
#print(gr3d.shape,gr3d)

oidx = numpy.where(atype == Otype)[0] # oxygen type index
#for i in range(natoms):
#    if atype[i] == Otype:
#        oidx.append(i)
print("Number of Otypes found = ",oidx.size)

for i in range(oidx.size-1):
    ii = oidx[i]
    r1 = pos[ii+1]-pos[ii] # H1-O
    r2 = pos[ii+2]-pos[ii] # H2-O
    v1 = (r1+r2)
    v1 /= numpy.linalg.norm(v1)
    v2 = numpy.cross(v1,r1)
    v2 /= numpy.linalg.norm(v2)
    v3 = numpy.cross(v1,v2)
    v3 /= numpy.linalg.norm(v3)
    # print(i,v1,v2,v3)
    #writevmd(ii,pos,v1,v2,v3)
    for j in range(i+1,oidx.size):
        jj = oidx[j]
        dr = pos[jj]-pos[ii]
        dx = numpy.dot(dr,v1)
        dy = numpy.dot(dr,v2)
        dz = numpy.dot(dr,v3)
        print(dx,dy,dz)
    
#open file to write
print("Creating file :",ofile)
f = open(ofile,"w")

#print header
hdr = "# opendx file for testing with 3dgr\n"
f.write(hdr)
hdr = "object 1 class gridpositions counts "+str(nx)+" "+str(ny)+" "+str(nz)+"\n"
f.write(hdr)
hdr = "origin " + str(rmin) + " " + str(rmin) + " " + str(rmin) +"\n"
f.write(hdr)
hdr = "delta "+str(delta[0])+" 0 0\ndelta 0 "+str(delta[1])+" 0\ndelta 0 0 "+str(delta[2])+"\n"
f.write(hdr)
hdr = "object 2 class gridconnections counts "+str(nx)+" "+str(ny)+" "+str(nz)+"\n"
f.write(hdr)
hdr = "object 3 class array type double rank 0 items "+str(nx*ny*nz)+" data follows\n"
f.write(hdr)
print("Wrote header")
exit(1)
