#!/usr/bin/python
# read lammps data and put back into box .

import sys
import re
import argparse
import numpy as np
from math import floor

def indexline(al,value):
    n = 0
    for x in al:
        n += 1
        #print x.rstrip()
        w = x.rstrip().split()
        if w:  # not empty )
            if w[0] == value: 
                break
    if n == len(al):
        print "Warning: can't fine", value
        n += 1
        #        exit(1)    
    return n
#----------------------------

parser = argparse.ArgumentParser(description="Periodic lammps init files")
parser.add_argument("infile",help="Lammps file name to process")
parser.add_argument("outfile",help="Lammps output file to create")
args = parser.parse_args()

#lammpsfile = "lammps.test.init"
lammpsfile = args.infile
ofile = args.outfile

print "Reading in Lammps file",lammpsfile

f = open(lammpsfile,"r")
lines = list(f)
f.close()
nlines = len(lines)

print "Read in ",nlines," lines"
f = open(ofile,"w")

l=0
#Find natoms
for x in lines:
    l += 1
    #print x.rstrip()
    w = x.rstrip().split()
    if len(w) > 1:  # list has at least 2 elements (and not empty ;-) )
        if w[1] == "atoms": # found keyword "atoms"
            natoms = int(w[0])
            if (l != 3):
                print "Atoms not on line 3! ",l
                exit(1)
            break
print "natoms in header = ",natoms

#Find atoms types
for x in lines:
    l += 1
    #print x.rstrip()
    if(re.search("atom types",x)):
        w = x.rstrip().split()
        atmtypes = int(w[0])
        break
print "Atoms types = ",atmtypes

# find box dimensions
l=0
lbox = -1
box = []
for x in lines:
    l += 1
    #print x.rstrip()
    w = x.rstrip().split()
    if (lbox > -1):  # box line
        lbox += 1
        box.extend([float(w[0]),float(w[1])])
        if (lbox > 2): # read in three lines!
            break
    if (not w and l>3):  # first empty line after line 3
        lbox += 1
print "box read in", box, l

for i in range(3):
    bmid = (box[i*2+1]-box[i*2])
    box[2*i] = -bmid/2
    box[2*i+1] = bmid/2
print "Writing new box: ",box

lines[l-3] = str(box[0]) + " " + str(box[1]) + " xlo xhi\n"
lines[l-2] = str(box[2]) + " " + str(box[3]) + " ylo yhi\n"
lines[l-1] = str(box[4]) + " " + str(box[5]) + " zlo zhi\n"

#find masses
nmass = -1
mass = []
for x in lines:
    if (nmass > -1): # reading in masses
        w = x.rstrip().split()
        if w: # not a blank line
            # print w,natm
            if not w[0].isdigit():
                # print w[0],"not digit!"
                break
            nmass += 1
            mass.append([int(w[0]),float(w[1])])
        
    if(re.search("Masses",x)):
        print "found Masses"
        nmass = 0
    
print "Masses found = ",nmass
#print mass, mass[0][1]
if (nmass != atmtypes):  # error check
    print "Something wrong, atom types do not match", atmtypes, " masses found !=", nmass
    exit()

#read natoms
natm = -1
atms = []
indx = []
gindx = []
hl = 0
nl = 0
for x in lines:
    w = x.rstrip().split()
    nl += 1
    if w: # not a blank line
        # print w,natm
        if (natm > -1):
            if not w[0].isdigit():
                # print w[0],"not digit!"
                break
            natm += 1
            nbs = []
            nbs.append(int(w[0])) # atom number
            nbs.append(int(w[1])) # molecule number
            nbs.append(int(w[2])) # atom type
            nbs.append(float(w[3])) # charge
            nbs.append(float(w[4])) # x
            nbs.append(float(w[5])) # y
            nbs.append(float(w[6])) # z
            nbs.append(int(w[7])) # cells in x
            nbs.append(int(w[8])) # cells in y
            nbs.append(int(w[9])) # cells in z
            atms.append(nbs) # atom line
            indx.append(int(w[0]))  # atom index
            gindx.append(int(w[1])) # group index
        if w[0] == "Atoms": # start counting Atoms!
            natm = 0
            hl = nl
print "natoms found = ",natm

if (natm != natoms):  # error check
    print "Something wrong, natoms do not match", natm, " found !=", natoms
    exit()

#print atms

indxs = np.argsort(indx) # index to sort array
gindxs = np.argsort(gindx) # group index to sort array
sgindx = sorted(set(gindx)) # sorted group index

#print indx, indxs
#print gindx, sgindx

#print header
for i in range(hl):
    f.write(lines[i])
f.write("\n")
print "Wrote header, getting molecule com and puting with box"
print "Molecule #:",
for i in range(gindx[-1]):
    print i,
    gatm = []
    for j in range(natoms):
        if(atms[j][1]==i+1): # get atoms with molecule number i+1
            gatm.append(j)

    tmass = 0
    com = np.zeros(3)
    pos = np.empty((len(gatm),3))
    for j in range(len(gatm)): # get positions and com
        k = gatm[j]
        gmass = mass[atms[k][2]-1][1]
        tmass += gmass
        pos[j][0] = atms[k][4]+2*box[1]*atms[k][7] # lammps writes number of box offsets
        pos[j][1] = atms[k][5]+2*box[3]*atms[k][8]
        pos[j][2] = atms[k][6]+2*box[5]*atms[k][9]
        com[0] += pos[j][0]*gmass
        com[1] += pos[j][1]*gmass
        com[2] += pos[j][2]*gmass
    for l in range(3):
        com[l] /= tmass;
    # print pos," ",com
    for l in range(3): # com becomes offset
        com[l] = floor(com[l]/(box[l*2+1]*2)+.5)*box[l*2+1]*2
    for j in range(len(gatm)): # get new positions and write out
        k = gatm[j] # atom index
        for l in range(3):
            pos[j][l] -= com[l]                
        sl = str(atms[k][0]) + " " + str(atms[k][1]) + " " + str(atms[k][2]) + " " + str(atms[k][3]) + " "
        sl += str(pos[j][0]) + " " + str(pos[j][1]) + " " + str(pos[j][2]) + " 0 0 0\n"
        f.write(sl)
f.write("\n")


print "Wrote atoms within box"

# Now for other data
nl = indexline(lines,"Atoms")
vl = indexline(lines,"Velocities")
bl = indexline(lines,"Bonds")
al = indexline(lines,"Angles")
dl = indexline(lines,"Dihedrals")
il = indexline(lines,"Impropers")

print nl,vl,bl,al,dl,il

if (vl > nlines):
    print "No Velocities"
else:
    f.write("Velocities\n\n")
    for i in range(vl+1,bl-1):
        f.write(lines[i])
f.write("\n")
print "Wrote Velocities"
                
f.write("Bonds\n\n")
for i in range(bl+1,al-1):
    f.write(lines[i])
f.write("\n")
print "Wrote Bonds"

f.write("Angles\n\n")
for i in range(al+1,dl-1):
    f.write(lines[i])
f.write("\n")
print "Wrote Angles"

f.write("Dihedrals\n\n")
for i in range(dl+1,il-1):
    f.write(lines[i])
f.write("\n")
print "Wrote Dihedrals"

f.write("Impropers\n\n")
for i in range(il+1,len(lines)):
    f.write(lines[i])
f.write("\n")
