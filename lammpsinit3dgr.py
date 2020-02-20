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
    natoms=-1
    for x in al:
        if(re.search("atoms",x)):
            #print (x.rstrip())
            w = x.rstrip().split()
            natoms = int(w[0])
            break
    if natoms == -1:
        print("Atoms keyword not found!")
        exit(1)
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
    if atmtypes == -1:
        print("Atoms types keyword not found!")
        exit(1)
    return(atmtypes)
#----------------------------
def getbox(lines):  # find box dimensions
    box = numpy.array([-1,-1,-1])
    for x in lines:
        w = x.rstrip().split()
        #print (x,w)
        if(re.search("xlo",x)):
            box[0] = float(w[1])-float(w[0])
        if(re.search("ylo",x)):
            box[1] = float(w[1])-float(w[0])
        if(re.search("zlo",x)):
            box[2] = float(w[1])-float(w[0])
        if (numpy.amin(box) > 0):
            print("Box found")
            break
    if (numpy.amin(box) < 0):
        print("Box not found!",box)
        exit(1)        
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
def readinitatoms(natoms,lines,pos,atype):
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
                #gindx = int(w[1]) # group index
                atype[indx-1] = int(w[2])
                #chg[indx-1] = float(w[3]) # atom charge
                pos[indx-1] = [float(w[4]),float(w[5]),float(w[6])]
                #boxoff[indx-1] = [int(w(7),int(w(8),int(w(9)]
            if w[0] == "Atoms": # start counting Atoms!
                natm = 0
                hl = nl
    print("natoms found = ",natm)
    if (natm != natoms):  # error check
        print("Something wrong, natoms do not match", natm, " found !=", natoms)
        exit()
#------------------------------------------
def getvectwater(ii,pos,v):
    r1 = pos[ii+1]-pos[ii] # H1-O
    r2 = pos[ii+2]-pos[ii] # H2-O
    v[0] = r2-r1 # H2-H1
    v[0] /= numpy.linalg.norm(v[0])
    v[2] = numpy.cross(v[0],(r1+r2))
    v[2] /= numpy.linalg.norm(v[2])
    v[1] = numpy.cross(v[2],v[0])
    v[1] /= numpy.linalg.norm(v[1])
#------------------------------------------
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
def bingr3d(pt,np,bins,gr3d):
    ptn = (pt-bins[0])/bins[2] # get bin numbers
    nr = ptn.astype(int) # get integer numbers
    #print(pt,ptn,nr)
    ptn -= nr # get fraction of bin
    for i in range(pt.shape[0]):
        if (numpy.amax(nr[i]-np)< -2 and numpy.amin(nr[i]+np)>0): # in range?
            gr3d[nr[i][0]][nr[i][1]][nr[i][2]] += (1-ptn[i][0])+(1-ptn[i][1])+(1-ptn[i][2])
    
            gr3d[nr[i][0]+1][nr[i][1]][nr[i][2]] += (ptn[i][0])+(1-ptn[i][1])+(1-ptn[i][2])
            gr3d[nr[i][0]][nr[i][1]+1][nr[i][2]] += (1-ptn[i][0])+(ptn[i][1])+(1-ptn[i][2])
            gr3d[nr[i][0]][nr[i][1]][nr[i][2]+1] += (1-ptn[i][0])+(1-ptn[i][1])+(ptn[i][2])
        
            gr3d[nr[i][0]+1][nr[i][1]+1][nr[i][2]] += (ptn[i][0])+(ptn[i][1])+(1-ptn[i][2])
            gr3d[nr[i][0]+1][nr[i][1]][nr[i][2]+1] += (ptn[i][0])+(1-ptn[i][1])+(ptn[i][2])
            gr3d[nr[i][0]][nr[i][1]+1][nr[i][2]+1] += (1-ptn[i][0])+(ptn[i][1])+(ptn[i][2])
        
            gr3d[nr[i][0]+1][nr[i][1]+1][nr[i][2]+1] += (ptn[i][0])+(ptn[i][1])+(ptn[i][2])
    
#----------------------------

parser = argparse.ArgumentParser(description="Read lammps init and traj file to create 3dgr of water")
parser.add_argument("infile1",help="Lammps init file name to process")
parser.add_argument("infile2",help="Lammps trajector file name to process")
parser.add_argument("outfile",help="Lammps output 3dgr file to create")
parser.add_argument('-max',type=float,dest='rmax',default=10,help="Max distancs in 3D g(r)")
parser.add_argument('-min',type=float,dest='rmin',default=-10,help="Min Distance in 3D g(r)")
parser.add_argument("-otype",type=int,dest='otype',default=1,help="Oxygen type in lammps file")
parser.add_argument("-npt",type=int,dest='npt',default=20,help="Number of point in each dimenstion")
args = parser.parse_args()

#lammpsfile = "lammps.test.init"
lammpsinitfile = args.infile1
lammpstrajfile = args.infile2
ofile = args.outfile
Otype = args.otype
rmin = args.rmin
rmax = args.rmax
npt = args.npt

print("Reading in Lammps init file",lammpsinitfile)
f = open(lammpsinitfile,"r")
#lines = list(f)
lines = f.readlines()
f.close()

#Grids
np = numpy.array([npt,npt,npt])
bins = numpy.zeros((3,3))
bins[0] = rmin
bins[1] = rmax
bins[2] = (bins[1]-bins[0])/np
print(bins)

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

pos = numpy.zeros((natoms,3))
atype = numpy.zeros((natoms))

readinitatoms(natoms,lines,pos,atype)

gr3d = numpy.zeros((np))
#print(gr3d.shape,gr3d)

oidx = numpy.where(atype == Otype)[0] # oxygen type index
#for i in range(natoms):
#    if atype[i] == Otype:
#        oidx.append(i)
print("Number of Otypes found = ",oidx.size)

v = numpy.zeros((3,3))
opos = numpy.take(pos,oidx,axis=0)
rpos = numpy.zeros((oidx.size,3))

for i in range(oidx.size):
    ii = oidx[i]
    getvectwater(ii,pos,v) # get verticies of water molecule
    print(i,oidx.size)
    # writevmd(ii,pos,v)
    rpos = opos-pos[ii] # distance between O and all other O's
    rpos = numpy.delete(rpos,i,axis=0) # remove self
    rpos -= numpy.floor(rpos/box+.5)*box # periodic boundary
    rpos = numpy.inner(rpos,v) # distance along each vecticies
    bingr3d(rpos,np,bins,gr3d) # great density plot

#open file to write
print("Creating file :",ofile)
f = open(ofile,"w")
#print header
hdr = "# opendx file for testing with 3dgr\n"
f.write(hdr)
hdr = "object 1 class gridpositions counts "+str(np[0])+" "+str(np[1])+" "+str(np[2])+"\n"
f.write(hdr)
hdr = "origin " + str(bins[0][0]) + " " + str(bins[0][1]) + " " + str(bins[0][2]) +"\n"
f.write(hdr)
hdr = "delta "+str(bins[2][0])+" 0 0\ndelta 0 "+str(bins[2][1])+" 0\ndelta 0 0 "+str(bins[2][2])+"\n"
f.write(hdr)
hdr = "object 2 class gridconnections counts "+str(np[0])+" "+str(np[1])+" "+str(np[2])+"\n"
f.write(hdr)
hdr = "object 3 class array type double rank 0 items "+str(np[0]*np[1]*np[2])+" data follows\n"
f.write(hdr)
#write data
fact=1./oidx.size
for i in range(np[0]):
    for j in range(np[1]):
        for k in range(np[2]):
            f.write(str(gr3d[i][j][k]*fact)+"\n")
f.write("\nobject density class field\n")
print("Wrote file",ofile)
print(numpy.amin(gr3d)*fact,numpy.amax(gr3d)*fact)
