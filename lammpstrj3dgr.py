#!/usr/bin/python3
# calculate g(r) of all types from lammps traj

import sys
import argparse
import numpy
import math
import time
import subprocess
import os.path
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
def getlammpsatypes(configname,natoms,atypes): # find type of each atom in lammps traj first config
    f = open(configname,'r')
    for i in range(9): # readover header
        lines = f.readline()
    for i in range(natoms): # readover  header
        words = f.readline().split()
        indx = int(words[0])-1 # atom index
        itype = int(words[1]) # get atom type
        atypes[indx] = itype
    f.close()
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
def getngrblist(configname): # find interacting atoms and types
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
    # should use fancy array indexing to get only what we want :-)
    # for i in range(natoms):
    #     if(ioffmat[atypes[j]][0] == -1):
    #          indexarray.append(j)
    #          jtypes.append(atypes[j])

    return(indexarray,jypes)
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
        if(atypes[num] != int(data[1])):
            print("ERROR! atom type changed on",num+1," Config:",config)
            exit(1)
    return 1 # read in configuration
# end readconfig
#------------------------------------------
def getvectwater(ii,pos,box,v):
    r1 = pos[ii+1]-pos[ii] # H1-O
    r2 = pos[ii+2]-pos[ii] # H2-O
    r1 -= numpy.floor(r1/box+.5)*box # periodic boundary
    r2 -= numpy.floor(r2/box+.5)*box # periodic boundary
    v[0] = r2-r1 # H2-H1
    v[0] /= numpy.linalg.norm(v[0])
    v[2] = numpy.cross(v[0],(r1+r2))
    v[2] /= numpy.linalg.norm(v[2])
    v[1] = numpy.cross(v[2],v[0])
    v[1] /= numpy.linalg.norm(v[1])
#------------------------------------------
def writevmd(i,pos,v):
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
    f.write("draw line {%f %f %f} {%f %f %f}\n" % (pos[i][0], pos[i][1], pos[i][2], pos[i][0]+v[0][0], pos[i][1]+v[0][1], pos[i][2]+v[0][2]))
    f.write("draw color %d\n"%(2))
    f.write("draw line {%f %f %f} {%f %f %f}\n" % (pos[i][0], pos[i][1], pos[i][2], pos[i][0]+v[1][0], pos[i][1]+v[1][1], pos[i][2]+v[1][2]))
    f.write("draw color %d\n"%(3))
    f.write("draw line {%f %f %f} {%f %f %f}\n" % (pos[i][0], pos[i][1], pos[i][2], pos[i][0]+v[2][0], pos[i][1]+v[2][1], pos[i][2]+v[2][2]))
    f.close()
#----------------------------
def bingr3d(pt,np,bins,gr3d):
    ptn = (pt-bins[0])/bins[2] # get bin numbers
    nr = ptn.astype(int) # get integer numbers
    ptn -= nr # get fraction of bin
    #print(pt,ptn,nr)
    for i in range(pt.shape[0]):
        if (numpy.amax(nr[i]-np)<-1 and numpy.amin(nr[i])>-1): # in range?
            gr3d[nr[i][0]][nr[i][1]][nr[i][2]] += (1-ptn[i][0])+(1-ptn[i][1])+(1-ptn[i][2])
    
            gr3d[nr[i][0]+1][nr[i][1]][nr[i][2]] += (ptn[i][0])+(1-ptn[i][1])+(1-ptn[i][2])
            gr3d[nr[i][0]][nr[i][1]+1][nr[i][2]] += (1-ptn[i][0])+(ptn[i][1])+(1-ptn[i][2])
            gr3d[nr[i][0]][nr[i][1]][nr[i][2]+1] += (1-ptn[i][0])+(1-ptn[i][1])+(ptn[i][2])
        
            gr3d[nr[i][0]+1][nr[i][1]+1][nr[i][2]] += (ptn[i][0])+(ptn[i][1])+(1-ptn[i][2])
            gr3d[nr[i][0]+1][nr[i][1]][nr[i][2]+1] += (ptn[i][0])+(1-ptn[i][1])+(ptn[i][2])
            gr3d[nr[i][0]][nr[i][1]+1][nr[i][2]+1] += (1-ptn[i][0])+(ptn[i][1])+(ptn[i][2])
        
            gr3d[nr[i][0]+1][nr[i][1]+1][nr[i][2]+1] += (ptn[i][0])+(ptn[i][1])+(ptn[i][2])

#-----------------------------------------------------------
def writedxfile(ofile,np,bins,gr3d,fact):
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
                if(i==np[0]-1 or j == np[1]-1 or k == np[1]-1):
                    f.write(str(gr3d[i][j][k]*fact*2.)+"\n") # scale last points by 2
                else:
                    f.write(str(gr3d[i][j][k]*fact)+"\n")
                    
    f.write("\nobject density class field\n")
    print("Wrote file",ofile)

#------------------------------------------------------------
# read command line arg using argparser
parser = argparse.ArgumentParser(description="Read in Lammps trajectory file and calculate g(r) for types")
parser.add_argument('input', help='input lammpstrj file')
parser.add_argument('output', help='output data file name')
parser.add_argument('-nbins',type=int,dest='nbins',default=20,help="Number of bins in each dimenstion")
parser.add_argument('-time',type=int,dest='time',default=10,help="Report and write progress every so many seconds")
parser.add_argument('-max',type=float,dest='rmax',default=10,help="Max distancs in g(r)")
parser.add_argument('-min',type=float,dest='rmin',default=-10,help="Min distance in g(r)")
parser.add_argument('-type1', type=int, dest='type1',default=1,help="Water Oxygen atoms type")
parser.add_argument('-type2', type=int, dest='type2',default=1,help="type to bin against water")

# read in arguments, now translate to variables
args = parser.parse_args()
print (args)
nbins = args.nbins
configname = args.input # 'test.lammpstrj'
outfile = args.output
itime = args.time
rmax = args.rmax
rmin = args.rmin
otype = args.type1
type2 = args.type2

#Grids
np = numpy.array([nbins,nbins,nbins])
bins = numpy.zeros((3,3))
bins[0] = rmin
bins[1] = rmax
bins[2] = (bins[1]-bins[0])/np
print(bins)

print("Processing",configname,"to",outfile)
natoms = getlammpsatoms(configname)
print("Number of atoms:",natoms)

types = getlammpstypes(configname,natoms)
ntypes =len(types.keys())
stypes = str(ntypes)
for k in sorted(types): # iterate over key, value pairs 
    stypes += " " + str(k) + ":" + str(types[k])
print("Number of atom types: ",stypes)

box = getlammpsbox(configname)
print("Initial box dimensions = ",box)

nconf = getlammpsnconf(configname,natoms) # count configurations

# check types?
    
# allocate arrays
pos = numpy.zeros((natoms,3))
atypes = numpy.zeros(natoms,dtype=int)
getlammpsatypes(configname,natoms,atypes)

gr3d = numpy.zeros((np))
#print(gr3d.shape,gr3d)

oidx = numpy.where(atypes == otype)[0] # oxygen type index
indx2 = numpy.where(atypes == type2)[0] # type2 index
noxy = oidx.size
ntype2 = indx2.size
print("Number of Otypes found = ",noxy,"Number of type 2 =",ntype2)
print(oidx,indx2)

v = numpy.zeros((3,3))
rpos = numpy.zeros((ntype2,3))

f = open(configname,'r')
nconfig = 0  # 
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

    pos2 = numpy.take(pos,indx2,axis=0)
    for i in range(noxy):
        ii = oidx[i]
        getvectwater(ii,pos,box,v) # get verticies of water molecule
        #print(i,noxy,ntype2,nconfig)
        #print(ii,pos[ii],pos[ii+1],pos[ii+2])
        #writevmd(ii,pos,v)
        rpos = pos2-pos[ii] # distance between O and all other O's
        if(otype == type2):
            rpos = numpy.delete(rpos,i,axis=0) # remove self
        rpos -= numpy.floor(rpos/box+.5)*box # periodic boundary
        rpos = numpy.inner(rpos,v) # distance along each vecticies
        bingr3d(rpos,np,bins,gr3d) # great density plot

        tnowi = time.time()
        if(itime < tnowi-tnow):
            runtime = tnowi-ttime
            etime = runtime*nconf/nconfig
            avol = svol/nconfig
            print('otype = {}, configs = {}/{} ~ {:2.1%}, time={:g}/{:g}  <vol> = {}'.format(i,noxy,nconfig,nconf,nconfig/nconf,runtime,etime,avol))
            fact = 1./(3*noxy*nconfig*bins[2][0]*bins[2][1]*bins[2][2])
            writedxfile(outfile,np,bins,gr3d,fact)
            #print(numpy.amin(gr3d)*fact,numpy.amax(gr3d)*fact)
            tnow = time.time()

avol = svol/nconfig
stdvol = math.sqrt((svol2/nconfig - avol*avol))
print("# configs read in:",nconfig," <vol> =",avol, "stdvol = ",stdvol)
print("Otypes:",otype,"with ",noxy,"molecules and type2",type2," with",ntype2)

fact = 1./(3*noxy*nconfig*bins[2][0]*bins[2][1]*bins[2][2])
writedxfile(outfile,np,bins,gr3d,fact)
print(numpy.amin(gr3d)*fact,numpy.amax(gr3d)*fact)
