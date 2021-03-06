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
def getngbrlist(configname,natoms,box,pos,atypes,ngbrs): # find interacting atoms and types
    f = open(configname,'r')
    readconfig(f,natoms,box,pos,atypes,0) # read one config to get atypes
    f.close()

    ingbr = [0]*len(ngbrs)
    for i in range(natoms):
        itype = atypes[i]-1
        #print(i,itype,ingbr[itype])
        ngbrs[itype][ingbr[itype]] = i
        ingbr[itype] += 1
        
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

    def boxlength(words):
        return float(words[1]) - float(words[0])
    box[0] = boxlength(f.readline().split()) # xlo xhi
    box[1] = boxlength(f.readline().split()) # ylo yhi
    box[2] = boxlength(f.readline().split()) # zlo zhi
    
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
# end readconfig

#------------------------------------------------------------
def distbin(apos,ipos,jbin): #bin the historgram
    dist = apos-ipos # get distances
    dist = dist-numpy.floor(dist/box+.5)*box # get periodic boundaries
    rdist = numpy.linalg.norm(dist,axis=1)
    ibin = ((rdist-rmin)/dx).astype(int) # get index
#    for j in range(len(rdist)):                
#        if(ibin[j]>=0 and ibin[j]<nbins-1):
#            frac = (rdist[j]-rmin)/dx - ibin[j]
#            hist[ibin[j]  ][jbin] += 1.-frac
#            hist[ibin[j]+1][jbin] += frac                   
    tibin = numpy.where((ibin<nbins-1)*(ibin>0))[0] # get where ibin is > 0 and < nbins-1
#    frac = (rdist[tibin] - rmin)/dx - ibin[tibin]
#    ibx = ibin[tibin]
#    hist[ibx][jbin] += 1 - fact
#    hist[ibx+1][jbin] += fact
#    print(tibin,frac)
    for j in tibin:
        frac1 = (rdist[j]-rmin)/dx - ibin[j]
        hist[ibin[j]  ][jbin] += 1.-frac1
        hist[ibin[j]+1][jbin] += frac1                   
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
parser.add_argument('-itypes', type=int, nargs="+", dest='itypes', help="iatoms types to processes")
parser.add_argument('-jtypes', type=int, nargs="+", dest='jtypes', help="jatoms types to processes")

# read in arguments, now translate to variables
args = parser.parse_args()
print (args)
nbins = args.nbins
configname = args.input # 'test.lammpstrj'
outfile = args.output
itime = args.time
rmax = args.rmax
rmin = args.rmin
aitypes = args.itypes
ajtypes = args.jtypes

print("Processing",configname,"to",outfile)
natoms = getlammpsatoms(configname)
print("Number of atoms:",natoms)

types = getlammpstypes(configname,natoms)
ntypes =len(types.keys())
ngbrs = [None]*ntypes # allocate empty arrays
stypes = str(ntypes)
for k in sorted(types): # iterate over key, value pairs 
    stypes += " " + str(k) + ":" + str(types[k])
    ngbrs[k-1] = numpy.empty(types[k],dtype=int)
print("Number of atom types: ",stypes)

box = getlammpsbox(configname)
print("Initial box dimensions = ",box)
nconf = getlammpsnconf(configname,natoms) # count configurations

# allocate arrays
pos = numpy.zeros((natoms,3))
atypes = numpy.zeros(natoms,dtype=int)
getngbrlist(configname,natoms,box,pos,atypes,ngbrs)

# check types
if (aitypes == None and ajtypes != None):
    print("Error: jtypes specified without itypes")
    exit(1)
    
itypes = []
jtypes = []
if aitypes == None:
    print("Using all types for itypes 1-",ntypes+1)
    for i in range(1,ntypes+1):
        itypes.append(i)
else :
    for it in aitypes:
        if it not in types:
            print("ERROR! itype not found",it,types)
            exit(1)
        itypes.append(it)

if ajtypes == None:
    print("Using all itypes for jtypes")
    jtypes = itypes
else :
    for jt in ajtypes:
        if jt not in types:
            print("ERROR! jtype not found",jt,types)
            exit(1)
        jtypes.append(jt)

print(itypes,jtypes)

# Create offset matrix for types we are processing
toff = 1
ioffmat = numpy.zeros((ntypes+1,ntypes+1),dtype=int)
pairs = []
for i in range(1,ntypes+1):
    if i in itypes: # are we processing this type?
        ioffmat[i][0] = -1
        for j in range(1,ntypes+1):
            if j in jtypes: # are we processing this type?
                ioffmat[0][j] = -1
                if(ioffmat[i][j] == 0):
                    ioffmat[i][j] = toff
                    ioffmat[j][i] = toff # diagonal matrix i->j is same as j->i
                    toff += 1
                    # print("i>=j",i,j,toff)
                    pairs.append(str(i)+"-"+str(j))
print("Offset Matrix\n",ioffmat)
print("Pairs\n",pairs)

dx = (rmax-rmin)/(nbins-1)
#toff = int(ntypes*(ntypes-1)/2+ntypes)+1
histnorm = numpy.zeros((nbins,toff)) 
hist = numpy.zeros((nbins,toff)) # types 1-1, 1-2 ... 1-n, 2-2.... 2-n, ... n-n
#hist[:,0] = numpy.arange(rmin,rmax,dx) # create bin values but can be off in rounding error so we do explicit.
for i in range(nbins): # create xrange for histograms
    r = i*dx+rmin
    hist[i][0] = r

# create norm array for historgram
fact = 1./(4.*math.pi*dx)
for i in range(nbins):
    histnorm[i][0] = 1
    xpt = (i+0.5)*dx+rmin
    for j in range(ntypes):
        for k in range(j,ntypes):
            joff = ioffmat[j+1][k+1]
            if(joff != 0):
                histnorm[i][joff] = fact/(types[j+1]*types[k+1]*xpt*xpt)
                # print(i,j,k,joff,fact,types[j+1],types[k+1],xpt)
                if(j==k) :
                    histnorm[i][joff] *= 2.  # diagonal terms need a factor of 2

histnorm[-1,1:] *= 2 # last value only is 1/2 volume so needs factor of 2 
# create header for file
hdr = ""
for i in range(len(sys.argv)):
    hdr += " " + sys.argv[i]
hdr += "\n nbins="+ str(nbins) + ", min="+str(rmin) + ", max=" + str(rmax) + ", dr=" + str(dx) + "\n Types:" + str(ntypes)
for k in sorted(types): # iterate over key, value pairs 
    hdr += " " + str(k) + ":" + str(types[k])
hdr += "\n r "
for p in pairs: # add pairs to header
    hdr += p+" "
#

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
    print(f'Config: {nconfig} [{box[0]:.3f} {box[1]:.3f} {box[2]:.3f}]',end='',flush=True)
    vol = box[0]*box[1]*box[2]
    svol += vol
    svol2 += vol*vol
    
    for it in range(ntypes):
        if(ioffmat[it+1][0] == -1):
            ipos = pos[ngbrs[it]]  # get itypes
            for jt in range(ntypes):
                if(ioffmat[0][jt+1] == -1 and jt>=it):
                    print(" ",it+1,"-",jt+1," ",sep='',end='',flush=True)
                    jbin = ioffmat[it+1][jt+1]
                    if(it==jt): # same type                        
                        for i in range(len(ipos)-1):
                            distbin(ipos[i+1:],ipos[i],jbin)
                    else: # different type
                        jpos = pos[ngbrs[jt]] # get jtypes
                        for i in range(len(ipos)):
                            distbin(jpos,ipos[i],jbin)
    print()
                        
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
stdvol = math.sqrt(abs(svol2/nconfig - avol*avol)) # abs for roundoff err close to zero
print("# configs read in:",nconfig," <vol> =",avol, "stdvol = ",stdvol)
hdr += "\n <vol> " + str(avol) + " sqrt(<(x-<x>)^2>) = " + str(stdvol)
histnorm[:,1:] *= avol/nconfig
savehist(outfile,hist,nconfig,hdr,histnorm)

