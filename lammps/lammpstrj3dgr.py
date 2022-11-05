#!/usr/bin/python3
# calculate g(dx,dy,dz) of all types from lammps traj from water

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
    data = f.readline().split() # xlo xhi
    box[0] = float(data[1])-float(data[0])
    data = f.readline().split() # ylo yhi
    box[1] = float(data[1])-float(data[0])
    data = f.readline().split() # zlo zhi
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

#    atyped = {"O":1,"H":2,"Na":3,"Cl":4,"He":5}

    for i in range(natoms): # readover header
        ntype = int(f.readline().split()[1]) # get atom type
#        ttype = f.readline().split()[0]
#        ntype = atyped[ttype]
        # print(ttype,ntype)
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
#    atyped = {"O":1,"H":2,"Na":3,"Cl":4,"He":5}
    for i in range(natoms): # readover  header
        words = f.readline().split()
        indx = int(words[0])-1 # atom index
        itype = int(words[1]) # get atom type
        atypes[indx] = itype
#        atypes[i] = atyped[words[0]]
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

    boxlo = numpy.zeros(3)
    data = f.readline().split() # xlo xhi
    boxlo[0] = float(data[0])
    box[0] = float(data[1])-float(data[0])

    data = f.readline().split() # ylo yhi
    boxlo[1] = float(data[0])
    box[1] = float(data[1])-float(data[0])

    data = f.readline().split() # zlo zhi
    boxlo[2] = float(data[0])
    box[2] = float(data[1])-float(data[0])
    
    line = f.readline() # read "ITEM: ATOMS id type xs ys zs" (scaled coord)
    if (line.rstrip() == "ITEM: ATOMS id type xs ys zs"):
        ibox = 1
    elif (line.rstrip() == "ITEM: ATOMS id type x y z ix iy iz"): # unscalled coord
        ibox = 2
    elif (line.rstrip() == "ITEM: ATOMS element x y z"): # unscalled w/ element
        ibox = 3
    else:
        print("ERROR! 9th line not \n\"ITEM: ATOMS id type xs ys zs\" or ")
        print("\"ITEM: ATOMS id type x y z ix iy iz\" or ")
        print("\"ITEM: ATOMS element x y z\"")
        exit(1)
        
    # loop over atoms
    for i in range(natoms):
        data = f.readline().split()
        if(ibox==3):
            num = i 
            pos[num] = [float(data[1]),float(data[2]),float(data[3])]
        else:
            num = int(data[0])-1 # lammps id's start with 1
            pos[num] = [float(data[2]),float(data[3]),float(data[4])]
            if(ibox==1):
                pos[num] *= box
            if(atypes[num] != int(data[1])):
                print("ERROR! atom type changed on",num+1," Config:",config)
                exit(1)
    pos -= boxlo
    pos -= numpy.floor(pos/box)*box # periodic boundary
    return 1 # read in configuration
# end readconfig
#------------------------------------------head
def ngbr(box,pos2,icell,ncell,cellatms,maxngbr): # get subcells from positions for neighbor calculations

    ntype2 = len(pos2)
    cellatms[:,:,:,0] = 0  # zero ngbr list
    for i in range(ntype2):
        ic = (pos2[i]/box*ncell).astype(int)  # cell index
        if((ic>=ncell).any() or (ic<0).any()):
            print("ERROR in ngbr",ic,i,pos[i],box,ncell)  # are we out of bounds?
            exit(1)
        icell[i] = ic # set cell id
        cellatms[ic[0]][ic[1]][ic[2]][0] += 1  # increase cell count by 1
        nl = cellatms[ic[0]][ic[1]][ic[2]][0] # get number of neighbors so far
        if(nl==maxngbr):
            maxngbr *= 2
            print("Warning: Max number of atoms/cell exceeded, increasing padding",nl,maxngbr,ic)
            exit()
            #cellatms = numpy.resize(cellatms,(ncell[0],ncell[1],ncell[2],maxngbr))
        cellatms[ic[0]][ic[1]][ic[2]][nl] = i  # cell neighbor index into cell list
#------------------------------------------
def getcellngbr(ncell):  # get cell indicies of ngbrs
    
    cellngbr = numpy.zeros((ncell[0],ncell[1],ncell[2],26,3),dtype=int)

    for l in range(ncell[0]):
        for m in range(ncell[1]):
            for n in range(ncell[2]):
                ingbr = 0
                for i in [-1,0,1]:
                    ii = l+i
                    if (ii<0):
                        ii = ncell[0]-1
                    if(ii==ncell[0]):
                        ii = 0
                    for j in [-1,0,1]:
                        jj = m+j
                        if (jj<0):
                            jj = ncell[1]-1
                        if(jj==ncell[1]):
                            jj = 0
                        for k in [-1,0,1]:
                            kk = k+n
                            if (kk<0):
                                kk  = ncell[2]-1
                            if(kk==ncell[2]):
                                kk = 0
                            if( not (i==0 and j==0 and k==0)): # don't include self
                                # print(ingbr,l,m,n,i,j,k,ii,jj,kk) 
                                cellngbr[l][m][n][ingbr] = [ii,jj,kk]
                                ingbr += 1
    return cellngbr
#------------------------------------------
def  getngbrindx(ival,idcell,ncell,cellatms,cellngbr,indl): # get the ngbr indicies
    maxl = len(indl)
    nl = cellatms[idcell[0]][idcell[1]][idcell[2]][0]
    slist = cellatms[idcell[0]][idcell[1]][idcell[2]][1:nl+1]
    if(ival == -1):  # different types will not have self
        tngbr = nl
    else:
        slist = slist[slist != ival] # remove self if there
        tngbr = len(slist)
        if(nl != tngbr+1):
            print("Error: More than 1 self atoms?",nl)
            exit(1)
    indl[0:tngbr] = slist
    for i in range(26):
        ingbr = cellngbr[idcell[0]][idcell[1]][idcell[2]][i]
        ii = ingbr[0]
        jj = ingbr[1]
        kk = ingbr[2]
        nl = cellatms[ii][jj][kk][0]
        if(tngbr+nl > maxl):
            print("Warning: expanding maxl",maxl,(tngbr+nl)*2,nl,tngbr)
            maxl = (tngbr+nl)*2
            indl = numpy.resize(indl,maxl)
        # print(ingbr,i,tngbr,nl,maxl,indl[tngbr:tngbr+nl],cellatms[ii][jj][kk])
        indl[tngbr:tngbr+nl] = cellatms[ii][jj][kk][1:nl+1]
        tngbr += nl 
    #print(asize,len(indl),indl)
    return tngbr

#------------------------------------------
def getvectwater(ii,pos,box,v):  # assumes water is O,H,H in index
    r1 = pos[ii+1]-pos[ii] # H1-OW
    r2 = pos[ii+2]-pos[ii] # H2-OW
    r1 -= numpy.floor(r1/box+.5)*box # periodic boundary
    r2 -= numpy.floor(r2/box+.5)*box # periodic boundary
    # get vectors that define molecular orientation
    v[0] = r2-r1 # H2-H1
    v[0] /= numpy.linalg.norm(v[0])
    v[2] = numpy.cross(v[0],(r1+r2))
    v[2] /= numpy.linalg.norm(v[2])
    v[1] = numpy.cross(v[2],v[0])
    v[1] /= numpy.linalg.norm(v[1])
#------------------------------------------
def getvect3atms(pos1,pos2,pos3,box,v):  # pos1 is origin
    r1 = pos2-pos1 # z vector (from pos2 and pos1)
    r2 = pos3-pos1 # x-y plane (from 
    r1 -= numpy.floor(r1/box+.5)*box # periodic boundary
    r2 -= numpy.floor(r2/box+.5)*box # periodic boundary
    # get vectors that define molecular orientation
    v[2] = r1 # z axis
    v[2] /= numpy.linalg.norm(v[2])
    v[1] = numpy.cross(r2,v[2]) # y axis
    v[1] /= numpy.linalg.norm(v[1])
    v[0] = numpy.cross(v[1],v[2])
    v[0] /= numpy.linalg.norm(v[0]) # x-axis
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
    ptn = pt[numpy.all(pt>(bins[0]-bins[2]),axis=1),:] # pull out only ones greater then -bin
    ptn = ptn[numpy.all(ptn<(bins[1]+bins[2]),axis=1),:] # pull out only ones less then bin
    ptn = (ptn-bins[0])/bins[2] # get bin numbers
    nr = numpy.floor(ptn).astype(int) # get integer bin numbers
    ptn -= nr # get fraction of bin
    np2=np-1
    for i in range(len(nr)): # loop over vectors, total sum = 12
        nri = nr[i]
        ptt = ptn[i]

        if(nri[0] > -1):
            if (nri[1] > -1):
                if (nri[2] > -1): # [0 0 0]
                    gr3d[nri[0]][nri[1]][nri[2]] += (1.-ptt[0])+(1.-ptt[1])+(1.-ptt[2])
                if(nri[2] < np2[2]): # [0 0 1]
                    gr3d[nri[0]][nri[1]][nri[2]+1] += (1.-ptt[0])+(1.-ptt[1])+(ptt[2])
            if(nri[1] <np2[1]):
                if(nri[2] > -1): #[0 1 0]
                    gr3d[nri[0]][nri[1]+1][nri[2]] += (1.-ptt[0])+(ptt[1])+(1.-ptt[2])
                if(nri[2] <np2[2]): #[0 1 1]
                    gr3d[nri[0]][nri[1]+1][nri[2]+1] += (1.-ptt[0])+(ptt[1])+(ptt[2])

        if (nri[0] < np2[0]):
            if (nri[1] > -1):
                if (nri[2] > -1): # [1 0 0]
                    gr3d[nri[0]+1][nri[1]][nri[2]] += (ptt[0])+(1.-ptt[1])+(1.-ptt[2])
                if (nri[2]<np2[2]): # [1 0 1]
                    gr3d[nri[0]+1][nri[1]][nri[2]+1] += (ptt[0])+(1.-ptt[1])+(ptt[2])
            if(nri[1]<np2[1]):
                if(nri[2] > -1): # [1 1 0]
                    gr3d[nri[0]+1][nri[1]+1][nri[2]] += (ptt[0])+(ptt[1])+(1.-ptt[2])
                if(nri[2] < np2[2]): # [1 1 1]
                    gr3d[nri[0]+1][nri[1]+1][nri[2]+1] += (ptt[0])+(ptt[1])+(ptt[2])
#-----------------------------------------------------------
def writedxfile(ofile,np,bins,gr3d,fact):
    #open file to write
    # print("Creating file :",ofile)
    f = open(ofile,"w")
    #print header
    hdr = "# opendx file for testing with 3dgr\n"
    f.write(hdr)
    hdr = "object 1 class gridpositions counts "+str(np[0])+" "+str(np[1])+" "+str(np[2])+"\n"
    f.write(hdr)
    hdr = "origin " + str(bins[0]) + " " + str(bins[0]) + " " + str(bins[0]) +"\n"
    f.write(hdr)
    hdr = "delta "+str(bins[2])+" 0 0\ndelta 0 "+str(bins[2])+" 0\ndelta 0 0 "+str(bins[2])+"\n"
    f.write(hdr)
    hdr = "object 2 class gridconnections counts "+str(np[0])+" "+str(np[1])+" "+str(np[2])+"\n"
    f.write(hdr)
    hdr = "object 3 class array type double rank 0 items "+str(np[0]*np[1]*np[2])+" data follows\n"
    f.write(hdr)
    #write data
    for i in range(np[0]):
        for j in range(np[1]):
            for k in range(np[2]): 
                #                if(i==np[0]-1 or j == np[1]-1 or k == np[1]-1 or i==0 or j==0 or k==0):
                #                    f.write(str(gr3d[i][j][k]*fact*2.)+"\n") # scale first and last points by 2 (as they are on the ends...
                #                else:
                #                    f.write(str(gr3d[i][j][k]*fact)+"\n")
                f.write(str(gr3d[i][j][k]*fact)+"\n")
                    
    f.write("\nobject density class field\n")
    print("Wrote file",ofile)

#------------------------------------------------------------
# read command line arg using argparser
parser = argparse.ArgumentParser(description="Read in Lammps trajectory file and calculate 3d g(x,y,z) for types relative to water")
parser.add_argument('input', help='input lammpstrj file')
parser.add_argument('output', help='output data file name')
parser.add_argument('-nbins',type=int,dest='nbins',default=20,help="Number of bins in each dimenstion")
parser.add_argument('-time',type=int,dest='time',default=10,help="Report and write progress every so many seconds")
parser.add_argument('-max',type=float,dest='rmax',default=10,help="Max distancs in g(r)")
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
rmin = -rmax
type1 = args.type1
type2 = args.type2

#Grids
np = numpy.array([nbins,nbins,nbins])
bins = numpy.zeros(3)
bins[0] = -rmax
bins[1] = rmax
bins[2] = 2.*rmax/(nbins-1) # should this be nbins-1?

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

indx1 = numpy.where(atypes == type1)[0] # type1 index
indx2 = numpy.where(atypes == type2)[0] # type2 index
ntype1 = indx1.size
ntype2 = indx2.size
print("Number of types 1 found = ",ntype1,"Number of type 2 =",ntype2)
#print(indx1,indx2)

# ngbr lists
icell = numpy.zeros((natoms,3),dtype=int)
ncell = numpy.array(box/rmax,dtype=int)
ncells = numpy.product(ncell)
cellngbr = getcellngbr(ncell)
maxngbr = int(ntype2/ncells*1.5+10)*2 # average density with a little padding
cellatms = numpy.zeros((ncell[0],ncell[1],ncell[2],maxngbr),dtype=int)
indl = numpy.zeros(maxngbr*27,dtype=int)
print("Neighbor Cell:",ncell)

v = numpy.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])  # default unit vectors
rpos = numpy.zeros((ntype2,3))

f = open(configname,'r')
nconfig = 0  # set counter
tnow = time.time()
ttime = tnow
svol = 0 # sum volume
svol2 = 0 # sum squared volume for stdev
print("Processing ",configname,"with",nconf,"configurations.")
while(readconfig(f,natoms,box,pos,atypes,nconfig)):
    #print(natoms,box,pos,atypes)
    nconfig += 1
    print("Processing config: ",nconfig,box)
    svol += box[0]*box[1]*box[2]
    svol2 += box[0]*box[1]*box[2]*box[0]*box[1]*box[2]

    pos2 = numpy.take(pos,indx2,axis=0)      # get type 2 positions
    ngbr(box,pos2,icell,ncell,cellatms,maxngbr) # sort type 2 pos into cells
    ival = -1
    
    for i in range(ntype1): # loop over all type 1 (Water Oxygens)
        ii = indx1[i]
        getvect3atms(pos[ii],pos[ii+1],pos[ii+2],box,v)
        #getvectwater(ii,pos,box,v) # get verticies of type1 molecule
        
        ic = (pos[ii]/box*ncell).astype(int)  # ii cell index of type1
        #writevmd(ii,pos,v)
        if(type1 == type2): # get self index
            ivals = numpy.where(indx2==ii)
            if(len(ivals) !=1):
                print("Error, found two of the same indices?",ii,indx2,ivals)
                exit(1)
            ival = ivals[0]
            # print(ii,pos[ii],ival,pos2[ival])
        ingbr = getngbrindx(ival,ic,ncell,cellatms,cellngbr,indl) #ngbr indxs of cell from pos2
        posn = numpy.take(pos2,indl[0:ingbr],axis=0) # position of neighbors
        rpos = posn-pos[ii] # distance between O and neighbors
        rpos -= numpy.floor(rpos/box+.5)*box # periodic boundary
        rpos = numpy.inner(rpos,v) # distance along each vecticies
        # print(ii,ic,pos[ii],posn,rpos)
        bingr3d(rpos,np,bins,gr3d) # create density plot

        tnowi = time.time()
        if(itime < tnowi-tnow):
            runtime = tnowi-ttime   # total time run
            fracO = i/ntype2          # fraction of atoms processes in this config
            perdone = (fracO+nconfig-1)/nconf  # percent done of total
            etime = runtime/perdone # estimated total time
            avol = svol/nconfig     # average volume
            print('ntype1 = {}/{}, configs = {}/{} ~ {:2.1%}, time={:g}/{:g}  <vol> = {}'.format(i,ntype1,nconfig,nconf,perdone,runtime,etime,avol))
            lgrid = (bins[2]) # length of grid box
            
            fact = 12*ntype1*ntype2*nconfig*bins[2]*bins[2]*bins[2]/(svol/nconfig)  # normalize...(12 is for grid) 
            writedxfile(outfile,np,bins,gr3d,1./fact)  # save dx file
            #print(numpy.amin(gr3d)*fact,numpy.amax(gr3d)*fact)
            tnow = time.time() # reset time since saved
# done processing configs, report final numbers.
avol = svol/nconfig
stdvol = math.sqrt((svol2/nconfig - avol*avol))
print("# configs read in:",nconfig," <vol> =",avol, "stdvol = ",stdvol)
print("types1:",type1,"with ",ntype1,"molecules and type2",type2," with",ntype2)

# save final dx file
fact = 12*ntype1*ntype2*nconfig*bins[2]*bins[2]*bins[2]/(avol)  # normalize...
writedxfile(outfile,np,bins,gr3d,1./fact)
print("Grid min:max",numpy.amin(gr3d),":",numpy.amax(gr3d))
print("Done!")
