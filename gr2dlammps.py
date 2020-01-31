#!/usr/bin/python

# python script to create g(r) from lammpstrj files 

import sys
import argparse
import numpy
from math import pi,sqrt

def histogramr(r,hist,hmin,dx,nbins):
        i = int((r - hmin) / dx)
        #        if (i < (nbins-1)) and (i >= 0):
        #                hist[i] += 1
        fr = (r - i*dx+hmin)/dx # fraction of bin width
        hist[i] += 1.0-fr
        hist[i+1] += fr

# print histogram
def histprint(grname,header,hist,nmin,dx,nbins,config,ni,nj,v):
	print 'Writing g(r) to ',grname,' with ',config,' configs read'
	fo = open (grname, 'w')
 	fo.write('#' + str(header) + ' ' + str(config) + '\n')
	prefact = config*ni*nj*(4.0*pi*dx)/v
        #print prefact,dx,hmin,config,ni,nj,pi,dx,v
        for i in range(nbins-1):
                #r = (i + 0.5) * dx + hmin # offset r
                r = i*dx + hmin 
                if (r==0):
                        y = 0
                else:
		        y = hist[i]/(prefact*r*r)
		s = str(r)+ ' ' + str(y)+'\n'
		fo.write(s)
	# end for
        
	y = 2.0*hist[-1]/(prefact*hmax*hmax) # last bin gets 1/2 density.
	s = str(hmax)+ ' ' + str(y)+'\n'
	fo.write(s)
	fo.close()

def readconfig(): # read lammps configuration
	global f
	global pos 
        global box
        global config
        global natoms
        global atypes
        
        line = f.readline() # read header or EOF
        if (line==""): # end of file! EOF!
                return 0
        if (line.rstrip() != "ITEM: TIMESTEP"):
                print "ERROR! first line of config is not \"ITEM: TIMESTEP\""
                exit(1)
        line = f.readline() # read timestep
        ts = int(line)

        line = f.readline() # read "ITEM: NUMBER OF ATOMS"
        if (line.rstrip() != "ITEM: NUMBER OF ATOMS"):
                print "ERROR! 3rd line of config is not \"ITEM: NUMBER OF ATOMS\""
                exit(1)
        line = f.readline() # natoms
        natms = int(line)
        if(config==0):
                natoms = natms
                pos = numpy.zeros(natoms*3)
                atypes = numpy.zeros(natoms,dtype=int)
        else:
                if(natoms != natms):
                        print "ERROR! natoms can not change!",natoms,natms
                        exit(1)
        
        line = f.readline() # read "ITEM: BOX BOUNDS pp pp pp"
        if (line.rstrip() != "ITEM: BOX BOUNDS pp pp pp"):
                print "ERROR! 3rd line of config is not \"ITEM: BOX BOUNDS pp pp pp\""
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
                print "ERROR! 3rd line of config is not \"ITEM: ATOMS id type xs ys zs\""
                exit(1)
        
        # loop over atoms
        for i in range(natoms):
                line = f.readline()
                data = line.split()
                num = int(data[0])-1 # lammps number start with 1
                pos[num*3+0] = float(data[2])
                pos[num*3+1] = float(data[3])
                pos[num*3+2] = float(data[4])
                if (config==0):
                        atypes[num] = int(data[1])
                else:
                        if(atypes[num] != int(data[1])):
                                print "ERROR! atom type change on",num
        config += 1
        return 1 # read in configuration
# end 

# read command line arg using argparser

parser = argparse.ArgumentParser(description=
				 "Calculate g(r) from Lammps trajectory.")
parser.add_argument('input', help='input lammpstrj file')
parser.add_argument('output', help='output g(r) file name')
parser.add_argument('-iatoms', type=int, dest='ia', required=True,
		    help="i atom type")
parser.add_argument('-jatoms', type=int, dest='ja', required=True,
		    help="j atom type")
parser.add_argument('-iwrite', type=int, dest='iwrite', default=50,
		    help="when to write g(r)")
parser.add_argument('-hmin', type=float, dest='hmin', default=0,
		    help="Min in histogram")
parser.add_argument('-hmax', type=float, dest='hmax', default=-1,
		    help="Max in histogram")
parser.add_argument('-nbins', type=int, dest='nbins', default=100,
		    help="Number of bins in histogram")

# readin arguments, now translate to variables
args = parser.parse_args()

configname = args.input # '0.05_test.lammpstrj'
grname = args.output # 'grlammps.dat'
nwrite = args.iwrite # how often to write out histogram

print "Called as", str(sys.argv),"\nUsing arguments",args

f = open (configname, 'r')

#histogram parameters
nbins = args.nbins
hist = numpy.zeros(nbins)
hmin = args.hmin

itype = int(args.ia) # itypes
jtype = int(args.ja) # jtypes

box = numpy.zeros(3)
abox = numpy.zeros(3)

config = 0  # read initial configuration
eof = readconfig()

if(args.hmax == -1): # histogram dimensions
        hmax = numpy.amin(box)/2.0
else:
        hmax = args.hmax
dx = (hmax - hmin)/float(nbins-1)

iap = numpy.array([],dtype=int)
jap = numpy.array([],dtype=int)

# parse out atom types indices we want
for i in range(natoms):
        if(atypes[i] == itype):
                iap = numpy.append(iap,i*3)
        if(atypes[i] == jtype):
                jap = numpy.append(jap,[i*3,i*3+1,i*3+2])

ni = iap.size
nj = jap.size/3
print "Found ",ni,"itypes and",nj,"jtypes."
r2 = numpy.zeros(nj)
ri = numpy.zeros(3)
rj = numpy.zeros(3)

while eof == 1: # loop until no more configs
	print "Config: ",config,box
        abox += box
	for i in range(ni): # loop over iatoms
                ii = iap[i]
		ri = pos[ii:ii+3] # slice out x,y,z of atom i
                # get rj coords into a nx3 vector and subtract to get dist
                rj = numpy.reshape(numpy.take(pos,jap),(nj,3))-ri # distances
                dr = (rj-numpy.floor(rj+.5))*box # peroidic boundaries in 3d 
                for i in range(nj):
		        r = sqrt(numpy.dot(dr[i],dr[i])) # get distant
                        if (r>hmin and r<hmax): 	# histogram distance
                                histogramr(r,hist,hmin,dx,nbins)

        if(( config % nwrite) == 0): # save histogram every nwrite configs
                v = (abox[0]/config)*(abox[1]/config)*(abox[2]/config)
                histprint(grname,sys.argv,hist,hmin,dx,nbins,
			  config,ni,nj,v)
	
	eof = readconfig() # read next config
#end while 

abox /= config
v = abox[0]*abox[1]*abox[2]
print 'Read in ',config,"configurations with average box",abox
histprint(grname,sys.argv,hist,hmin,dx,nbins,config,ni,nj,v)

exit(0)
