#!/usr/bin/python

# python script to create g(r) from lammpstrj files 

import sys
import argparse
import numpy
from math import sqrt,pi,floor

def anint(x):
        return floor(x+0.5)

def histogram(r,hist,hmin,dx,nbins):
        i = int((r - hmin) / dx)
#	print i,r,hmin,nbins
        if (i < nbins) and (i >= 0):
                hist[i] += 1

# print histogram
def histprint(grname,header,hist,nmin,dx,nbins,config,ni,nj,v):
	print 'writing g(r) to ',grname,' with ',config,' configs read'
	fo = open (grname, 'w')
 	fo.write('#' + str(header) + ' ' + str(config) + '\n')
	prefact = config*ni*nj*(4.0*pi*dx)/v
        print prefact,dx,hmin,config,ni,nj,pi,dx,v
        for i in range(nbins):
                r = (i + 0.5) * dx + hmin
                factor = prefact*(r**2) 
		y = hist[i]/factor
		s = str(r)+ ' ' + str(y)+'\n'
		fo.write(s)
	# end for
	fo.close()

def readconfig(): # read lammps configuration
	global f
	global pos 
        global box
        global config
        global natoms
        global atypes
        
        line = f.readline() # read header "ITEM: TIMESTEP"
        if (line==""): # end of file!
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
        return 1
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
parser.add_argument('-self', default=False, action="store_true", dest="s",
		    help="Include self molecule")
parser.add_argument('-iwrite', type=int, dest='iwrite', default=50,
		    help="when to write g(r)")

# readin arguments, now translate to variables
args = parser.parse_args()

configname = args.input # '0.05_test.lammpstrj'
grname = args.output # 'grlammps.dat'
iself = args.s # do we include same molecule, default = False
nwrite = args.iwrite # how often to write out histogram

print "Called as", str(sys.argv)

f = open (configname, 'r')

#histogram parameters
nbins = 100
hist = [0] * nbins
hmin = 0.0
hmax = 30.0
dx = (hmax - hmin) / float(nbins)

itype = int(args.ia) # itypes
jtype = int(args.ja) # jtypes

ri = numpy.zeros(3)
rj = numpy.zeros(3)
box = numpy.zeros(3)

config = 0  # read initial configuration
eof = readconfig()

iap = numpy.array([],dtype=int)
jap = numpy.array([],dtype=int)

# parse out atom types indices we want
for i in range(natoms):
        if(atypes[i] == itype):
                iap = numpy.append(iap,i*3)
        if(atypes[i] == jtype):
                jap = numpy.append(jap,[i*3,i*3+1,i*3+2])

print "Found ",iap.size,"itypes and",jap.size/3,"jtypes.",eof

while eof == 1: # loop until no more configs
	print "Config: ",config,box
	for i in range(iap.size):
                ii = iap[i]
		ri = pos[ii:ii+3] # slice out x,y,z of atom i
                rj = numpy.reshape(numpy.take(pos,jap),(-1,3))-ri # distances
                dr = (rj-numpy.floor(rj+.5))*box # peroidic boundaries 
                r2 = numpy.zeros(dr.shape[0])
                for i in range(r2.size):
		        r2[i] = numpy.dot(dr[i],dr[i])
			# histogram distance
                        if (r2[i]>hmin):
                                histogram(sqrt(r2[i]),hist,hmin,dx,nbins) 
                #       print r

        if(( config % nwrite) == 0): # save histogram every nwrite configs
                v = box[0]*box[1]*box[2]
                histprint(grname,sys.argv,hist,hmin,dx,nbins,
			  config,iap.size,jap.size/3.,v)
	
	eof = readconfig() # read next config
#end while 

v = box[0]*box[1]*box[2]
print 'Read in ',config,' configurations'
histprint(grname,sys.argv,hist,hmin,dx,nbins,config,iap.size,jap.size/3.,v)

exit(0)
