#!/usr/bin/python

# python script to create g(r) from lammpstrj files 

import sys
import argparse
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
        for i in range(nbins):
                r = (i + 0.5) * dx + hmin
                factor = prefact*(r**2) 
		y = hist[i]/factor
		s = str(r)+ ' ' + str(y)+'\n'
		fo.write(s)
	# end for
	fo.close()

def splitconf(line): # split line every 8 characters
	a = []
	n = len(line) - 1
	for i in range(0,n,8):
		a.append(float(line[i:i+8]))
	return a

def readconfig(): # read lammps configuration
	global f
	global pos 
        global box
	i = 0
	pos=[]
        while i < numlines:
                line = f.readline()
                if line == '':
                        return 0
                i = i + 1
                a = splitconf(line) 
                pos.extend(a)

        line = f.readline()
	if line == '':
		print 'Error: Premeture ending on config file'
		exit(1)
        a = line.split()
        if float(a[0]) != box:
                print 'Error at box', a[0], box 
                exit(1)
        return 1
# end 

# read command line arg using argparser

parser = argparse.ArgumentParser(description=
				 "Calculate g(r) from amber configs.")
parser.add_argument('input', help='input amber config file')
parser.add_argument('output', help='output g(r) file name')
parser.add_argument('-iatoms', type=int, nargs="+", dest='ia', required=True,
		    help="i atom list")
parser.add_argument('-jatoms', type=int, nargs="+", dest='ja', required=True,
		    help="j atom list")
parser.add_argument('xbox', type=float, help="x box length")
parser.add_argument('natoms', type=int, help="number of atoms in system")
parser.add_argument('-self', default=False, action="store_true", dest="s",
		    help="Include self molecule")
parser.add_argument('-iwrite', type=int, dest='iwrite', default=50,
		    help="when to write g(r)")

# readin arguments, now translate to variables
args = parser.parse_args()

configname = args.input # 'tbac-water512-nvtprod.conf'
grname = args.output # 'tbacnitclgr.dat'
box = args.xbox # 66.492
natoms = args.natoms # 29184 
iself = args.s # do we include same molecule, default = False
nwrite = args.iwrite # how often to write out histogram

# calculate number of lines/config from number of atoms
numlines = natoms * 3 / 10
if (natoms * 3) % 10 != 0: # not exact match so add line to get remainder
	numlines = numlines + 1

print "Called as", str(sys.argv)

f = open (configname, 'r')
line = f.readline() # read header

n = 512 # number of tbac's 
v = box ** 3 # volume of box

#histogram parameters
nbins = 100
hist = [0] * nbins
hmin = 0.0
hmax = 30.0
dx = (hmax - hmin) / float(nbins)

ilist = args.ia # list of iatoms
jlist = args.ja # list of jatoms

im = len(ilist) 
jm = len(jlist) 
ri = [0]*3
rj = [0]*3

config = 1  # read initial configuration
eof = readconfig()

while eof == 1: # loop until no more configs
	print config
	for i in range(n):
                for k in range(im):
                        ii = i * 54 * 3 + (ilist[k] * 3)
			ri = pos[ii:ii+3] # slice out x,y,z of atom i
                        for j in range(n):
				if i==j and iself==False: # don't use same mol
					# print i,j,iself
					continue # go to next iteration of j
                                for l in range(jm):
                                        jj = j * 54 * 3 + (jlist[l] * 3)
                                        rj = pos[jj:jj+3] # slice j's x,y,z
					r2 = 0
					for kk in range(3) : # periodic image
						dr = rj[kk] - ri[kk]
						dr -= box*anint(dr/box)
						r2 += dr*dr
					# histogram distance
                                        histogram(sqrt(r2),hist,hmin,dx,nbins) 
                #       print r

        if( config % nwrite) == 0: # save histogram every nwrite configs
                histprint(grname,sys.argv,hist,hmin,dx,nbins,
			  config,n*im,n*jm,v)
	
	config = config + 1 # read next config
	eof = readconfig()
#end while 

config = config - 1 # reset config to last config
print 'read in ',config,' configurations'
histprint(grname,sys.argv,hist,hmin,dx,nbins,config,n*im,n*jm,v)

exit(0)
