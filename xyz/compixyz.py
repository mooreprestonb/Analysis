#!/usr/bin/python3

# read in xyz and center each pi atom
import sys
import numpy

def readconfxyz(fp,n,apos,header):
    line = fp.readline() # read header or EOF
    if (line==""): # end of file! EOF!
        return 0
    natm = int(line.split()[0])
    if(natm != n):
        print("Error! num atoms don't match",natm,n)
    natm = int(line.split()[0])
    header[0] = line
    header[1] = fp.readline() # read comment line

    for i in range(n):
        line = fp.readline().split()
        apos[i][0] = float(line[1])
        apos[i][1] = float(line[2])
        apos[i][2] = float(line[3])
        
    return 1

nbeads = 8
#filename = "test2.xyz"
filename = sys.argv[1]
if(len(sys.argv)==3):
    nbeads = int(sys.argv[2])
    
fp = open(filename,'r')

line = fp.readline()
natoms = int(line.split()[0])
if(natoms%nbeads):
    print("Error: natoms in not multiple of nbeads")
    exit(1)
npath = int(natoms/nbeads)
#lines = [None]*natoms
header = [None]*2
apos = numpy.empty((natoms,3),dtype=float)
momn = numpy.zeros((2,3),dtype=float)

#print(natoms,nbeads,npath)
fp.seek(0,0) # rewind

nconf = 0
while(readconfxyz(fp,natoms,apos,header)):
    nconf += 1
    print(header[0],header[1],end='')
    #print(apos)
    for j in range(npath):
        npos = apos[j::npath] # slice out individual paths
        avg = numpy.mean(npos,axis=0)
        npos -= avg
#        print(numpy.mean(npos,axis=0))
        for i in range(nbeads):
            print("Ar",npos[i][0],npos[i][1],npos[i][2])
        momn[0] = numpy.sum(npos,axis=0) # should be zero 
        momn[1] = numpy.sum(npos**2,axis=0) # variance

print(nconf)
print(momn)
