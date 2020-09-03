#!/usr/bin/python3

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

#filename = "test2.xyz"
filename = sys.argv[1]
fp = open(filename,'r')

line = fp.readline()
natoms = int(line.split()[0])
#lines = [None]*natoms
header = [None]*2
apos = numpy.empty((natoms,3),dtype=float)
#print(natoms)
fp.seek(0,0)

nconf = 0
while(readconfxyz(fp,natoms,apos,header)):
    nconf += 1
    print(header[0],header[1],end='')
    #print(apos)
    avg = numpy.mean(apos,axis=0)
    npos = apos-avg
    for i in range(natoms):
        print("Ar",npos[i][0],npos[i][1],npos[i][2])
    
#print(nconf)
