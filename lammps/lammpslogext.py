#!/usr/bin/python3
# extract thermo data from lammps log files

import sys

fi = open(sys.argv[1],"r") # read in log file
lines = fi.readlines() 
fi.close()

fo = open(sys.argv[2],"w") # file with therm data
nstep = 0
data = 0
for l in lines:
    words = l.split()
    if(len(words)>0 and words[0]=="Step"):
        if(nstep ==0):
            nstep = 1
            fo = open(sys.argv[2],"w") # file with therm data
        else:
            nstep += 1
            fo.close()
            fo = open(sys.argv[2]+"."+str(nstep),"w") # file with therm data
        fo.write("# " + l)
        data = 1
        continue # don't do rest of loop

    if(data == 1):
        if(words[0].isdigit()==False):
            data = 0
        else:
            fo.write(l)
