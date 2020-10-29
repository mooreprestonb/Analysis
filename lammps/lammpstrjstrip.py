#!/usr/bin/python3
# code to read lammps trj files and strip out every so often

import sys

infile = "0.05_test.lammpstrj"
outfile = "0.05_strip.lammpstrj"
stride = 10 

def readconfiglines(): #subroutine to read configuration lines
    lines[0] = lammpsfp.readline()
    if (lines[0].strip() != "ITEM: TIMESTEP"):
        return 0 # not read in... End of file...
    for i in range(1,nlines): # already read in first line
        lines[i] = lammpsfp.readline()
    return 1

print(sys.argv[0],infile,outfile,stride)

#readheader

lammpsfp = open(infile,"r")
head = [ next(lammpsfp) for i in range(9) ] # readin first 9 lines
lammpsfp.seek(0) # rewind file

# error check header
if(head[0].strip() != "ITEM: TIMESTEP"):
    print("Error: not a lammps traj file?")
    exit(1)
        
if(head[2].strip() != "ITEM: NUMBER OF ATOMS"):
    print("Error: now atoms line found")
    exit(1)

natoms = int(head[3]) # set number of atoms
nlines = natoms+9     # number of lines per configurations
lines = [None]*nlines # allocate space

lammpsoutfp = open(outfile,"w") # open output file

nconf = 0
while(readconfiglines()):  # loop over configurations
    if(nconf%stride==0): # only write out every stride times
        for line in lines:
            lammpsoutfp.write(line)
        print("Config",nconf,"output") # report config written
    nconf += 1

lammpsoutfp.close()
lammpsfp.close()

print(nconf," Configurations read in\nDone")
