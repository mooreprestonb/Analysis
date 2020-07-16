#!/usr/bin/python3

import numpy
import re

print("Hello")

file1 = "avg-receptor-5vai-100ns.pdb"
file2 = "betavals-5vai.dat"
file3 = "new.pdb"

fin = open(file2,"r")
beta = fin.readlines()
fin.close()

fin = open(file1,"r")
pdb = fin.readlines()
fin.close()

resid = 0
fout = open(file3,"w")
for line in pdb:
    if(line.split()[0] != "ATOM"):
        fout.write(line)
    else:
        val = line.split()
        if(resid+1==int(val[5])):
            resid += 1
            bval = beta[resid].split()[1]
            #print(resid,bval)
#        nline = re.sub(r"0.00  0.00 ","0.00  %s " % bval[:4],line.rstrip()) 
        nline = re.sub(r"0.00  0.00 ","0.00  %s " % bval[:4],line) 
        fout.write(nline)

fout.close()
print("goodbye!")
