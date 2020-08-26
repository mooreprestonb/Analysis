#!/usr/bin/python
# sort by column 1

import sys
import argparse
import numpy as np
from operator import itemgetter

dat = np.loadtxt("salt.atoms.3")  # read in x y .... data 

#print dat
#print
dat1 = sorted(dat,key=itemgetter(0))

indx = dat1[1:]

for x in dat1:
    # print x
    print('{:5d} {:5d} {:5d} {} {} {} {} {:3d} {:3d} {:3d}'.format(int(x[0]),int(x[1]),int(x[2]),x[3],x[4],x[5],x[6],int(x[7]),int(x[8]),int(x[9])))
#print


for x in lines:
    w = x.rstrip().split()
    nl += 1
    if w: # not a blank line
        # print w,natm
        if natm > -1:
            if not w[0].isdigit():
                # print w[0],"not digit!"
                break
            natm += 1
            atms.append(w)
            indx.append(int(w[0]))
        if w[0] == "Atoms": # start counting Atoms!
            natm = 0
            hl = nl
#print "natoms found = ",natm

print indx
