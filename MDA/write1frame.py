#!/usr/bin/python3

import sys
import os.path
import numpy as np
import math
import MDAnalysis as mda

PRM = "tgsn166.prmtop"
TRJ = "05_Prod4_tgsn166.nc"

if(not os.path.isfile(PRM)):
    print(PRM," does not exist!")
    exit(1)
if(not os.path.isfile(TRJ)):
    print('\"',TRJ,"\" does not exist!")
    exit(2)
    
u = mda.Universe(PRM, TRJ) # define universe
print(mda.__version__)
print(u)

print("#",sys.argv,PRM,TRJ)
step = 0

aa = u.select_atoms("all")
aa.write("frame1.nc")

exit(1)

