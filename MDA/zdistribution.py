#!/usr/bin/python3

import sys
import numpy as np
import math
import MDAnalysis as mda

u = mda.Universe("tgsn166.prmtop", "05_Prod4_tgsn166.nc") # define universe
#print(mda.__version__)
#print(u)
rad2deg = 180/math.pi
# th = acos(a.b/|a||b|)
residstart=5
residend=35

step = 0
for ts in u.trajectory:
    step += 1 # increment traj step
    CA = u.select_atoms("protein and name CA") # get alpha carbons
    r = CA.positions[residstart:residend] # get positions
    z = r[:,2] # grab z coordinates of C_alpha carbons
    v = r[-1]-r[0]  # last - first pos (end to end)
    lv = math.sqrt(np.sum(v**2))
    th1 = math.acos(v[2]/lv)*rad2deg
    print(step,*z,lv,th1) # output result
