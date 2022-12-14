#!/usr/bin/python3

import sys
import os.path
import numpy as np
import math
import MDAnalysis as mda

PRM = "tgsn166.prmtop"
TRJ = "05_Prod4_tgsn166.nc"

#PRM = "/home/hgibbs/research/moore/vpu/tgsn166.prmtop"
#TRJ = "/home/hgibbs/research/moore/vpu/04_Prod_v21l.nc"

def get_MOI(): # get moment of intertia
    # get com
    x,y,z = r.T
    tm = len(x) # all masses are the same (ie carbon)
    xcom = sum(x)/tm  
    ycom = sum(y)/tm
    zcom = sum(z)/tm
    x -= xcom # subtract center of mass (COM)
    y -= ycom
    z -= zcom
    Ixx = np.sum(y**2 + z**2) # get MOI tensor
    Iyy = np.sum(x**2 + z**2)
    Izz = np.sum(x**2 + y**2)
    Ixy = -np.sum(x * y)
    Iyz = -np.sum(y * z)
    Ixz = -np.sum(x * z)
    I = np.array([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])
    return I

if(not os.path.isfile(PRM)):
    print(PRM," does not exist!")
    exit(1)
if(not os.path.isfile(TRJ)):
    print('\"',TRJ,"\" does not exist!")
    exit(2)
    
u = mda.Universe(PRM, TRJ) # define universe
#print(mda.__version__)
#print(u)
rad2deg = 180/math.pi
# th = acos(a.b/|a||b|)
residstart=5
residend=35

print("#",sys.argv,PRM,TRJ)
step = 0

for ts in u.trajectory:
    step += 1 # increment traj step
    CA = u.select_atoms("protein and name CA") # get alpha carbons
    r = CA.positions[residstart:residend] # get positions
    v = r[-1]-r[0]  # last - first pos (end to end)
    lv = math.sqrt(np.sum(v**2))
    th1 = math.acos(v[2]/lv)*rad2deg
    # use Moment of inertia instead of end to end
    I = get_MOI()
    di, ev = np.linalg.eig(I) # principle compents (eigenvectors and eigenvalues)
    si = np.argsort(di) # sorted index
    di = di[si] # sort eigval
    ev = ev[:,si] # sorted eigvec
    #print(step,di[0],ev[:,0])
    th = math.acos(ev[2][0])*rad2deg # get angle
    if(th>90): # make sure angle is between 0 and 90
        th = 180-th
    print(step,di[0],th,th1) # output result
