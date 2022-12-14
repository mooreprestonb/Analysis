#!/usr/bin/python3

import sys
import os.path
import numpy as np
import math
import MDAnalysis as mda
from MDAnalysis.analysis import helix_analysis as hel

PRM = "tgsn166.prmtop"
TRJ = "05_Prod4_tgsn166.nc"
#TRJ = "frame1.nc"
TRJOUT = "Trantr.nc"
ANGOUT = "Angleres.dat"

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
print(mda.__version__)
print(u)

print("#",sys.argv,PRM,TRJ)
step = 0

aa = u.select_atoms("all")
CA = u.select_atoms("protein and name CA") # get alpha carbons

rad2deg = 180/math.pi
# th = acos(a.b/|a||b|)
residstart=5-1  # offset is 1
residend=35 # offset is but want to be inclusive
rlen = residend-residstart
CAH=CA[residstart:residend] # CA atoms in helix to use for computations

#print("Running helix")
#h = hel.HELANAL(CAH,ref_axis=[0,0,1])
#print(h)
#print(h.summary['global_tilts'])  # this didn't work...

f = open(ANGOUT,"w")
anc = np.zeros(residend-residstart)

with mda.Writer(TRJOUT,aa.n_atoms) as W:
    for ts in u.trajectory:
        step += 1 # increment traj step
        com = CAH.center_of_mass()  # center COM to origin
        aa.translate(-com)
        r = CAH.positions # get positions after translation
        v = r[-1]-r[0]  # last - first pos (end to end)
        lv = math.sqrt(np.sum(v**2)) # length of vecotor 
        th1 = math.acos(v[2]/lv)*rad2deg # degrees

        # use Moment of inertia instead of end to end
        I = get_MOI()
        di, ev = np.linalg.eig(I) # principle compents (evect and eigenvalues)
        si = np.argsort(di) # sorted index
        di = di[si] # sort eigval
        ev = ev[:,si] # sorted eigvec (already normalized)
        pv = ev[:,0] # get principle vector of MOI
        if (pv[2] < 0): # make z positive
            pv = -1.*pv  # e.g. have principle compent pointing up.
            #print(pv)
        th = math.acos(pv[2])*rad2deg # get angle
        pvl = math.sqrt(pv[0]**2+pv[1]**2) # length of chosen vector in xy-planae
        Rth = math.acos(pv[1]/pvl)*rad2deg # rotation angle in degrees to yz plane
        aa.rotateby(Rth,axis=[0,0,1],point=[0,0,0]) #rotate around z by specified degrees
        print(step,di[0],th,th1,Rth)
        W.write(aa)

        r = CAH.positions # get positions after rotation
        # use Moment of inertia again after rotation
        I = get_MOI()
        di, ev = np.linalg.eig(I) # principle compents (evect and eigenvalues)
        si = np.argsort(di) # sorted index
        di = di[si] # sort eigval
        ev = ev[:,si] # sorted eigvec (already normalized)
        pv = ev[:,0] # get principle vector of MOI
        pv /= math.sqrt(np.sum(pv**2))  # normalize (should already be)
        if (pv[2] < 0): # make z positive
            pv = -1.*pv  # e.g. have principle compent pointing up.
            #print(pv)
        
        # angle of Calpha from pv which goes through origin
        out = str(step) + " " + str(th) + " " 
        for ir in range(rlen):  
            pr = np.dot(r[ir],pv)*pv  # projection along axis to closest c-alpha
            vp = r[ir]-pr # vector from axis to c-alpha point at 90 deg.
            vpl = math.sqrt(np.sum(vp**2))  # length of vector
            v1 = np.cross(pv,[0,0,1]) # define coordinate system
            vpz = np.cross(pv,v1) # define vector pointing up, looking up vector
            vpzl = math.sqrt(np.sum(vpz**2))  # length of vector
            sng = np.sign(np.dot(vp,v1)) # positive or negative from y axis
            athc = np.dot(vp,vpz)/(vpl*vpzl)
            athc = max(athc,-1)  # had rounding error... "out of bounds"
            athc = min(athc,1)
            thc = sng*math.acos(athc)
            # to print the coordinate looking up the moi vector
            #print(ir,thc*rad2deg,vpl*math.sin(thc),vpl*math.cos(thc))
            if(step > 1):  # unwrap angle by tau (2 pi)
                if(abs(anc[ir]-thc)> math.pi/2.0):
                    if(thc<0):
                        thc += math.tau
                    else:
                        thc -= math.tau
            anc[ir] = thc
            out += str(thc*rad2deg) + " " 
        f.write(out+"\n")
print("Done, configs = ",step)
