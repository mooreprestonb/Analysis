#!/usr/bin/python3
# routines that read in lammps trajectories files

import numpy
import subprocess

#------------------------------------------------------------
def getlammpsbox(configname): # box in lammps traj
    f = open(configname,'r')

    box = numpy.zeros(3)
    line = f.readline() # read header or EOF
    if (line==""): # end of file! EOF!
        return 0
    if (line.rstrip() != "ITEM: TIMESTEP"):
        print("ERROR! first line not \"ITEM: TIMESTEP\"")
        exit(1)
        
    line = f.readline() # read timestep
    ts = int(line)

    line = f.readline() # read "ITEM: NUMBER OF ATOMS"
    if (line.rstrip() != "ITEM: NUMBER OF ATOMS"):
        print("ERROR! 3rd line not \"ITEM: NUMBER OF ATOMS\"")
        exit(1)
    line = f.readline() # natoms
    natms = int(line)
            
    line = f.readline() # read "ITEM: BOX BOUNDS pp pp pp"
    if (line.rstrip()[:16] != "ITEM: BOX BOUNDS"):
        print("ERROR! 3rd line not \"ITEM: BOX BOUNDS pp pp pp\"")
        exit(1)
    data = f.readline().split() # xlo xhi
    box[0] = float(data[1])-float(data[0])
    data = f.readline().split() # ylo yhi
    box[1] = float(data[1])-float(data[0])
    data = f.readline().split() # zlo zhi
    box[2] = float(data[1])-float(data[0])

    f.close()
    return box

#-----find number of configs (quicker with "wc" but not as portable) ------
def getlammpsnconf(configname,natoms): 
    lines = 0
    print("Counting configures of ",configname)
    blines = subprocess.check_output(["wc","-l",configname])
    # lines = int((str(blines,'utf-8').split())[0])
    lines = int((str(blines.decode()).split())[0])
#    for line in open(configname): # above is 5* faster, but this always works
#        lines += 1
    nconf = int(lines/(natoms+9))
    print ("# of configs = ",nconf)
    return (nconf)
#------------------------------------------------------------
def getlammpstypes(configname,natoms): # find # types in lammps traj
    types={}
    f = open(configname,'r')

    atyped = {"O":1,"H":2,"Na":3,"Cl":4,"He":5}

    for i in range(9): # readover header
        line = f.readline()
    if (line.rstrip() == "ITEM: ATOMS id type xs ys zs"):
        ibox = 1
    elif (line.rstrip() == "ITEM: ATOMS id type x y z ix iy iz"): # unscalled coord
        ibox = 2
    elif (line.rstrip() == "ITEM: ATOMS element x y z"): # unscalled w/ element
        ibox = 3
    else:
        print("ERROR! 9th line not \n\"ITEM: ATOMS id type xs ys zs\" or ")
        print("\"ITEM: ATOMS id type x y z ix iy iz\" or ")
        print("\"ITEM: ATOMS element x y z\"")
        exit(1)

    for i in range(natoms): # readover header
        if(ibox == 3): # atoms identified by elements
            ttype = f.readline().split()[0]
            ntype = atyped[ttype]
            # print(ttype,ntype)
        else: # atoms identified by id and type (scalled or unscalled)
            ntype = int(f.readline().split()[1]) # get atom type
        if ntype in types: 
            types[ntype] += 1 # add to types
        else :
            types[ntype] = 1 # add types if doesn't exist 
    f.close()
    return types
#------------------------------------------------------------ 
def getlammpsatypes(configname,natoms,atypes): # find type of each atom in lammps traj first config
    f = open(configname,'r')

    atyped = {"O":1,"H":2,"Na":3,"Cl":4,"He":5}
    for i in range(9): # readover header
        line = f.readline()
    if (line.rstrip() == "ITEM: ATOMS id type xs ys zs"):
        ibox = 1
    elif (line.rstrip() == "ITEM: ATOMS id type x y z ix iy iz"): # unscalled coord
        ibox = 2
    elif (line.rstrip() == "ITEM: ATOMS element x y z"): # unscalled w/ element
        ibox = 3
    else:
        print("ERROR! 9th line not \n\"ITEM: ATOMS id type xs ys zs\" or ")
        print("\"ITEM: ATOMS id type x y z ix iy iz\" or ")
        print("\"ITEM: ATOMS element x y z\"")
        exit(1)

    for i in range(natoms): # loop over atoms
        words = f.readline().split()
        if(ibox == 3):
            atypes[i] = atyped[words[0]] # if type is elements
        else :
            indx = int(words[0])-1 # atom index
            itype = int(words[1]) # get atom type
            atypes[indx] = itype  # if type is numbers
    f.close()
#------------------------------------------------------------ 
def getlammpsatoms(configname): # find # atoms in lammps header
    f = open(configname,'r')

    line = f.readline() # read header or EOF
    if (line==""): # end of file! EOF!
        print("EOF while reading header?!?")
        return 0
    if (line.rstrip() != "ITEM: TIMESTEP"):
        print("ERROR! first line of config is not \"ITEM: TIMESTEP\"")
        exit(1)
    line = f.readline() # read timestep
    ts = int(line)

    line = f.readline() # read "ITEM: NUMBER OF ATOMS"
    if (line.rstrip() != "ITEM: NUMBER OF ATOMS"):
        print("ERROR! 3rd line of config is not \"ITEM: NUMBER OF ATOMS\"")
        exit(1)
    line = f.readline() # natoms
    natms = int(line)
    f.close()
    return(natms)
#------------------------------------------------------------
def readconfig(f,natoms,box,pos,atypes,config): # read lammps configuration

    line = f.readline() # read header or EOF
    if (line==""): # end of file! EOF!
        return 0
    if (line.rstrip() != "ITEM: TIMESTEP"):
        print("ERROR! first line not \"ITEM: TIMESTEP\"")
        exit(1)
        
    line = f.readline() # read timestep
    ts = int(line)

    line = f.readline() # read "ITEM: NUMBER OF ATOMS"
    if (line.rstrip() != "ITEM: NUMBER OF ATOMS"):
        print("ERROR! 3rd line not \"ITEM: NUMBER OF ATOMS\"")
        exit(1)
    line = f.readline() # natoms
    natms = int(line)
    if(natoms != natms):
        print("ERROR! natoms can not change!",natoms,natms," config",config)
        exit(1)
            
    line = f.readline() # read "ITEM: BOX BOUNDS pp pp pp"
    if (line.rstrip()[:16] != "ITEM: BOX BOUNDS"):
        print("ERROR! 5th line not \"ITEM: BOX BOUNDS pp pp pp\"")
        exit(1)

    boxlo = numpy.zeros(3)
    data = f.readline().split() # xlo xhi
    boxlo[0] = float(data[0])
    box[0] = float(data[1])-float(data[0])

    data = f.readline().split() # ylo yhi
    boxlo[1] = float(data[0])
    box[1] = float(data[1])-float(data[0])

    data = f.readline().split() # zlo zhi
    boxlo[2] = float(data[0])
    box[2] = float(data[1])-float(data[0])
    
    line = f.readline() # read "ITEM: ATOMS id type xs ys zs" (scaled coord)
    if (line.rstrip() == "ITEM: ATOMS id type xs ys zs"):
        ibox = 1
    elif (line.rstrip() == "ITEM: ATOMS id type x y z ix iy iz"): # unscalled coord
        ibox = 2
    elif (line.rstrip() == "ITEM: ATOMS element x y z"): # unscalled w/ element
        ibox = 3
    else:
        print("ERROR! 9th line not \n\"ITEM: ATOMS id type xs ys zs\" or ")
        print("\"ITEM: ATOMS id type x y z ix iy iz\" or ")
        print("\"ITEM: ATOMS element x y z\"")
        exit(1)
        
    # loop over atoms
    for i in range(natoms):
        data = f.readline().split()
        if(ibox==3):
            num = i 
            pos[num] = [float(data[1]),float(data[2]),float(data[3])]
        else:
            num = int(data[0])-1 # lammps id's start with 1
            pos[num] = [float(data[2]),float(data[3]),float(data[4])]
            if(ibox==1):
                pos[num] *= box
            if(atypes[num] != int(data[1])):
                print("ERROR! atom type changed on",num+1," Config:",config)
                exit(1)
    pos -= boxlo
    pos -= numpy.floor(pos/box)*box # periodic boundary
    return 1 # read in configuration
# end readconfig
