#!/usr/bin/python3
# calculate g(dx,dy,dz) of all types from lammps traj from water

import sys
import argparse
import numpy
import math
import time
import os.path
import getlammpstrj
import networkx as nx
import matplotlib.pyplot as plt

#------------------------------------------------------------
# read command line arg using argparser
parser = argparse.ArgumentParser(description="Read in Lammps trajectory file and calculate neighbors and clusters")
parser.add_argument('input', help='input lammpstrj file')
parser.add_argument('-time',type=int,dest='time',default=10,help="Report and write progress every so many seconds")
parser.add_argument('-max',type=float,dest='rmax',default=2.5,help="Max Neighbor distancs")
parser.add_argument('-type1', type=int, dest='type1',default=1,help="Atoms type for neighbor")
parser.add_argument('-type2', type=int, dest='type2',default=1,help="Neigbor atom type")

# read in arguments, now translate to variables
args = parser.parse_args()
print (args)
configname = args.input # 'test.lammpstrj'
itime = args.time
rmax = args.rmax
type1 = args.type1
type2 = args.type2

rmax2 = rmax*rmax

outfile = "out.out"
print("Processing",configname,"to",outfile)
natoms = getlammpstrj.getlammpsatoms(configname)
print("Number of atoms:",natoms)

types = getlammpstrj.getlammpstypes(configname,natoms)
ntypes =len(types.keys())
stypes = str(ntypes)
for k in sorted(types): # iterate over key, value pairs 
    stypes += " " + str(k) + ":" + str(types[k])
print("Number of atom types: ",stypes)

box = getlammpstrj.getlammpsbox(configname)
print("Initial box dimensions = ",box)

nconf = getlammpstrj.getlammpsnconf(configname,natoms) # count configurations

# check types?
    
# allocate arrays
pos = numpy.zeros((natoms,3))
atypes = numpy.zeros(natoms,dtype=int)
getlammpstrj.getlammpsatypes(configname,natoms,atypes)

#print(gr3d.shape,gr3d)
G = nx.Graph()

indx1 = numpy.where(atypes == type1)[0] # type1 index
indx2 = numpy.where(atypes == type2)[0] # type2 index
ntype1 = indx1.size
ntype2 = indx2.size
print("Number of types 1 found = ",ntype1,"Number of type 2 =",ntype2)
#print(indx1,indx2)

f = open(configname,'r')
nconfig = 0  # set counter
tnow = time.time()
ttime = tnow
single = 0
svol = 0 # sum volume
svol2 = 0 # sum squared volume for stdev
specs = {}

print("Processing ",configname,"with",nconf,"configurations.")
while(getlammpstrj.readconfig(f,natoms,box,pos,atypes,nconfig)):
    #print(natoms,box,pos,atypes)
    nconfig += 1
    G.clear()
    
    print("Processing config: ",nconfig,box)
    svol += box[0]*box[1]*box[2]
    svol2 += box[0]*box[1]*box[2]*box[0]*box[1]*box[2]

    pos1 = numpy.take(pos,indx1,axis=0) # get type 1 positions
    for i in range(len(pos1)):
        # print(i,indx1[i],pos[indx1[i]],pos1[i])
        G.add_node(indx1[i],pos = (pos1[i]),atype=type1)
        # G.add_nodes_from([(indx1[i],{"pos": (pos1[i])})])
    #print(G)
    pos2 = numpy.take(pos,indx2,axis=0)      # get type 2 positions
    for i in range(len(pos2)):
        # print(i,indx2[i],pos[indx2[i]])
        G.add_node(indx2[i],pos=(pos2[i]),atype=type2)

    for nodei in G.nodes():
        posi = G.nodes[nodei].get('pos')
        for nodej in G.nodes():
            if(nodei < nodej):
                posj = G.nodes[nodej].get('pos')
                rpos = posj-posi # distance between neighbors
                rpos -= numpy.floor(rpos/box+.5)*box # periodic boundary
                d2 = numpy.sum(numpy.square(rpos)) # distance
                if(d2<rmax2):
                    G.add_edge(nodei,nodej)
    #G.add_edges_from(nx.geometric_edges(G,radius=rmax))
    #print(G)

    nnodes = G.number_of_nodes() 
    while(nnodes>0):
        nodei = list(G.nodes)[0] # pick first node in list
        nodeis = list(nx.node_connected_component(G,nodei))
        keyt = []
        for ns in nodeis:
            keyt.append(G.nodes[ns].get('atype'))
        keyt.sort()
        key = str(keyt[0])
        for i in range(1,len(keyt)):
            key += "-"+str(keyt[i])
        if key in specs:
            specs[key] += 1
        else:
            specs[key] = 1
        G.remove_nodes_from(nodeis)
        nnodes = G.number_of_nodes() 

    # nx.draw(G,with_labels=True)
    # atom_types = nx.get_node_attributes(G,'atype')
    # nx.draw(G,labels=atom_types)
    # plt.show()
    # exit(1)

    tnowi = time.time()
    if(itime < tnowi-tnow):
        runtime = tnowi-ttime   # total time run
        perdone = (nconfig-1)/nconf  # percent done of total
        etime = runtime/perdone # estimated total time
        avol = svol/nconfig     # average volume
        print('configs = {}/{} ~ {:2.1%}, time={:g}/{:g}  <vol> = {}'.format(nconfig,nconf,perdone,runtime,etime,avol))
        print(specs)
        tnow = time.time() # reset time since saved

# done processing configs, report final numbers.
print(specs)
avol = svol/nconfig
stdvol = math.sqrt((svol2/nconfig - avol*avol))
print("# configs read in:",nconfig," <vol> =",avol, "stdvol = ",stdvol)
print("types1:",type1,"with",ntype1,"and type2:",type2," with",ntype2)
print("Done!")
