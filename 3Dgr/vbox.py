#! /usr/bin/python3
#
# Writes a Blender Voxel accepting this format:
# ----------------------------
# X Y Z
# data data data data data
# ----------------------------
# Where X Y Z are the dimentions (separated by spaces)
# data is a float between 0.0 and 1.0
#
# Data is read by Blender in this manner:
# Z * resY * resX + Y * resX + X
#
# So write your data line like this for a 3 x 3 voxel box:
# (0, 0, 0) (1, 0, 0) (2, 0, 0)
# (0, 1, 0) (1, 1, 0) (2, 1, 0)
# (0, 2, 0) (1, 2, 0) (2, 2, 0)
#
# (0, 0, 1) (1, 0, 1) (2, 0, 1)
# (0, 1, 1) (1, 1, 1) (2, 1, 1)
# (0, 2, 1) (1, 2, 1) (2, 2, 1)
#
# (0, 0, 2) (1, 0, 2) (2, 0, 2)
# (0, 1, 2) (1, 1, 2) (2, 1, 2)
# (0, 2, 2) (1, 2, 2) (2, 2, 2)
#
# Remember to write these all on 1 line in a constant row

import sys, os
from struct import *

def main():
   if (len(sys.argv) != 2):
	   print("Usage: error!",sys.argv)
	   return
   
   SOURCE = sys.argv[1]
   
   if not os.path.exists(SOURCE):
	   print("Can't find source file " + SOURCE)
	   return
   
   source = open(SOURCE, "r")
   
   destination = open('./voxels.bvox', "wb")
   
   lines = source.readlines()
   
   X, Y, Z = lines[0].split(' ')
   
   destination.write(pack('i', int(X)))
   destination.write(pack('i', int(Y)))
   destination.write(pack('i', int(Z)))
   destination.write(pack('i', int(50))) # 1 frame

   print(lines[1].split(" "))
   for value in lines[1].split(" "):
        if(value != "\n"):
                destination.write(pack('f', float(value)))
        
   destination.close()
   source.close()
	   
   print("Success.")

#------------------------------------------------------
main()
