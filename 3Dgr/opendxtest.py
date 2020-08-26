#/usr/bin/python3

nx = 20
ny = 10
nz = 15
xorg = -5
yorg = -5
zorg = -5

values = ""
print("# opendx file for testing with vmd")
print("object 1 class gridpositions counts",nx,ny,nz)
print("origin ",xorg,yorg,zorg)
print("delta .5 0 0")
print("delta 0 .5 0")
print("delta 0 0 .5")
print("object 2 class gridconnections counts",nx,ny,nz)
print("object 3 class array type double rank 0 items",nx*ny*nz," data follows")
n = 0
for i in range(nx):
    for j in range(ny):
        for k in range(nz):
            print(float(i)*float(j)*float(k))
            n += 1
print("object density class field")

