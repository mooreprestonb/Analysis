#!/usr/bin/python3

import numpy
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import scipy.interpolate

data = numpy.loadtxt('NaCldens.dat',usecols=range(7))
temp = numpy.array([10,25,40,55,70,85])
mol = data[:,0]
data = numpy.delete(data,(0),axis=1)
#print(temp,mol,data,temp[2],mol[10],data[10][2])

data2 = numpy.loadtxt('wghtperNaCl.dat',usecols=range(8))
temp2 = numpy.array([0,10,25,40,60,80,100])
perc = data2[:,0]
data2 = numpy.delete(data2,(0),axis=1)
print("Weight %:",perc)
for i in range(len(perc)):
    perc[i] = ((perc[i]*1000)/(22.990+35.453))/(100-perc[i])
print("Molality:",perc)

xv, yv = numpy.meshgrid(temp,mol)
xv2, yv2 = numpy.meshgrid(temp2,perc)

# ScatterPlot.
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(xv, yv, data)
ax.scatter(xv2,yv2,data2)
#ax.plot_wireframe(xv,yv,data,rstride=len(mol),cstride=len(temp))
#ax.plot_wireframe(xv,yv,data,color='b')
#ax.plot_wireframe(xv2,yv2,data2,color='r')

ax.set_xlabel('temp')
ax.set_ylabel('molality')
ax.set_zlabel('density')
plt.show()

#print(xv,yv,data)

#spline = scipy.interpolate.SmoothBivariateSpline(xv.flatten(),yv.flatten(),data.flatten())
#lin = scipy.interpolate.interp2d(temp,mol,data)
cub = scipy.interpolate.interp2d(temp,mol,data,kind='cubic')
cub2 = scipy.interpolate.interp2d(temp2,perc,data2,kind='cubic')
#print(mxy)
while(True):
    temp = float(input("Enter temp: "))
    molality = float(input("Enter molality: "))
    #    mxy = spline.ev(temp,molality)
    #lxy = lin(temp,molality)
    cxy = cub(temp,molality)
    cxy2 = cub2(temp,molality)
    print(temp,molality,cxy,cxy2)
