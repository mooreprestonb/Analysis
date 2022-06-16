#!/usr/bin/python3

# fit to quad function
import numpy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from scipy.optimize import curve_fit

def fun2d(x,y,*p):
    z = 0
    for i in range(order):
        for j in range(order):
            z += p[i*order+j] * x**i * y**j
    return(z)

def fit2d(M,*p):
    x,y = M  # unpack 2d array
    z = numpy.zeros(x.shape)
    for i in range(order):
        for j in range(order):
            z += p[i*order+j] * x**i * y**j
    return(z)

order = 4  # how well to fit data 
c = numpy.ones(order*order)

data = numpy.loadtxt('NaCldens.dat',usecols=range(7))
temp = numpy.array([10,25,40,55,70,85])
mol = data[:,0]
data = numpy.delete(data,(0),axis=1)

xv, yv = numpy.meshgrid(temp,mol)
xdata=numpy.vstack((xv.ravel(),yv.ravel()))
zdata=data.ravel()

popt,pconv = curve_fit(fit2d,xdata,zdata,p0=c)
print(popt)

# data2
data2 = numpy.loadtxt('wghtperNaCl.dat',usecols=range(8))
temp2 = numpy.array([0,10,25,40,60,80,100])
perc = data2[:,0]
data2 = numpy.delete(data2,(0),axis=1)
print("Weight %:",perc)
for i in range(len(perc)):
    perc[i] = ((perc[i]*1000)/(22.990+35.453))/(100-perc[i])
print("Molality:",perc)

xv2, yv2 = numpy.meshgrid(temp2,perc)
xdata2=numpy.vstack((xv2.ravel(),yv2.ravel()))
zdata2=data2.ravel()
popt2,pconv2 = curve_fit(fit2d,xdata2,zdata2,p0=c)
print(popt2)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(xv, yv, data,color='b')
#ax.scatter(xv2,yv2,data2)
ax.scatter(xv,yv,fit2d(xdata,*popt),color='g')
ax.scatter(xv2,yv2,fit2d(xdata2,*popt2),color='r')

ax.set_xlabel('temp')
ax.set_ylabel('molality')
ax.set_zlabel('density')
plt.show()

while(True):
    temp = float(input("Enter temp: "))
    molality = float(input("Enter molality: "))
    print(temp,molality,fun2d(temp,molality,*popt),fun2d(temp,molality,*popt2))
