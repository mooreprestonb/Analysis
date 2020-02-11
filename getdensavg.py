#!/usr/bin/python3

import sys
import numpy
import matplotlib.pyplot as plt

print (sys.argv)
data = numpy.loadtxt(sys.argv[1])
h = data[1][0] - data[0][0]
x = data[:,0]
dens = data[:,1]
densp = (dens[:-2]-dens[2:])/(2*h)

avg = numpy.mean(densp)
std = numpy.std(densp)

iface1=0
iface2=0
interface1 = []
interface2 = []
# find interface
for i in range(len(densp)):
    if abs(densp[i]) > std/16.:
        print (i,densp[i],x[i+1],dens[i+1],iface1,iface2)
        if(iface1==0) :
            iface1=1
            interface1.append(i)
        elif(interface1[-1] == i-1):
            interface1.append(i)
        elif(iface2==0):
            iface2 = 1
            interface2.append(i)
        elif(iface2==1 and interface2[-1] == i-1):
            interface2.append(i)
        else:
            print("Error found more then one interface!")
            exit(1)

#print(iface1,interface1,iface2,interface2)
if(iface2==0):
    print("Error didn't find two interfaces!")
    exit(1)

#plt.plot(x, dens, 'bo', label="Density")
#plt.plot(x[1:-1], densp, 'r-', label="d Density/dx")
#plt.legend(loc='best')
#plt.show()

avg1 = numpy.mean(dens[:interface1[0]-2])
std1 = numpy.std(dens[:interface1[0]-2])
avg2 = numpy.mean(dens[interface1[-1]+3:interface2[0]-1])
std2 = numpy.mean(dens[interface1[-1]+3:interface2[0]-1])
avg3 = numpy.mean(dens[interface2[-1]+5:])
std3= numpy.mean(dens[interface2[-1]+5:])

print(dens[0],dens[interface1[0]-2])
print(dens[interface1[-1]+3],dens[interface2[0]-1])
print(dens[interface2[-1]+5],dens[-1])

print(avg1,avg3,std1,std3,avg2,std2)


#numpy.savetxt("hist.dat",numpy.column_stack((xhist,hist)))



