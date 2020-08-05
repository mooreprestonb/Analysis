#!/usr/bin/python3
# smooth data

import numpy
import scipy
from scipy import signal

data = numpy.loadtxt("gr.dat")

#print(data)


data1 = numpy.array(data)
# First, design the Butterworth filter
N  = 3    # Filter order
Wn = .3 # Cutoff frequency
B, A = signal.butter(N, Wn, output='ba')

data1[:,1] = signal.filtfilt(B,A, data1[:,1])
data1[:,2] = signal.filtfilt(B,A, data1[:,2])

#data1[:,1] = signal.savgol_filter(data1[:,1],7,3)
#data1[:,2] = signal.savgol_filter(data1[:,2],7,3)

numpy.savetxt("grs.dat",data1)
#print(data1)
