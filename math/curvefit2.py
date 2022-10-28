#!/usr/bin/python3
# fit to polynomial function

import sys              # use to get command name from system
import argparse         # parses commandline arguments
import numpy            # numerical routines
from scipy.optimize import curve_fit  # fitting routines
import matplotlib.pyplot as plt       # ploting routines 

def fpoly(x,*p): # polynomial function f = p[0] + p[1]*x + p[2]*x**2 + ...
    f = 0
    for i in range(args.order+1):
        f += p[i]*x**i
    return(f)

# set command line arguments (for help use % curvefit -h)
parser = argparse.ArgumentParser(description="Fit x y data in file to polynomial")
parser.add_argument('input', help='input file (data in x y format')
parser.add_argument('output',nargs='?',help='output file for fitted function')
parser.add_argument('-xcolumn', type=int, dest='xc', default=1,
		    help="column to use for x data for fitting (Default 1)")
parser.add_argument('-ycolumn', type=int, dest='yc', default=2,
		    help="column to use for y data for fitting (Default 2)")
parser.add_argument('-o', type=int, dest='order', default=2,
		    help="order of polynomial to fit (Default 2 e.g. quadratic)")
parser.add_argument('-p', default=False, action="store_true", dest="p",
		    help="Show graphs (Default off)")
parser.add_argument('-var', default=False, action="store_true", dest="v",
		    help="Show covariance matrix (Default off)")

args = parser.parse_args()  # get command line arguments

data = numpy.loadtxt(args.input,usecols=(args.xc-1,args.yc-1)) # load in data from input file
xdata = data[:,0] # slice out x data
ydata = data[:,1] # slice out y data
#print (xdata,ydata) 

p = numpy.ones(args.order+1) # initial guess of all ones for polynomial
popt,pconv = curve_fit(fpoly,xdata,ydata,p0=p) # fit to poly
print("Order=",args.order,"\nCoef =",popt) # print out coeff and covariance matrix
if(args.v):
    print("Covar matrix =\n",pconv)

fit = fpoly(xdata,*popt)  # create array from fitted function
res = ydata-fit           # create residual (error) between fit function and data

# create header for file
head = str(sys.argv) + " x y fit res(y-fit)\n"
strpoly = "f = " + '{:g}'.format(popt[0])
for i in range(args.order):
    strpoly += "+" + '{:g}'.format(popt[i+1]) + "*x^" + str(i+1)
head += strpoly + " " + str(popt)

# output with x, y, fit, and residual
if(args.output):
    numpy.savetxt(args.output,numpy.column_stack((xdata,ydata,fit,res)),fmt='%12.8g',header=head)
else:
    print(head)
    print(numpy.column_stack((xdata,ydata,fit,res)))

# show graphs if requested
if(args.p):
    plt.figure()
    plt.subplot(211)
    plt.plot(xdata,ydata,'bo',label='data') # load in data
    plt.plot(xdata,fit,'r-',label='fit, order='+str(args.order))  # fit data
    plt.ylabel('y data')
    plt.legend()
    
    plt.subplot(212)  # residual
    plt.plot(xdata,res,label='residual')
    plt.ylabel('(data-fit)')
    plt.legend()
    
    plt.xlabel('x data')
    
    plt.show()
