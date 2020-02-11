#!/usr/bin/python3
# <examples/doc_model_gaussian.py>

import matplotlib.pyplot as plt
from numpy import exp, linspace, random
from scipy.optimize import curve_fit

def gaussian(x, amp, cen, wid):
    return amp * exp(-(x-cen)**2 / wid)

x = linspace(-10, 10, 101)
y = gaussian(x, 2.33, 0.21, 1.51) + random.normal(0, 0.2, x.size)

ivals = [1, 0, 1]  # for [amp, cen, wid]
yi = gaussian(x,ivals[0],ivals[1],ivals[2])
fvals, covar = curve_fit(gaussian, x, y, p0=ivals)
print('best_vals: {}'.format(fvals))

yr = gaussian(x,fvals[0],fvals[1],fvals[2])
plt.plot(x, y, 'bo')
plt.plot(x, yi, 'g-', label='init')
plt.plot(x, yr, 'r-', label='fit')
plt.legend(loc='best')
plt.show()
# <end examples/doc_model_gaussian.py>
