import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import gca
import matplotlib.font_manager
import sys
from scipy import interpolate
from scipy.special import comb

import matplotlib.pylab as pylab
params = {'legend.fontsize': 'xx-large',
         'axes.labelsize': 'xx-large',
         'axes.titlesize': 'xx-large',
         'xtick.labelsize': 'xx-large',
         'ytick.labelsize': 'xx-large'}
pylab.rcParams.update(params)

# python3 plots/paper/semi_conv.py

def bernstein_poly(i, n, t):
    """
     The Bernstein polynomial of n, i as a function of t
    """

    return comb(n, i) * ( t**(n-i) ) * (1 - t)**i

def bezier_curve(points, nTimes=1000):
    """
       Given a set of control points, return the
       bezier curve defined by the control points.

       points should be a list of lists, or list of tuples
       such as [ [1,1], 
                 [2,3], 
                 [4,5], ..[Xn, Yn] ]
        nTimes is the number of time steps, defaults to 1000

        See http://processingjs.nihongoresources.com/bezierinfo/
    """

    nPoints = len(points)
    xPoints = np.array([p[0] for p in points])
    yPoints = np.array([p[1] for p in points])

    t = np.linspace(0.0, 1.0, nTimes)

    polynomial_array = np.array([ bernstein_poly(i, nPoints-1, t) for i in range(0, nPoints)   ])

    xvals = np.dot(xPoints, polynomial_array)
    yvals = np.dot(yPoints, polynomial_array)

    return xvals, yvals

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# x = [0, 2, 5, 9, 14, 20, 25, 30, 35, 40, 45]
# y = [100, 70, 40, 20, 10, 11, 12.5, 14, 15.5, 17.5, 19]
x = [0, 2, 5, 9, 14, 20, 25, 30, 35, 40, 45, 55, 70]
y = [100, 70, 40, 20, 8, 10, 12.5, 14, 15.5, 17.5, 19, 20.5, 22]
y = [i/100 for i in y]

points = list(zip(x, y))
xvals, yvals = bezier_curve(points, nTimes=1000)

fig = plt.figure(figsize=(6, 6))

# plt.scatter(x, y)
x_new = np.linspace(min(x), max(x), 100)
bspline = interpolate.make_interp_spline(x, y)
y_new = bspline(x_new)
# plt.plot(x_new, y_new, linewidth=2)
plt.plot(xvals, yvals, linewidth=2)
plt.grid()
plt.xlim([min(x), max(x)])
plt.ylim([0, max(y)+0.1])
plt.ylabel(r'$\|x^{(k)}-\overline{x}\|$')
plt.xlabel(r'Iterations')
# plt.tick_params(left=False,
#                 bottom=False,
#                 labelleft=False,
#                 labelbottom=False)
# plt.tick_params(labelleft=False,
#                 labelbottom=False)

filename_fig = "semi_conv"

plt.show()
fig.savefig("plots/paper/pdf/"+filename_fig+".pdf", bbox_inches='tight')
fig.savefig("plots/paper/png/"+filename_fig+".png", bbox_inches='tight')