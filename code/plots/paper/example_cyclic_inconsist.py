import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import gca
import matplotlib.font_manager
import sys
from skspatial.objects import Line, Point
from skspatial.plotting import plot_2d

import matplotlib.pylab as pylab
params = {'legend.fontsize': 'x-large',
         'axes.labelsize': 'xx-large',
         'axes.titlesize': 'xx-large',
         'xtick.labelsize': 'xx-large',
         'ytick.labelsize': 'xx-large'}
pylab.rcParams.update(params)

# python3 plots/paper/example_cyclic_inconsist.py

x_init = [4, 5]

x_ls = [2.5, 3.2]

x_init_point = Point(x_init)

l1 = Line.from_points([0,-1], [3,7])
l2 = Line.from_points([7,0], [0,7])
l3 = Line.from_points([0,1], [6,4])
# l4 = Line.from_points([2,0], [0,4])

x1_point = l1.project_point(x_init_point)
x1 = np.array(x1_point)
x2_point = l2.project_point(x1_point)
x2 = np.array(x2_point)
x3_point = l3.project_point(x2_point)
x3 = np.array(x3_point)
x4_point = l1.project_point(x3_point)
x4 = np.array(x4_point)
x5_point = l2.project_point(x4_point)
x5 = np.array(x5_point)
x6_point = l3.project_point(x5_point)
x6 = np.array(x6_point)
x7_point = l1.project_point(x6_point)
x7 = np.array(x7_point)
x8_point = l2.project_point(x7_point)
x8 = np.array(x8_point)
x9_point = l3.project_point(x8_point)
x9 = np.array(x9_point)
x9_point = l1.project_point(x9_point)
x10 = np.array(x9_point)
x10_point = l2.project_point(x9_point)
x11 = np.array(x10_point)
x11_point = l3.project_point(x10_point)
x12 = np.array(x11_point)
x12_point = l1.project_point(x11_point)
x13 = np.array(x12_point)
x13_point = l2.project_point(x12_point)
x14 = np.array(x13_point)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

fig = plt.figure(figsize=(6, 6))
plt.plot([0,3], [-1,7], color='red', label=r'$\langle a_1,x \rangle = b_1$')
plt.plot([7,0], [0,7], color='green', label=r'$\langle a_2,x \rangle = b_2$')
plt.plot([0,6], [1,4], color='blue', label=r'$\langle a_3,x \rangle = b_3$')
# plt.plot([2,0], [0,4], color='magenta', label=r'$\langle a_4,x \rangle = b_4$')
plt.plot([x_init[0], x1[0]], [x_init[1], x1[1]], color='grey', linestyle='--')
plt.plot([x1[0], x2[0]], [x1[1], x2[1]], color='grey', linestyle='--')
plt.plot([x2[0], x3[0]], [x2[1], x3[1]], color='grey', linestyle='--')
plt.plot([x3[0], x4[0]], [x3[1], x4[1]], color='grey', linestyle='--')
plt.plot([x4[0], x5[0]], [x4[1], x5[1]], color='grey', linestyle='--')
plt.plot([x5[0], x6[0]], [x5[1], x6[1]], color='grey', linestyle='--')
plt.plot([x6[0], x7[0]], [x6[1], x7[1]], color='grey', linestyle='--')
plt.plot([x7[0], x8[0]], [x7[1], x8[1]], color='grey', linestyle='--')
plt.plot([x8[0], x9[0]], [x8[1], x9[1]], color='grey', linestyle='--')
plt.plot([x9[0], x10[0]], [x9[1], x10[1]], color='grey', linestyle='--')
plt.plot([x10[0], x11[0]], [x10[1], x11[1]], color='grey', linestyle='--')
plt.plot([x11[0], x12[0]], [x11[1], x12[1]], color='grey', linestyle='--')
plt.plot([x12[0], x13[0]], [x12[1], x13[1]], color='grey', linestyle='--')
plt.plot([x13[0], x14[0]], [x13[1], x14[1]], color='grey', linestyle='--')
plt.scatter(x_init[0], x_init[1], color='c', label=r'$x^{(0)}$')
plt.scatter(x1[0], x1[1], color='black', label=r'$x^{(k)}$', zorder=2)
plt.scatter(x2[0], x2[1], color='black', zorder=2)
plt.scatter(x3[0], x3[1], color='black', zorder=2)
plt.scatter(x4[0], x4[1], color='black', zorder=2)
plt.scatter(x5[0], x5[1], color='black', zorder=2)
plt.scatter(x6[0], x6[1], color='black', zorder=2)
plt.scatter(x7[0], x7[1], color='black', zorder=2)
plt.scatter(x8[0], x8[1], color='black', zorder=2)
plt.scatter(x9[0], x9[1], color='black', zorder=2)
plt.scatter(x10[0], x10[1], color='black', zorder=2)
plt.scatter(x11[0], x11[1], color='black', zorder=2)
plt.scatter(x12[0], x12[1], color='black', zorder=2)
plt.scatter(x13[0], x13[1], color='black', zorder=2)
plt.scatter(x14[0], x14[1], color='black', zorder=2)
plt.scatter(x_ls[0], x_ls[1], color='magenta', label=r'$x_{LS}$')
plt.grid()
plt.axis('equal')
plt.xlim([0, 6])
plt.ylim([0, 6])
plt.legend(loc='best')
plt.ylabel(r'$y$')
plt.xlabel(r'$x$')

filename_fig = "example_cyclic_inconsist"

plt.show()
fig.savefig("plots/paper/pdf/"+filename_fig+".pdf", bbox_inches='tight')
fig.savefig("plots/paper/png/"+filename_fig+".png", bbox_inches='tight')
