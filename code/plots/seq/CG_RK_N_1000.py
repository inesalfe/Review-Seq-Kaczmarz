import numpy as np
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import gca
import matplotlib.font_manager
import sys

# python3 plots/seq/CG_RK_N_1000.py

x_1000 = [4000, 20000, 40000, 80000, 160000]

filename = "outputs/seq/cg_N_1000.txt";

with open(filename) as f:
	lines = f.read().splitlines()

file_size = len(lines)

time = []
for i in range(5):
	time.append(float(lines[i].split()[3]))

cg_time_1000 = time

filename = "outputs/seq/cgls_N_1000.txt";

with open(filename) as f:
	lines = f.read().splitlines()

file_size = len(lines)

time = []
for i in range(5):
	time.append(float(lines[i].split()[3]))

cgls_time_1000 = time

filename = "outputs/seq/RK.txt";

with open(filename) as f:
	lines = f.read().splitlines()

file_size = len(lines)

time = []
for i in range(file_size):
	time.append(float(lines[i].split()[2]))

indices = (10, 16, 24, 33, 43)
RK_time_1000 = [time[i] for i in indices]

import matplotlib.pylab as pylab
params = {'legend.fontsize': 'xx-large',
         'axes.labelsize': 'xx-large',
         'axes.titlesize': 'xx-large',
         'xtick.labelsize': 'xx-large',
         'ytick.labelsize': 'xx-large'}
pylab.rcParams.update(params)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plot_title = r'Overdetermined systems with $n=1000$ sampled from $N(\mu,\sigma)$'

fig = plt.figure(figsize=(10, 7))
plt.plot(x_1000, cg_time_1000, color='green', linewidth=2)
plt.scatter(x_1000, cg_time_1000, color='green', label=r'CG', marker='^', s=75)
plt.plot(x_1000, cgls_time_1000, color='orange', linewidth=2)
plt.scatter(x_1000, cgls_time_1000, color='orange', label=r'CGLS', marker='s', s=75)
plt.plot(x_1000, RK_time_1000, color='red', linewidth=2)
plt.scatter(x_1000, RK_time_1000, color='red', label=r'RK', marker='o', s=75)
plt.grid()
plt.legend()
# plt.title(plot_title)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$m$')
plt.ylabel(r'Total Time (s)')

filename_fig = "CG_RK_N_1000_time"

plt.show()
fig.savefig("plots/seq/pdf/"+filename_fig+".pdf", bbox_inches='tight')
fig.savefig("plots/seq/png/"+filename_fig+".png", bbox_inches='tight')