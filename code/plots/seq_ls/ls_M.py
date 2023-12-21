import numpy as np
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import gca
import matplotlib.font_manager
import sys

# python3 plots/seq_ls/ls_M.py

filename = "outputs/seq_ls/REK_ls.txt";

with open(filename) as f:
	lines = f.read().splitlines()

file_size = len(lines)

time_REK = []
it_REK = []
for i in range(file_size):
	time_REK.append(float(lines[i].split()[2]))
	it_REK.append(int(lines[i].split()[3]))

filename = "outputs/seq_ls/RGS_ls.txt";

with open(filename) as f:
	lines = f.read().splitlines()

file_size = len(lines)

time_RGS = []
it_RGS = []
for i in range(file_size):
	time_RGS.append(float(lines[i].split()[2]))
	it_RGS.append(int(lines[i].split()[3]))

filename = "outputs/seq_ls/cgls_ls.txt";

with open(filename) as f:
	lines = f.read().splitlines()

file_size = len(lines)

time_CGLS = []
it_CGLS = []
for i in range(file_size):
	time_CGLS.append(float(lines[i].split()[3]))
	it_CGLS.append(int(lines[i].split()[2]))

indices = (11, 12, 13, 14, 15, 16, 17, 18)
it_CGLS_20000 = [it_CGLS[i] for i in indices]
time_CGLS_20000 = [time_CGLS[i] for i in indices]
indices = (19, 20, 21, 22, 23, 24, 25, 26, 27)
it_CGLS_40000 = [it_CGLS[i] for i in indices]
time_CGLS_40000 = [time_CGLS[i] for i in indices]
indices = (28, 29, 30, 31, 32, 33, 34, 35, 36)
it_CGLS_80000 = [it_CGLS[i] for i in indices]
time_CGLS_80000 = [time_CGLS[i] for i in indices]

indices = (11, 12, 13, 14, 15, 16, 17, 18)
it_REK_20000 = [it_REK[i] for i in indices]
time_REK_20000 = [time_REK[i] for i in indices]
it_RGS_20000 = [it_RGS[i] for i in indices]
time_RGS_20000 = [time_RGS[i] for i in indices]
x_20000 = [50, 100, 200, 500, 750, 1000, 2000, 4000]

indices = (19, 20, 21, 22, 23, 24, 25, 26, 27)
it_REK_40000 = [it_REK[i] for i in indices]
time_REK_40000 = [time_REK[i] for i in indices]
it_RGS_40000 = [it_RGS[i] for i in indices]
time_RGS_40000 = [time_RGS[i] for i in indices]
x_40000 = [50, 100, 200, 500, 750, 1000, 2000, 4000, 10000]

indices = (28, 29, 30, 31, 32, 33, 34, 35, 36)
it_REK_80000 = [it_REK[i] for i in indices]
time_REK_80000 = [time_REK[i] for i in indices]
it_RGS_80000 = [it_RGS[i] for i in indices]
time_RGS_80000 = [time_RGS[i] for i in indices]
x_80000 = [50, 100, 200, 500, 750, 1000, 2000, 4000, 10000]

import matplotlib.pylab as pylab
params = {'legend.fontsize': 'xx-large',
         'axes.labelsize': 'xx-large',
         'axes.titlesize': 'xx-large',
         'xtick.labelsize': 'xx-large',
         'ytick.labelsize': 'xx-large'}
pylab.rcParams.update(params)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

fig = plt.figure(figsize=(10, 7))
plt.plot(x_20000, time_RGS_20000, color='red', linewidth=2, label=r'RGS')
plt.scatter(x_20000, time_RGS_20000, color='red')
plt.plot(x_20000, time_REK_20000, color='blue', linewidth=2, label=r'REK')
plt.scatter(x_20000, time_REK_20000, color='blue')
plt.plot(x_20000, time_CGLS_20000, color='orange', linewidth=2, label=r'CGLS')
plt.scatter(x_20000, time_CGLS_20000, color='orange')
plt.grid()
plt.legend()
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$n$')
plt.ylabel(r'Total Time (s)')

filename_fig = "plots/seq_ls/ls_M_time_20000"

plt.show()
fig.savefig(filename_fig+".pdf", bbox_inches='tight')
fig.savefig(filename_fig+".png", bbox_inches='tight')

fig = plt.figure(figsize=(10, 7))
plt.plot(x_40000, time_RGS_40000, color='red', linewidth=2, label=r'RGS')
plt.scatter(x_40000, time_RGS_40000, color='red')
plt.plot(x_40000, time_REK_40000, color='blue', linewidth=2, label=r'REK')
plt.scatter(x_40000, time_REK_40000, color='blue')
plt.plot(x_40000, time_CGLS_40000, color='orange', linewidth=2, label=r'CGLS')
plt.scatter(x_40000, time_CGLS_40000, color='orange')
plt.grid()
plt.legend()
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$n$')
plt.ylabel(r'Total Time (s)')

filename_fig = "plots/seq_ls/ls_M_time_40000"

plt.show()
fig.savefig(filename_fig+".pdf", bbox_inches='tight')
fig.savefig(filename_fig+".png", bbox_inches='tight')

fig = plt.figure(figsize=(10, 7))
plt.plot(x_80000, time_RGS_80000, color='red', linewidth=2, label=r'RGS')
plt.scatter(x_80000, time_RGS_80000, color='red')
plt.plot(x_80000, time_REK_80000, color='blue', linewidth=2, label=r'REK')
plt.scatter(x_80000, time_REK_80000, color='blue')
plt.plot(x_80000, time_CGLS_80000, color='orange', linewidth=2, label=r'CGLS')
plt.scatter(x_80000, time_CGLS_80000, color='orange')
plt.grid()
plt.legend()
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$n$')
plt.ylabel(r'Total Time (s)')

filename_fig = "plots/seq_ls/ls_M_time_80000"

plt.show()
fig.savefig(filename_fig+".pdf", bbox_inches='tight')
fig.savefig(filename_fig+".png", bbox_inches='tight')

fig = plt.figure(figsize=(10, 7))
plt.plot(x_20000, it_RGS_20000, color='red', linewidth=2, label=r'RGS')
plt.scatter(x_20000, it_RGS_20000, color='red')
plt.plot(x_20000, it_REK_20000, color='blue', linewidth=2, label=r'REK')
plt.scatter(x_20000, it_REK_20000, color='blue')
plt.plot(x_20000, it_CGLS_20000, color='orange', linewidth=2, label=r'CGLS')
plt.scatter(x_20000, it_CGLS_20000, color='orange')
plt.grid()
plt.legend()
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$n$')
plt.ylabel(r'Iterations')

filename_fig = "plots/seq_ls/ls_M_it_20000"

plt.show()
fig.savefig(filename_fig+".pdf", bbox_inches='tight')
fig.savefig(filename_fig+".png", bbox_inches='tight')

fig = plt.figure(figsize=(10, 7))
plt.plot(x_40000, it_RGS_40000, color='red', linewidth=2, label=r'RGS')
plt.scatter(x_40000, it_RGS_40000, color='red')
plt.plot(x_40000, it_REK_40000, color='blue', linewidth=2, label=r'REK')
plt.scatter(x_40000, it_REK_40000, color='blue')
plt.plot(x_40000, it_CGLS_40000, color='orange', linewidth=2, label=r'CGLS')
plt.scatter(x_40000, it_CGLS_40000, color='orange')
plt.grid()
plt.legend()
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$n$')
plt.ylabel(r'Iterations')

filename_fig = "plots/seq_ls/ls_M_it_40000"

plt.show()
fig.savefig(filename_fig+".pdf", bbox_inches='tight')
fig.savefig(filename_fig+".png", bbox_inches='tight')

fig = plt.figure(figsize=(10, 7))
plt.plot(x_80000, it_RGS_80000, color='red', linewidth=2, label=r'RGS')
plt.scatter(x_80000, it_RGS_80000, color='red')
plt.plot(x_80000, it_REK_80000, color='blue', linewidth=2, label=r'REK')
plt.scatter(x_80000, it_REK_80000, color='blue')
plt.plot(x_80000, it_CGLS_80000, color='orange', linewidth=2, label=r'CGLS')
plt.scatter(x_80000, it_CGLS_80000, color='orange')
plt.grid()
plt.legend()
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$n$')
plt.ylabel(r'Iterations')

filename_fig = "plots/seq_ls/ls_M_it_80000"

plt.show()
fig.savefig(filename_fig+".pdf", bbox_inches='tight')
fig.savefig(filename_fig+".png", bbox_inches='tight')

fig = plt.figure(figsize=(10, 7))
plt.plot(x_20000, time_RGS_20000, color='orange', linewidth=2, label=r'RGS - $m=20000$')
plt.scatter(x_20000, time_RGS_20000, color='orange')
plt.plot(x_20000, time_REK_20000, color='orange', linestyle='--', linewidth=2, label=r'REK - $m=20000$')
plt.scatter(x_20000, time_REK_20000, color='orange')
plt.plot(x_20000, time_CGLS_20000, color='orange', linestyle='-.', linewidth=2, label=r'CGLS - $m=20000$')
plt.scatter(x_20000, time_CGLS_20000, color='orange')
plt.plot(x_40000, time_RGS_40000, color='red', linewidth=2, label=r'RGS - $m=40000$')
plt.scatter(x_40000, time_RGS_40000, color='red')
plt.plot(x_40000, time_REK_40000, color='red', linestyle='--', linewidth=2, label=r'REK - $m=40000$')
plt.scatter(x_40000, time_REK_40000, color='red')
plt.plot(x_40000, time_CGLS_40000, color='red', linestyle='-.', linewidth=2, label=r'CGLS - $m=40000$')
plt.scatter(x_40000, time_CGLS_40000, color='red')
plt.plot(x_80000, time_RGS_80000, color='blue', linewidth=2, label=r'RGS - $m=80000$')
plt.scatter(x_80000, time_RGS_80000, color='blue')
plt.plot(x_80000, time_REK_80000, color='blue', linestyle='--', linewidth=2, label=r'REK - $m=80000$')
plt.scatter(x_80000, time_REK_80000, color='blue')
plt.plot(x_80000, time_CGLS_80000, color='blue', linestyle='-.', linewidth=2, label=r'CGLS - $m=80000$')
plt.scatter(x_80000, time_CGLS_80000, color='blue')
plt.grid()
plt.legend()
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$n$')
plt.ylabel(r'Total Time (s)')

filename_fig = "plots/seq_ls/ls_M_time_all"

plt.show()
fig.savefig(filename_fig+".pdf", bbox_inches='tight')
fig.savefig(filename_fig+".png", bbox_inches='tight')

fig = plt.figure(figsize=(10, 7))
plt.plot(x_20000, it_RGS_20000, color='orange', linewidth=2, label=r'RGS - $m=20000$')
plt.scatter(x_20000, it_RGS_20000, color='orange')
plt.plot(x_20000, it_REK_20000, color='orange', linestyle='--', linewidth=2, label=r'REK - $m=20000$')
plt.scatter(x_20000, it_REK_20000, color='orange')
plt.plot(x_20000, it_CGLS_20000, color='orange', linestyle='-.', linewidth=2, label=r'CGLS - $m=20000$')
plt.scatter(x_20000, it_CGLS_20000, color='orange')
plt.plot(x_40000, it_RGS_40000, color='red', linewidth=2, label=r'RGS - $m=40000$')
plt.scatter(x_40000, it_RGS_40000, color='red')
plt.plot(x_40000, it_REK_40000, color='red', linestyle='--', linewidth=2, label=r'REK - $m=40000$')
plt.scatter(x_40000, it_REK_40000, color='red')
plt.plot(x_40000, it_CGLS_40000, color='red', linestyle='-.', linewidth=2, label=r'CGLS - $m=40000$')
plt.scatter(x_40000, it_CGLS_40000, color='red')
plt.plot(x_80000, it_RGS_80000, color='blue', linewidth=2, label=r'RGS - $m=80000$')
plt.scatter(x_80000, it_RGS_80000, color='blue')
plt.plot(x_80000, it_REK_80000, color='blue', linestyle='--', linewidth=2, label=r'REK - $m=80000$')
plt.scatter(x_80000, it_REK_80000, color='blue')
plt.plot(x_80000, it_CGLS_80000, color='blue', linestyle='-.', linewidth=2, label=r'CGLS - $m=80000$')
plt.scatter(x_80000, it_CGLS_80000, color='blue')
plt.grid()
plt.legend()
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$n$')
plt.ylabel(r'Iterations')

filename_fig = "plots/seq_ls/ls_M_it_all"

plt.show()
fig.savefig(filename_fig+".pdf", bbox_inches='tight')
fig.savefig(filename_fig+".png", bbox_inches='tight')