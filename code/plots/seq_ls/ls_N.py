import numpy as np
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import gca
import matplotlib.font_manager
import sys

# python3 plots/seq_ls/ls_N.py

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

indices = (10, 16, 24, 33)
it_CGLS_1000 = [it_CGLS[i] for i in indices]
time_CGLS_1000 = [time_CGLS[i] for i in indices]
indices = (18, 26, 35)
it_CGLS_4000 = [it_CGLS[i] for i in indices]
time_CGLS_4000 = [time_CGLS[i] for i in indices]
indices = (27, 36)
it_CGLS_10000 = [it_CGLS[i] for i in indices]
time_CGLS_10000 = [time_CGLS[i] for i in indices]

indices = (10, 16, 24, 33)
it_REK_1000 = [it_REK[i] for i in indices]
time_REK_1000 = [time_REK[i] for i in indices]
it_RGS_1000 = [it_RGS[i] for i in indices]
time_RGS_1000 = [time_RGS[i] for i in indices]
x_1000 = [4000, 20000, 40000, 80000]

indices = (18, 26, 35)
it_REK_4000 = [it_REK[i] for i in indices]
time_REK_4000 = [time_REK[i] for i in indices]
it_RGS_4000 = [it_RGS[i] for i in indices]
time_RGS_4000 = [time_RGS[i] for i in indices]
x_4000 = [20000, 40000, 80000]

indices = (27, 36)
it_REK_10000 = [it_REK[i] for i in indices]
time_REK_10000 = [time_REK[i] for i in indices]
it_RGS_10000 = [it_RGS[i] for i in indices]
time_RGS_10000 = [time_RGS[i] for i in indices]
x_10000 = [40000, 80000]

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
plt.plot(x_1000, time_RGS_1000, color='red', linewidth=2, label=r'RGS')
plt.scatter(x_1000, time_RGS_1000, color='red')
plt.plot(x_1000, time_REK_1000, color='blue', linewidth=2, label=r'REK')
plt.scatter(x_1000, time_REK_1000, color='blue')
plt.plot(x_1000, time_CGLS_1000, color='orange', linewidth=2, label=r'CGLS')
plt.scatter(x_1000, time_CGLS_1000, color='orange')
plt.grid()
plt.legend()
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$m$')
plt.ylabel(r'Total Time (s)')

filename_fig = "plots/seq_ls/ls_N_time_1000"

plt.show()
fig.savefig(filename_fig+".pdf", bbox_inches='tight')
fig.savefig(filename_fig+".png", bbox_inches='tight')

fig = plt.figure(figsize=(10, 7))
plt.plot(x_4000, time_RGS_4000, color='red', linewidth=2, label=r'RGS')
plt.scatter(x_4000, time_RGS_4000, color='red')
plt.plot(x_4000, time_REK_4000, color='blue', linewidth=2, label=r'REK')
plt.scatter(x_4000, time_REK_4000, color='blue')
plt.plot(x_4000, time_CGLS_4000, color='orange', linewidth=2, label=r'CGLS')
plt.scatter(x_4000, time_CGLS_4000, color='orange')
plt.grid()
plt.legend()
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$m$')
plt.ylabel(r'Total Time (s)')

filename_fig = "plots/seq_ls/ls_N_time_4000"

plt.show()
fig.savefig(filename_fig+".pdf", bbox_inches='tight')
fig.savefig(filename_fig+".png", bbox_inches='tight')

fig = plt.figure(figsize=(10, 7))
plt.plot(x_10000, time_RGS_10000, color='red', linewidth=2, label=r'RGS')
plt.scatter(x_10000, time_RGS_10000, color='red')
plt.plot(x_10000, time_REK_10000, color='blue', linewidth=2, label=r'REK')
plt.scatter(x_10000, time_REK_10000, color='blue')
plt.plot(x_10000, time_CGLS_10000, color='orange', linewidth=2, label=r'CGLS')
plt.scatter(x_10000, time_CGLS_10000, color='orange')
plt.grid()
plt.legend()
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$m$')
plt.ylabel(r'Total Time (s)')

filename_fig = "plots/seq_ls/ls_N_time_10000"

plt.show()
fig.savefig(filename_fig+".pdf", bbox_inches='tight')
fig.savefig(filename_fig+".png", bbox_inches='tight')

fig = plt.figure(figsize=(10, 7))
plt.plot(x_1000, it_RGS_1000, color='red', linewidth=2, label=r'RGS')
plt.scatter(x_1000, it_RGS_1000, color='red')
plt.plot(x_1000, it_REK_1000, color='blue', linewidth=2, label=r'REK')
plt.scatter(x_1000, it_REK_1000, color='blue')
plt.plot(x_1000, it_CGLS_1000, color='orange', linewidth=2, label=r'CGLS')
plt.scatter(x_1000, it_CGLS_1000, color='orange')
plt.grid()
plt.legend()
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$m$')
plt.ylabel(r'Iterations')

filename_fig = "plots/seq_ls/ls_N_it_1000"

plt.show()
fig.savefig(filename_fig+".pdf", bbox_inches='tight')
fig.savefig(filename_fig+".png", bbox_inches='tight')

fig = plt.figure(figsize=(10, 7))
plt.plot(x_4000, it_RGS_4000, color='red', linewidth=2, label=r'RGS')
plt.scatter(x_4000, it_RGS_4000, color='red')
plt.plot(x_4000, it_REK_4000, color='blue', linewidth=2, label=r'REK')
plt.scatter(x_4000, it_REK_4000, color='blue')
plt.plot(x_4000, it_CGLS_4000, color='orange', linewidth=2, label=r'CGLS')
plt.scatter(x_4000, it_CGLS_4000, color='orange')
plt.grid()
plt.legend()
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$m$')
plt.ylabel(r'Iterations')

filename_fig = "plots/seq_ls/ls_N_it_4000"

plt.show()
fig.savefig(filename_fig+".pdf", bbox_inches='tight')
fig.savefig(filename_fig+".png", bbox_inches='tight')

fig = plt.figure(figsize=(10, 7))
plt.plot(x_10000, it_RGS_10000, color='red', linewidth=2, label=r'RGS')
plt.scatter(x_10000, it_RGS_10000, color='red')
plt.plot(x_10000, it_REK_10000, color='blue', linewidth=2, label=r'REK')
plt.scatter(x_10000, it_REK_10000, color='blue')
plt.plot(x_10000, it_CGLS_10000, color='orange', linewidth=2, label=r'CGLS')
plt.scatter(x_10000, it_CGLS_10000, color='orange')
plt.grid()
plt.legend()
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$m$')
plt.ylabel(r'Iterations')

filename_fig = "plots/seq_ls/ls_N_it_10000"

plt.show()
fig.savefig(filename_fig+".pdf", bbox_inches='tight')
fig.savefig(filename_fig+".png", bbox_inches='tight')

fig = plt.figure(figsize=(10, 7))
plt.plot(x_1000, time_RGS_1000, color='orange', linewidth=2, label=r'RGS - $n=1000$')
plt.scatter(x_1000, time_RGS_1000, color='orange')
plt.plot(x_1000, time_REK_1000, color='orange', linestyle='--', linewidth=2, label=r'REK - $n=1000$')
plt.scatter(x_1000, time_REK_1000, color='orange')
plt.plot(x_1000, time_CGLS_1000, color='orange', linestyle='-.', linewidth=2, label=r'CGLS - $n=1000$')
plt.scatter(x_1000, time_CGLS_1000, color='orange')
plt.plot(x_4000, time_RGS_4000, color='red', linewidth=2, label=r'RGS - $n=4000$')
plt.scatter(x_4000, time_RGS_4000, color='red')
plt.plot(x_4000, time_REK_4000, color='red', linestyle='--', linewidth=2, label=r'REK - $n=4000$')
plt.scatter(x_4000, time_REK_4000, color='red')
plt.plot(x_4000, time_CGLS_4000, color='red', linestyle='-.', linewidth=2, label=r'CGLS - $n=4000$')
plt.scatter(x_4000, time_CGLS_4000, color='red')
plt.plot(x_10000, time_RGS_10000, color='blue', linewidth=2, label=r'RGS - $n=10000$')
plt.scatter(x_10000, time_RGS_10000, color='blue')
plt.plot(x_10000, time_REK_10000, color='blue', linestyle='--', linewidth=2, label=r'REK - $n=10000$')
plt.scatter(x_10000, time_REK_10000, color='blue')
plt.plot(x_10000, time_CGLS_10000, color='blue', linestyle='-.', linewidth=2, label=r'CGLS - $n=10000$')
plt.scatter(x_10000, time_CGLS_10000, color='blue')
plt.grid()
plt.legend()
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$m$')
plt.ylabel(r'Total Time (s)')

filename_fig = "plots/seq_ls/ls_N_time_all"

plt.show()
fig.savefig(filename_fig+".pdf", bbox_inches='tight')
fig.savefig(filename_fig+".png", bbox_inches='tight')

fig = plt.figure(figsize=(10, 7))
plt.plot(x_1000, it_RGS_1000, color='orange', linewidth=2, label=r'RGS - $n=1000$')
plt.scatter(x_1000, it_RGS_1000, color='orange')
plt.plot(x_1000, it_REK_1000, color='orange', linestyle='--', linewidth=2, label=r'REK - $n=1000$')
plt.scatter(x_1000, it_REK_1000, color='orange')
plt.plot(x_1000, it_CGLS_1000, color='orange', linestyle='-.', linewidth=2, label=r'CGLS - $n=1000$')
plt.scatter(x_1000, it_CGLS_1000, color='orange')
plt.plot(x_4000, it_RGS_4000, color='red', linewidth=2, label=r'RGS - $n=4000$')
plt.scatter(x_4000, it_RGS_4000, color='red')
plt.plot(x_4000, it_REK_4000, color='red', linestyle='--', linewidth=2, label=r'REK - $n=4000$')
plt.scatter(x_4000, it_REK_4000, color='red')
plt.plot(x_4000, it_CGLS_4000, color='red', linestyle='-.', linewidth=2, label=r'CGLS - $n=4000$')
plt.scatter(x_4000, it_CGLS_4000, color='red')
plt.plot(x_10000, it_RGS_10000, color='blue', linewidth=2, label=r'RGS - $n=10000$')
plt.scatter(x_10000, it_RGS_10000, color='blue')
plt.plot(x_10000, it_REK_10000, color='blue', linestyle='--', linewidth=2, label=r'REK - $n=10000$')
plt.scatter(x_10000, it_REK_10000, color='blue')
plt.plot(x_10000, it_CGLS_10000, color='blue', linestyle='-.', linewidth=2, label=r'CGLS - $n=10000$')
plt.scatter(x_10000, it_CGLS_10000, color='blue')
plt.grid()
plt.legend()
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$m$')
plt.ylabel(r'Iterations')

filename_fig = "plots/seq_ls/ls_N_it_all"

plt.show()
fig.savefig(filename_fig+".pdf", bbox_inches='tight')
fig.savefig(filename_fig+".png", bbox_inches='tight')