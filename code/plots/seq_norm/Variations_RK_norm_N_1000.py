import numpy as np
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import gca
import matplotlib.font_manager
import sys

# python3 plots/seq_norm/Variations_RK_norm_N_1000.py

x_1000 = [4000, 20000, 40000, 80000]

filename = "outputs/seq_norm/GSSRK_norm_N_1000.txt";

with open(filename) as f:
	lines = f.read().splitlines()

file_size = len(lines)

it = []
time = []
for i in range(4):
	it.append(int(lines[i].split()[3]))
	time.append(float(lines[i].split()[2]))

GSSRK_it_1000 = it
GSSRK_time_1000 = time

filename = "outputs/seq_norm/GRK_norm_N_1000.txt";

with open(filename) as f:
	lines = f.read().splitlines()

file_size = len(lines)

it = []
time = []
for i in range(4):
	it.append(int(lines[i].split()[3]))
	time.append(float(lines[i].split()[2]))

GRK_it_1000 = it
GRK_time_1000 = time

filename = "outputs/seq_norm/RK_norm.txt";

with open(filename) as f:
	lines = f.read().splitlines()

file_size = len(lines)

it = []
time = []
for i in range(file_size):
	it.append(int(lines[i].split()[3]))
	time.append(float(lines[i].split()[2]))

indices = (10, 16, 24, 33)
RK_it_1000 = [it[i] for i in indices]
RK_time_1000 = [time[i] for i in indices]

filename = "outputs/seq_norm/NSSRK_norm_N_1000.txt";

with open(filename) as f:
	lines = f.read().splitlines()

file_size = len(lines)

it = []
time = []
for i in range(4):
	it.append(int(lines[i].split()[3]))
	time.append(float(lines[i].split()[2]))

NSSRK_it_1000 = it
NSSRK_time_1000 = time

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
plt.scatter(x_1000, GSSRK_time_1000, color='blue')
plt.plot(x_1000, GSSRK_time_1000, color='blue', linewidth=2, label=r'GSSRK')
plt.scatter(x_1000, GRK_time_1000, color='orange')
plt.plot(x_1000, GRK_time_1000, color='orange', linewidth=2, label=r'GRK')
plt.scatter(x_1000, NSSRK_time_1000, color='red')
plt.plot(x_1000, NSSRK_time_1000, color='red', linewidth=2, label=r'NSSRK')
plt.scatter(x_1000, RK_time_1000, color='green')
plt.plot(x_1000, RK_time_1000, color='green', linewidth=2, label=r'RK')
plt.grid()
plt.legend()
# plt.title(plot_title)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$m$')
plt.ylabel(r'Total Time (s)')

filename_fig = "Variations_RK_N_1000_time_norm"

plt.show()
fig.savefig("plots/seq_norm/pdf/"+filename_fig+".pdf", bbox_inches='tight')
fig.savefig("plots/seq_norm/png/"+filename_fig+".png", bbox_inches='tight')

plot_title = r'Overdetermined systems with $n=1000$ sampled from $N(\mu,\sigma)$'

fig = plt.figure(figsize=(10, 7))
plt.scatter(x_1000, GSSRK_it_1000, color='blue')
plt.plot(x_1000, GSSRK_it_1000, color='blue', linewidth=2, label=r'GSSRK')
plt.scatter(x_1000, GRK_it_1000, color='orange')
plt.plot(x_1000, GRK_it_1000, color='orange', linewidth=2, label=r'GRK')
plt.scatter(x_1000, NSSRK_it_1000, color='red')
plt.plot(x_1000, NSSRK_it_1000, color='red', linewidth=2, label=r'NSSRK')
plt.scatter(x_1000, RK_it_1000, color='green')
plt.plot(x_1000, RK_it_1000, color='green', linewidth=2, label=r'RK')
plt.grid()
plt.legend()
# plt.title(plot_title)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$m$')
plt.ylabel(r'Iterations')

filename_fig = "Variations_RK_N_1000_it_norm"

plt.show()
fig.savefig("plots/seq_norm/pdf/"+filename_fig+".pdf", bbox_inches='tight')
fig.savefig("plots/seq_norm/png/"+filename_fig+".png", bbox_inches='tight')