import numpy as np
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import gca
import matplotlib.font_manager
import sys

# python3 plots/seq/Variations_RK_N_1000.py

x_1000 = [4000, 20000, 40000, 80000]

filename = "outputs/seq/GSSRK_N_1000.txt";

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

filename = "outputs/seq/GRK_N_1000.txt";

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

filename = "outputs/seq/RK.txt";

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

filename = "outputs/seq/NSSRK_N_1000.txt";

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
plt.scatter(x_1000, GSSRK_time_1000, color='blue', label=r'GSSRK', marker='s', s=75)
plt.plot(x_1000, GSSRK_time_1000, color='blue', linewidth=2)
plt.scatter(x_1000, GRK_time_1000, color='orange', label=r'GRK', marker='o', s=75)
plt.plot(x_1000, GRK_time_1000, color='orange', linewidth=2)
plt.scatter(x_1000, NSSRK_time_1000, color='red', label=r'NSSRK', marker='X', s=75)
plt.plot(x_1000, NSSRK_time_1000, color='red', linewidth=2)
plt.scatter(x_1000, RK_time_1000, color='green', label=r'RK', marker='^', s=75)
plt.plot(x_1000, RK_time_1000, color='green', linewidth=2)
plt.grid()
plt.legend()
# plt.title(plot_title)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$m$')
plt.ylabel(r'Total Time (s)')

filename_fig = "Variations_RK_N_1000_time"

plt.show()
fig.savefig("plots/seq/pdf/"+filename_fig+".pdf", bbox_inches='tight')
fig.savefig("plots/seq/png/"+filename_fig+".png", bbox_inches='tight')

plot_title = r'Overdetermined systems with $n=1000$ sampled from $N(\mu,\sigma)$'

fig = plt.figure(figsize=(10, 7))
plt.scatter(x_1000, GSSRK_it_1000, color='blue', label=r'GSSRK', marker='s', s=75)
plt.plot(x_1000, GSSRK_it_1000, color='blue', linewidth=2)
plt.scatter(x_1000, GRK_it_1000, color='orange', label=r'GRK', marker='o', s=75)
plt.plot(x_1000, GRK_it_1000, color='orange', linewidth=2)
plt.scatter(x_1000, NSSRK_it_1000, color='red', label=r'NSSRK', marker='X', s=75)
plt.plot(x_1000, NSSRK_it_1000, color='red', linewidth=2)
plt.scatter(x_1000, RK_it_1000, color='green', label=r'RK', marker='^', s=75)
plt.plot(x_1000, RK_it_1000, color='green', linewidth=2)
plt.grid()
plt.legend()
# plt.title(plot_title)
plt.yscale('log')
plt.xscale('log')
x_ticks = [4000, 10000, 20000, 40000, 80000]
plt.xticks(x_ticks)
ax = fig.gca()
ax.set_xticklabels(x_ticks)
from matplotlib.ticker import LogFormatterSciNotation
plt.gca().xaxis.set_major_formatter(LogFormatterSciNotation(base=10,minor_thresholds=(10,10)))
plt.xlabel(r'$m$')
plt.ylabel(r'Iterations')

filename_fig = "Variations_RK_N_1000_it"

plt.show()
fig.savefig("plots/seq/pdf/"+filename_fig+".pdf", bbox_inches='tight')
fig.savefig("plots/seq/png/"+filename_fig+".png", bbox_inches='tight')