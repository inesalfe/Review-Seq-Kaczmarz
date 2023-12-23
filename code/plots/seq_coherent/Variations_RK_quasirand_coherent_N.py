import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import gca
import matplotlib.font_manager
import sys

# python3 plots/seq_coherent/Variations_RK_quasirand_coherent_N.py

filename = "outputs/seq_coherent/RK_coherent.txt";

with open(filename) as f:
    lines = f.read().splitlines()

file_size = len(lines)

it = []
time = []
for i in range(5):
    it.append(int(lines[i].split()[3]))
    time.append(float(lines[i].split()[2]))

RK_it_1000 = it
RK_time_1000 = time

filename = "outputs/seq_coherent/RK_quasirand_halton_coherent.txt";

with open(filename) as f:
    lines = f.read().splitlines()

file_size = len(lines)

it = []
time = []
for i in range(5):
    it.append(int(lines[i].split()[3]))
    time.append(float(lines[i].split()[2]))

SRK_halton_it_1000 = it
SRK_halton_time_1000 = time

filename = "outputs/seq_coherent/RK_quasirand_sobol_coherent.txt";

with open(filename) as f:
    lines = f.read().splitlines()

file_size = len(lines)

it = []
time = []
for i in range(5):
    it.append(int(lines[i].split()[3]))
    time.append(float(lines[i].split()[2]))

SRK_sobol_it_1000 = it
SRK_sobol_time_1000 = time

filename = "outputs/seq_coherent/RK_rand_coherent.txt";

with open(filename) as f:
    lines = f.read().splitlines()

file_size = len(lines)

it = []
time = []
for i in range(5):
    it.append(int(lines[i].split()[3]))
    time.append(float(lines[i].split()[2]))

RK_rand_it_1000 = it
RK_rand_time_1000 = time

filename = "outputs/seq_coherent/RK_norep_rand_noshuffle_coherent.txt";

with open(filename) as f:
    lines = f.read().splitlines()

file_size = len(lines)

it = []
time = []
for i in range(5):
    it.append(int(lines[i].split()[3]))
    time.append(float(lines[i].split()[2]))

RK_norep_rand_noshuffle_it_1000 = it
RK_norep_rand_noshuffle_time_1000 = time

x_1000 = [4000, 20000, 40000, 80000, 160000]

import matplotlib.pylab as pylab
params = {'legend.fontsize': 'xx-large',
         'axes.labelsize': 'xx-large',
         'axes.titlesize': 'xx-large',
         'xtick.labelsize': 'xx-large',
         'ytick.labelsize': 'xx-large'}
pylab.rcParams.update(params)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plot_title = r'Overdetermined systems sampled from $N(\mu,\sigma)$ with $n=1000$'

fig = plt.figure(figsize=(10, 7))
plt.scatter(x_1000, RK_time_1000, linewidth=1.5, color='green', label=r'RK', marker='^', s=75)
plt.plot(x_1000, RK_time_1000, linewidth=1.5, color='green')
plt.scatter(x_1000, SRK_halton_time_1000, linewidth=1.5, color='orange', label=r'CK', marker='o', s=75)
plt.plot(x_1000, SRK_halton_time_1000, linewidth=1.5, color='orange')
plt.scatter(x_1000, RK_rand_time_1000, linewidth=1.5, color='blue', label=r'SRK', marker='s', s=75)
plt.plot(x_1000, RK_rand_time_1000, linewidth=1.5, color='blue')
plt.scatter(x_1000, RK_norep_rand_noshuffle_time_1000, linewidth=1.5, color='red', label=r'SRKWOR', marker='X', s=75)
plt.plot(x_1000, RK_norep_rand_noshuffle_time_1000, linewidth=1.5, color='red')
plt.scatter(x_1000, SRK_sobol_time_1000, linewidth=1.5, color='black', label=r'SRK-Sobol', marker='o', s=75)
plt.plot(x_1000, SRK_sobol_time_1000, linewidth=1.5, color='black')
plt.grid()
plt.legend(loc='upper right')
# plt.title(plot_title)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$m$')
plt.ylabel(r'Total Time (s)')

filename_fig = "Variations_RK_quasirand_N_time_coherent_1000"

# plt.show()
fig.savefig("plots/seq_coherent/pdf/"+filename_fig+".pdf", bbox_inches='tight')
fig.savefig("plots/seq_coherent/png/"+filename_fig+".png", bbox_inches='tight')
plt.close()

plot_title = r'Overdetermined systems sampled from $N(\mu,\sigma)$ with $n=1000$'

fig = plt.figure(figsize=(10, 7))
plt.scatter(x_1000, RK_it_1000, linewidth=1.5, color='green', label=r'RK', marker='^', s=75)
plt.plot(x_1000, RK_it_1000, linewidth=1.5, color='green')
plt.scatter(x_1000, SRK_halton_it_1000, linewidth=1.5, color='orange', label=r'SRK-Halton', marker='o', s=75)
plt.plot(x_1000, SRK_halton_it_1000, linewidth=1.5, color='orange')
plt.scatter(x_1000, RK_rand_it_1000, linewidth=1.5, color='blue', label=r'SRK', marker='s', s=75)
plt.plot(x_1000, RK_rand_it_1000, linewidth=1.5, color='blue')
plt.scatter(x_1000, RK_norep_rand_noshuffle_it_1000, linewidth=1.5, color='red', label=r'SRKWOR', marker='X', s=75)
plt.plot(x_1000, RK_norep_rand_noshuffle_it_1000, linewidth=1.5, color='red')
plt.scatter(x_1000, SRK_sobol_it_1000, linewidth=1.5, color='black', label=r'SRK-Sobol', marker='o', s=75)
plt.plot(x_1000, SRK_sobol_it_1000, linewidth=1.5, color='black')
plt.grid()
plt.legend(loc='upper right')
# plt.title(plot_title)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$m$')
plt.ylabel(r'Iterations')

filename_fig = "Variations_RK_quasirand_N_it_coherent_1000"

# plt.show()
fig.savefig("plots/seq_coherent/pdf/"+filename_fig+".pdf", bbox_inches='tight')
fig.savefig("plots/seq_coherent/png/"+filename_fig+".png", bbox_inches='tight')
plt.close()