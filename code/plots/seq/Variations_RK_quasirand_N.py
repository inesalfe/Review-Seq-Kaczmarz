import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import gca
import matplotlib.font_manager
import sys

# python3 plots/seq/Variations_RK_quasirand_N.py

filename = "outputs/seq/RK.txt";

with open(filename) as f:
    lines = f.read().splitlines()

file_size = len(lines)

it = []
time = []
for i in range(file_size):
    it.append(int(lines[i].split()[3]))
    time.append(float(lines[i].split()[2]))

indices = (0, 5, 11, 19, 28, 38)
RK_it_50 = [it[i] for i in indices]
RK_time_50 = [time[i] for i in indices]
indices = (1, 6, 12, 20, 29, 39)
RK_it_100 = [it[i] for i in indices]
RK_time_100 = [time[i] for i in indices]
indices = (2, 7, 13, 21, 30, 40)
RK_it_200 = [it[i] for i in indices]
RK_time_200 = [time[i] for i in indices]
indices = (3, 8, 14, 22, 31, 41)
RK_it_500 = [it[i] for i in indices]
RK_time_500 = [time[i] for i in indices]
indices = (4, 9, 15, 23, 32, 42)
RK_it_750 = [it[i] for i in indices]
RK_time_750 = [time[i] for i in indices]
indices = (10, 16, 24, 33, 43)
RK_it_1000 = [it[i] for i in indices]
RK_time_1000 = [time[i] for i in indices]
indices = (17, 25, 34, 44)
RK_it_2000 = [it[i] for i in indices]
RK_time_2000 = [time[i] for i in indices]
indices = (18, 26, 35, 45)
RK_it_4000 = [it[i] for i in indices]
RK_time_4000 = [time[i] for i in indices]
indices = (27, 36, 46)
RK_it_10000 = [it[i] for i in indices]
RK_time_10000 = [time[i] for i in indices]
indices = (37, 47)
RK_it_20000 = [it[i] for i in indices]
RK_time_20000 = [time[i] for i in indices]

filename = "outputs/seq/RK_quasirand_halton.txt";

with open(filename) as f:
    lines = f.read().splitlines()

file_size = len(lines)

it = []
time = []
for i in range(file_size):
    it.append(int(lines[i].split()[3]))
    time.append(float(lines[i].split()[2]))

indices = (0, 5, 11, 19, 28, 38)
RK_quasirand_halton_it_50 = [it[i] for i in indices]
RK_quasirand_halton_time_50 = [time[i] for i in indices]
indices = (1, 6, 12, 20, 29, 39)
RK_quasirand_halton_it_100 = [it[i] for i in indices]
RK_quasirand_halton_time_100 = [time[i] for i in indices]
indices = (2, 7, 13, 21, 30, 40)
RK_quasirand_halton_it_200 = [it[i] for i in indices]
RK_quasirand_halton_time_200 = [time[i] for i in indices]
indices = (3, 8, 14, 22, 31, 41)
RK_quasirand_halton_it_500 = [it[i] for i in indices]
RK_quasirand_halton_time_500 = [time[i] for i in indices]
indices = (4, 9, 15, 23, 32, 42)
RK_quasirand_halton_it_750 = [it[i] for i in indices]
RK_quasirand_halton_time_750 = [time[i] for i in indices]
indices = (10, 16, 24, 33, 43)
RK_quasirand_halton_it_1000 = [it[i] for i in indices]
RK_quasirand_halton_time_1000 = [time[i] for i in indices]
indices = (17, 25, 34, 44)
RK_quasirand_halton_it_2000 = [it[i] for i in indices]
RK_quasirand_halton_time_2000 = [time[i] for i in indices]
indices = (18, 26, 35, 45)
RK_quasirand_halton_it_4000 = [it[i] for i in indices]
RK_quasirand_halton_time_4000 = [time[i] for i in indices]
indices = (27, 36, 46)
RK_quasirand_halton_it_10000 = [it[i] for i in indices]
RK_quasirand_halton_time_10000 = [time[i] for i in indices]
indices = (37, 47)
RK_quasirand_halton_it_20000 = [it[i] for i in indices]
RK_quasirand_halton_time_20000 = [time[i] for i in indices]

filename = "outputs/seq/RK_quasirand_sobol.txt";

with open(filename) as f:
    lines = f.read().splitlines()

file_size = len(lines)

it = []
time = []
for i in range(file_size):
    it.append(int(lines[i].split()[3]))
    time.append(float(lines[i].split()[2]))

indices = (0, 5, 11, 19, 28, 38)
RK_quasirand_sobol_it_50 = [it[i] for i in indices]
RK_quasirand_sobol_time_50 = [time[i] for i in indices]
indices = (1, 6, 12, 20, 29, 39)
RK_quasirand_sobol_it_100 = [it[i] for i in indices]
RK_quasirand_sobol_time_100 = [time[i] for i in indices]
indices = (2, 7, 13, 21, 30, 40)
RK_quasirand_sobol_it_200 = [it[i] for i in indices]
RK_quasirand_sobol_time_200 = [time[i] for i in indices]
indices = (3, 8, 14, 22, 31, 41)
RK_quasirand_sobol_it_500 = [it[i] for i in indices]
RK_quasirand_sobol_time_500 = [time[i] for i in indices]
indices = (4, 9, 15, 23, 32, 42)
RK_quasirand_sobol_it_750 = [it[i] for i in indices]
RK_quasirand_sobol_time_750 = [time[i] for i in indices]
indices = (10, 16, 24, 33, 43)
RK_quasirand_sobol_it_1000 = [it[i] for i in indices]
RK_quasirand_sobol_time_1000 = [time[i] for i in indices]
indices = (17, 25, 34, 44)
RK_quasirand_sobol_it_2000 = [it[i] for i in indices]
RK_quasirand_sobol_time_2000 = [time[i] for i in indices]
indices = (18, 26, 35, 45)
RK_quasirand_sobol_it_4000 = [it[i] for i in indices]
RK_quasirand_sobol_time_4000 = [time[i] for i in indices]
indices = (27, 36, 46)
RK_quasirand_sobol_it_10000 = [it[i] for i in indices]
RK_quasirand_sobol_time_10000 = [time[i] for i in indices]
indices = (37, 47)
RK_quasirand_sobol_it_20000 = [it[i] for i in indices]
RK_quasirand_sobol_time_20000 = [time[i] for i in indices]

filename = "outputs/seq/RK_norep_rand_noshuffle.txt";

with open(filename) as f:
    lines = f.read().splitlines()

file_size = len(lines)

it = []
time = []
for i in range(file_size):
    it.append(int(lines[i].split()[3]))
    time.append(float(lines[i].split()[2]))

indices = (0, 5, 11, 19, 28, 38)
RK_norep_rand_noshuffle_it_50 = [it[i] for i in indices]
RK_norep_rand_noshuffle_time_50 = [time[i] for i in indices]
indices = (1, 6, 12, 20, 29, 39)
RK_norep_rand_noshuffle_it_100 = [it[i] for i in indices]
RK_norep_rand_noshuffle_time_100 = [time[i] for i in indices]
indices = (2, 7, 13, 21, 30, 40)
RK_norep_rand_noshuffle_it_200 = [it[i] for i in indices]
RK_norep_rand_noshuffle_time_200 = [time[i] for i in indices]
indices = (3, 8, 14, 22, 31, 41)
RK_norep_rand_noshuffle_it_500 = [it[i] for i in indices]
RK_norep_rand_noshuffle_time_500 = [time[i] for i in indices]
indices = (4, 9, 15, 23, 32, 42)
RK_norep_rand_noshuffle_it_750 = [it[i] for i in indices]
RK_norep_rand_noshuffle_time_750 = [time[i] for i in indices]
indices = (10, 16, 24, 33, 43)
RK_norep_rand_noshuffle_it_1000 = [it[i] for i in indices]
RK_norep_rand_noshuffle_time_1000 = [time[i] for i in indices]
indices = (17, 25, 34, 44)
RK_norep_rand_noshuffle_it_2000 = [it[i] for i in indices]
RK_norep_rand_noshuffle_time_2000 = [time[i] for i in indices]
indices = (18, 26, 35, 45)
RK_norep_rand_noshuffle_it_4000 = [it[i] for i in indices]
RK_norep_rand_noshuffle_time_4000 = [time[i] for i in indices]
indices = (27, 36, 46)
RK_norep_rand_noshuffle_it_10000 = [it[i] for i in indices]
RK_norep_rand_noshuffle_time_10000 = [time[i] for i in indices]
indices = (37, 47)
RK_norep_rand_noshuffle_it_20000 = [it[i] for i in indices]
RK_norep_rand_noshuffle_time_20000 = [time[i] for i in indices]

x_50 = [2000, 4000, 20000, 40000, 80000, 160000]
x_100 = [2000, 4000, 20000, 40000, 80000, 160000]
x_200 = [2000, 4000, 20000, 40000, 80000, 160000]
x_500 = [2000, 4000, 20000, 40000, 80000, 160000]
x_750 = [2000, 4000, 20000, 40000, 80000, 160000]
x_1000 = [4000, 20000, 40000, 80000, 160000]
x_2000 = [20000, 40000, 80000, 160000]
x_4000 = [20000, 40000, 80000, 160000]
x_10000 = [40000, 80000, 160000]
x_20000 = [80000, 160000]

import matplotlib.pylab as pylab
params = {'legend.fontsize': 'xx-large',
         'axes.labelsize': 'xx-large',
         'axes.titlesize': 'xx-large',
         'xtick.labelsize': 'xx-large',
         'ytick.labelsize': 'xx-large'}
pylab.rcParams.update(params)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plot_title = r'Overdetermined systems sampled from $N(\mu,\sigma)$ with $n=50$'

fig = plt.figure(figsize=(10, 7))
plt.scatter(x_50, RK_time_50, linewidth=1.5, color='green', label=r'RK')
plt.plot(x_50, RK_time_50, linewidth=1.5, color='green')
plt.scatter(x_50, RK_quasirand_halton_time_50, linewidth=1.5, color='orange', label=r'SRK - Halton')
plt.plot(x_50, RK_quasirand_halton_time_50, linewidth=1.5, color='orange')
plt.scatter(x_50, RK_quasirand_sobol_time_50, linewidth=1.5, color='blue', label=r'SRK - Sobol')
plt.plot(x_50, RK_quasirand_sobol_time_50, linewidth=1.5, color='blue')
plt.scatter(x_50, RK_norep_rand_noshuffle_time_50, linewidth=1.5, color='red', label=r'SRKWOR')
plt.plot(x_50, RK_norep_rand_noshuffle_time_50, linewidth=1.5, color='red')
plt.grid()
plt.legend(loc='upper right')
# plt.title(plot_title)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$m$')
plt.ylabel(r'Total Time (s)')

filename_fig = "Variations_RK_quasirand_N_time_50"

# plt.show()
fig.savefig("plots/seq/pdf/"+filename_fig+".pdf", bbox_inches='tight')
fig.savefig("plots/seq/png/"+filename_fig+".png", bbox_inches='tight')
plt.close()

plot_title = r'Overdetermined systems sampled from $N(\mu,\sigma)$ with $n=50$'

fig = plt.figure(figsize=(10, 7))
plt.scatter(x_50, RK_it_50, linewidth=1.5, color='green', label=r'RK')
plt.plot(x_50, RK_it_50, linewidth=1.5, color='green')
plt.scatter(x_50, RK_quasirand_halton_it_50, linewidth=1.5, color='orange', label=r'SRK - Halton')
plt.plot(x_50, RK_quasirand_halton_it_50, linewidth=1.5, color='orange')
plt.scatter(x_50, RK_quasirand_sobol_it_50, linewidth=1.5, color='blue', label=r'SRK - Sobol')
plt.plot(x_50, RK_quasirand_sobol_it_50, linewidth=1.5, color='blue')
plt.scatter(x_50, RK_norep_rand_noshuffle_it_50, linewidth=1.5, color='red', label=r'SRKWOR')
plt.plot(x_50, RK_norep_rand_noshuffle_it_50, linewidth=1.5, color='red')
plt.grid()
plt.legend(loc='upper right')
# plt.title(plot_title)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$m$')
plt.ylabel(r'Iterations')

filename_fig = "Variations_RK_quasirand_N_it_50"

# plt.show()
fig.savefig("plots/seq/pdf/"+filename_fig+".pdf", bbox_inches='tight')
fig.savefig("plots/seq/png/"+filename_fig+".png", bbox_inches='tight')
plt.close()

plot_title = r'Overdetermined systems sampled from $N(\mu,\sigma)$ with $n=100$'

fig = plt.figure(figsize=(10, 7))
plt.scatter(x_100, RK_time_100, linewidth=1.5, color='green', label=r'RK')
plt.plot(x_100, RK_time_100, linewidth=1.5, color='green')
plt.scatter(x_100, RK_quasirand_halton_time_100, linewidth=1.5, color='orange', label=r'SRK - Halton')
plt.plot(x_100, RK_quasirand_halton_time_100, linewidth=1.5, color='orange')
plt.scatter(x_100, RK_quasirand_sobol_time_100, linewidth=1.5, color='blue', label=r'SRK - Sobol')
plt.plot(x_100, RK_quasirand_sobol_time_100, linewidth=1.5, color='blue')
plt.scatter(x_100, RK_norep_rand_noshuffle_time_100, linewidth=1.5, color='red', label=r'SRKWOR')
plt.plot(x_100, RK_norep_rand_noshuffle_time_100, linewidth=1.5, color='red')
plt.grid()
plt.legend(loc='upper right')
# plt.title(plot_title)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$m$')
plt.ylabel(r'Total Time (s)')

filename_fig = "Variations_RK_quasirand_N_time_100"

# plt.show()
fig.savefig("plots/seq/pdf/"+filename_fig+".pdf", bbox_inches='tight')
fig.savefig("plots/seq/png/"+filename_fig+".png", bbox_inches='tight')
plt.close()

plot_title = r'Overdetermined systems sampled from $N(\mu,\sigma)$ with $n=100$'

fig = plt.figure(figsize=(10, 7))
plt.scatter(x_100, RK_it_100, linewidth=1.5, color='green', label=r'RK')
plt.plot(x_100, RK_it_100, linewidth=1.5, color='green')
plt.scatter(x_100, RK_quasirand_halton_it_100, linewidth=1.5, color='orange', label=r'SRK - Halton')
plt.plot(x_100, RK_quasirand_halton_it_100, linewidth=1.5, color='orange')
plt.scatter(x_100, RK_quasirand_sobol_it_100, linewidth=1.5, color='blue', label=r'SRK - Sobol')
plt.plot(x_100, RK_quasirand_sobol_it_100, linewidth=1.5, color='blue')
plt.scatter(x_100, RK_norep_rand_noshuffle_it_100, linewidth=1.5, color='red', label=r'SRKWOR')
plt.plot(x_100, RK_norep_rand_noshuffle_it_100, linewidth=1.5, color='red')
plt.grid()
plt.legend(loc='upper right')
# plt.title(plot_title)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$m$')
plt.ylabel(r'Iterations')

filename_fig = "Variations_RK_quasirand_N_it_100"

# plt.show()
fig.savefig("plots/seq/pdf/"+filename_fig+".pdf", bbox_inches='tight')
fig.savefig("plots/seq/png/"+filename_fig+".png", bbox_inches='tight')
plt.close()

plot_title = r'Overdetermined systems sampled from $N(\mu,\sigma)$ with $n=200$'

fig = plt.figure(figsize=(10, 7))
plt.scatter(x_200, RK_time_200, linewidth=1.5, color='green', label=r'RK')
plt.plot(x_200, RK_time_200, linewidth=1.5, color='green')
plt.scatter(x_200, RK_quasirand_halton_time_200, linewidth=1.5, color='orange', label=r'SRK - Halton')
plt.plot(x_200, RK_quasirand_halton_time_200, linewidth=1.5, color='orange')
plt.scatter(x_200, RK_quasirand_sobol_time_200, linewidth=1.5, color='blue', label=r'SRK - Sobol')
plt.plot(x_200, RK_quasirand_sobol_time_200, linewidth=1.5, color='blue')
plt.scatter(x_200, RK_norep_rand_noshuffle_time_200, linewidth=1.5, color='red', label=r'SRKWOR')
plt.plot(x_200, RK_norep_rand_noshuffle_time_200, linewidth=1.5, color='red')
plt.grid()
plt.legend(loc='upper right')
# plt.title(plot_title)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$m$')
plt.ylabel(r'Total Time (s)')

filename_fig = "Variations_RK_quasirand_N_time_200"

# plt.show()
fig.savefig("plots/seq/pdf/"+filename_fig+".pdf", bbox_inches='tight')
fig.savefig("plots/seq/png/"+filename_fig+".png", bbox_inches='tight')
plt.close()

plot_title = r'Overdetermined systems sampled from $N(\mu,\sigma)$ with $n=200$'

fig = plt.figure(figsize=(10, 7))
plt.scatter(x_200, RK_it_200, linewidth=1.5, color='green', label=r'RK')
plt.plot(x_200, RK_it_200, linewidth=1.5, color='green')
plt.scatter(x_200, RK_quasirand_halton_it_200, linewidth=1.5, color='orange', label=r'SRK - Halton')
plt.plot(x_200, RK_quasirand_halton_it_200, linewidth=1.5, color='orange')
plt.scatter(x_200, RK_quasirand_sobol_it_200, linewidth=1.5, color='blue', label=r'SRK - Sobol')
plt.plot(x_200, RK_quasirand_sobol_it_200, linewidth=1.5, color='blue')
plt.scatter(x_200, RK_norep_rand_noshuffle_it_200, linewidth=1.5, color='red', label=r'SRKWOR')
plt.plot(x_200, RK_norep_rand_noshuffle_it_200, linewidth=1.5, color='red')
plt.grid()
plt.legend(loc='upper right')
# plt.title(plot_title)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$m$')
plt.ylabel(r'Iterations')

filename_fig = "Variations_RK_quasirand_N_it_200"

# plt.show()
fig.savefig("plots/seq/pdf/"+filename_fig+".pdf", bbox_inches='tight')
fig.savefig("plots/seq/png/"+filename_fig+".png", bbox_inches='tight')
plt.close()

plot_title = r'Overdetermined systems sampled from $N(\mu,\sigma)$ with $n=500$'

fig = plt.figure(figsize=(10, 7))
plt.scatter(x_500, RK_time_500, linewidth=1.5, color='green', label=r'RK')
plt.plot(x_500, RK_time_500, linewidth=1.5, color='green')
plt.scatter(x_500, RK_quasirand_halton_time_500, linewidth=1.5, color='orange', label=r'SRK - Halton')
plt.plot(x_500, RK_quasirand_halton_time_500, linewidth=1.5, color='orange')
plt.scatter(x_500, RK_quasirand_sobol_time_500, linewidth=1.5, color='blue', label=r'SRK - Sobol')
plt.plot(x_500, RK_quasirand_sobol_time_500, linewidth=1.5, color='blue')
plt.scatter(x_500, RK_norep_rand_noshuffle_time_500, linewidth=1.5, color='red', label=r'SRKWOR')
plt.plot(x_500, RK_norep_rand_noshuffle_time_500, linewidth=1.5, color='red')
plt.grid()
plt.legend(loc='upper right')
# plt.title(plot_title)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$m$')
plt.ylabel(r'Total Time (s)')

filename_fig = "Variations_RK_quasirand_N_time_500"

# plt.show()
fig.savefig("plots/seq/pdf/"+filename_fig+".pdf", bbox_inches='tight')
fig.savefig("plots/seq/png/"+filename_fig+".png", bbox_inches='tight')
plt.close()

plot_title = r'Overdetermined systems sampled from $N(\mu,\sigma)$ with $n=500$'

fig = plt.figure(figsize=(10, 7))
plt.scatter(x_500, RK_it_500, linewidth=1.5, color='green', label=r'RK')
plt.plot(x_500, RK_it_500, linewidth=1.5, color='green')
plt.scatter(x_500, RK_quasirand_halton_it_500, linewidth=1.5, color='orange', label=r'SRK - Halton')
plt.plot(x_500, RK_quasirand_halton_it_500, linewidth=1.5, color='orange')
plt.scatter(x_500, RK_quasirand_sobol_it_500, linewidth=1.5, color='blue', label=r'SRK - Sobol')
plt.plot(x_500, RK_quasirand_sobol_it_500, linewidth=1.5, color='blue')
plt.scatter(x_500, RK_norep_rand_noshuffle_it_500, linewidth=1.5, color='red', label=r'SRKWOR')
plt.plot(x_500, RK_norep_rand_noshuffle_it_500, linewidth=1.5, color='red')
plt.grid()
plt.legend(loc='upper right')
# plt.title(plot_title)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$m$')
plt.ylabel(r'Iterations')

filename_fig = "Variations_RK_quasirand_N_it_500"

# plt.show()
fig.savefig("plots/seq/pdf/"+filename_fig+".pdf", bbox_inches='tight')
fig.savefig("plots/seq/png/"+filename_fig+".png", bbox_inches='tight')
plt.close()

plot_title = r'Overdetermined systems sampled from $N(\mu,\sigma)$ with $n=750$'

fig = plt.figure(figsize=(10, 7))
plt.scatter(x_750, RK_time_750, linewidth=1.5, color='green', label=r'RK')
plt.plot(x_750, RK_time_750, linewidth=1.5, color='green')
plt.scatter(x_750, RK_quasirand_halton_time_750, linewidth=1.5, color='orange', label=r'SRK - Halton')
plt.plot(x_750, RK_quasirand_halton_time_750, linewidth=1.5, color='orange')
plt.scatter(x_750, RK_quasirand_sobol_time_750, linewidth=1.5, color='blue', label=r'SRK - Sobol')
plt.plot(x_750, RK_quasirand_sobol_time_750, linewidth=1.5, color='blue')
plt.scatter(x_750, RK_norep_rand_noshuffle_time_750, linewidth=1.5, color='red', label=r'SRKWOR')
plt.plot(x_750, RK_norep_rand_noshuffle_time_750, linewidth=1.5, color='red')
plt.grid()
plt.legend(loc='upper right')
# plt.title(plot_title)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$m$')
plt.ylabel(r'Total Time (s)')

filename_fig = "Variations_RK_quasirand_N_time_750"

# plt.show()
fig.savefig("plots/seq/pdf/"+filename_fig+".pdf", bbox_inches='tight')
fig.savefig("plots/seq/png/"+filename_fig+".png", bbox_inches='tight')
plt.close()

plot_title = r'Overdetermined systems sampled from $N(\mu,\sigma)$ with $n=750$'

fig = plt.figure(figsize=(10, 7))
plt.scatter(x_750, RK_it_750, linewidth=1.5, color='green', label=r'RK')
plt.plot(x_750, RK_it_750, linewidth=1.5, color='green')
plt.scatter(x_750, RK_quasirand_halton_it_750, linewidth=1.5, color='orange', label=r'SRK - Halton')
plt.plot(x_750, RK_quasirand_halton_it_750, linewidth=1.5, color='orange')
plt.scatter(x_750, RK_quasirand_sobol_it_750, linewidth=1.5, color='blue', label=r'SRK - Sobol')
plt.plot(x_750, RK_quasirand_sobol_it_750, linewidth=1.5, color='blue')
plt.scatter(x_750, RK_norep_rand_noshuffle_it_750, linewidth=1.5, color='red', label=r'SRKWOR')
plt.plot(x_750, RK_norep_rand_noshuffle_it_750, linewidth=1.5, color='red')
plt.grid()
plt.legend(loc='upper right')
# plt.title(plot_title)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$m$')
plt.ylabel(r'Iterations')

filename_fig = "Variations_RK_quasirand_N_it_750"

# plt.show()
fig.savefig("plots/seq/pdf/"+filename_fig+".pdf", bbox_inches='tight')
fig.savefig("plots/seq/png/"+filename_fig+".png", bbox_inches='tight')
plt.close()

plot_title = r'Overdetermined systems sampled from $N(\mu,\sigma)$ with $n=1000$'

fig = plt.figure(figsize=(10, 7))
plt.scatter(x_1000, RK_time_1000, linewidth=1.5, color='green', label=r'RK')
plt.plot(x_1000, RK_time_1000, linewidth=1.5, color='green')
plt.scatter(x_1000, RK_quasirand_halton_time_1000, linewidth=1.5, color='orange', label=r'SRK - Halton')
plt.plot(x_1000, RK_quasirand_halton_time_1000, linewidth=1.5, color='orange')
plt.scatter(x_1000, RK_quasirand_sobol_time_1000, linewidth=1.5, color='blue', label=r'SRK - Sobol')
plt.plot(x_1000, RK_quasirand_sobol_time_1000, linewidth=1.5, color='blue')
plt.scatter(x_1000, RK_norep_rand_noshuffle_time_1000, linewidth=1.5, color='red', label=r'SRKWOR')
plt.plot(x_1000, RK_norep_rand_noshuffle_time_1000, linewidth=1.5, color='red')
plt.grid()
plt.legend(loc='upper right')
# plt.title(plot_title)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$m$')
plt.ylabel(r'Total Time (s)')

filename_fig = "Variations_RK_quasirand_N_time_1000"

# plt.show()
fig.savefig("plots/seq/pdf/"+filename_fig+".pdf", bbox_inches='tight')
fig.savefig("plots/seq/png/"+filename_fig+".png", bbox_inches='tight')
plt.close()

plot_title = r'Overdetermined systems sampled from $N(\mu,\sigma)$ with $n=1000$'

fig = plt.figure(figsize=(10, 7))
plt.scatter(x_1000, RK_it_1000, linewidth=1.5, color='green', label=r'RK')
plt.plot(x_1000, RK_it_1000, linewidth=1.5, color='green')
plt.scatter(x_1000, RK_quasirand_halton_it_1000, linewidth=1.5, color='orange', label=r'SRK - Halton')
plt.plot(x_1000, RK_quasirand_halton_it_1000, linewidth=1.5, color='orange')
plt.scatter(x_1000, RK_quasirand_sobol_it_1000, linewidth=1.5, color='blue', label=r'SRK - Sobol')
plt.plot(x_1000, RK_quasirand_sobol_it_1000, linewidth=1.5, color='blue')
plt.scatter(x_1000, RK_norep_rand_noshuffle_it_1000, linewidth=1.5, color='red', label=r'SRKWOR')
plt.plot(x_1000, RK_norep_rand_noshuffle_it_1000, linewidth=1.5, color='red')
plt.grid()
plt.legend(loc='upper right')
# plt.title(plot_title)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$m$')
plt.ylabel(r'Iterations')

filename_fig = "Variations_RK_quasirand_N_it_1000"

# plt.show()
fig.savefig("plots/seq/pdf/"+filename_fig+".pdf", bbox_inches='tight')
fig.savefig("plots/seq/png/"+filename_fig+".png", bbox_inches='tight')
plt.close()

plot_title = r'Overdetermined systems sampled from $N(\mu,\sigma)$ with $n=2000$'

fig = plt.figure(figsize=(10, 7))
plt.scatter(x_2000, RK_time_2000, linewidth=1.5, color='green', label=r'RK')
plt.plot(x_2000, RK_time_2000, linewidth=1.5, color='green')
plt.scatter(x_2000, RK_quasirand_halton_time_2000, linewidth=1.5, color='orange', label=r'SRK - Halton')
plt.plot(x_2000, RK_quasirand_halton_time_2000, linewidth=1.5, color='orange')
plt.scatter(x_2000, RK_quasirand_sobol_time_2000, linewidth=1.5, color='blue', label=r'SRK - Sobol')
plt.plot(x_2000, RK_quasirand_sobol_time_2000, linewidth=1.5, color='blue')
plt.scatter(x_2000, RK_norep_rand_noshuffle_time_2000, linewidth=1.5, color='red', label=r'SRKWOR')
plt.plot(x_2000, RK_norep_rand_noshuffle_time_2000, linewidth=1.5, color='red')
plt.grid()
plt.legend(loc='upper right')
# plt.title(plot_title)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$m$')
plt.ylabel(r'Total Time (s)')

filename_fig = "Variations_RK_quasirand_N_time_2000"

# plt.show()
fig.savefig("plots/seq/pdf/"+filename_fig+".pdf", bbox_inches='tight')
fig.savefig("plots/seq/png/"+filename_fig+".png", bbox_inches='tight')
plt.close()

plot_title = r'Overdetermined systems sampled from $N(\mu,\sigma)$ with $n=2000$'

fig = plt.figure(figsize=(10, 7))
plt.scatter(x_2000, RK_it_2000, linewidth=1.5, color='green', label=r'RK')
plt.plot(x_2000, RK_it_2000, linewidth=1.5, color='green')
plt.scatter(x_2000, RK_quasirand_halton_it_2000, linewidth=1.5, color='orange', label=r'SRK - Halton')
plt.plot(x_2000, RK_quasirand_halton_it_2000, linewidth=1.5, color='orange')
plt.scatter(x_2000, RK_quasirand_sobol_it_2000, linewidth=1.5, color='blue', label=r'SRK - Sobol')
plt.plot(x_2000, RK_quasirand_sobol_it_2000, linewidth=1.5, color='blue')
plt.scatter(x_2000, RK_norep_rand_noshuffle_it_2000, linewidth=1.5, color='red', label=r'SRKWOR')
plt.plot(x_2000, RK_norep_rand_noshuffle_it_2000, linewidth=1.5, color='red')
plt.grid()
plt.legend(loc='upper right')
# plt.title(plot_title)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$m$')
plt.ylabel(r'Iterations')

filename_fig = "Variations_RK_quasirand_N_it_2000"

# plt.show()
fig.savefig("plots/seq/pdf/"+filename_fig+".pdf", bbox_inches='tight')
fig.savefig("plots/seq/png/"+filename_fig+".png", bbox_inches='tight')
plt.close()

plot_title = r'Overdetermined systems sampled from $N(\mu,\sigma)$ with $n=4000$'

fig = plt.figure(figsize=(10, 7))
plt.scatter(x_4000, RK_time_4000, linewidth=1.5, color='green', label=r'RK')
plt.plot(x_4000, RK_time_4000, linewidth=1.5, color='green')
plt.scatter(x_4000, RK_quasirand_halton_time_4000, linewidth=1.5, color='orange', label=r'SRK - Halton')
plt.plot(x_4000, RK_quasirand_halton_time_4000, linewidth=1.5, color='orange')
plt.scatter(x_4000, RK_quasirand_sobol_time_4000, linewidth=1.5, color='blue', label=r'SRK - Sobol')
plt.plot(x_4000, RK_quasirand_sobol_time_4000, linewidth=1.5, color='blue')
plt.scatter(x_4000, RK_norep_rand_noshuffle_time_4000, linewidth=1.5, color='red', label=r'SRKWOR')
plt.plot(x_4000, RK_norep_rand_noshuffle_time_4000, linewidth=1.5, color='red')
plt.grid()
plt.legend(loc='upper right')
# plt.title(plot_title)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$m$')
plt.ylabel(r'Total Time (s)')

filename_fig = "Variations_RK_quasirand_N_time_4000"

# plt.show()
fig.savefig("plots/seq/pdf/"+filename_fig+".pdf", bbox_inches='tight')
fig.savefig("plots/seq/png/"+filename_fig+".png", bbox_inches='tight')
plt.close()

plot_title = r'Overdetermined systems sampled from $N(\mu,\sigma)$ with $n=4000$'

fig = plt.figure(figsize=(10, 7))
plt.scatter(x_4000, RK_it_4000, linewidth=1.5, color='green', label=r'RK')
plt.plot(x_4000, RK_it_4000, linewidth=1.5, color='green')
plt.scatter(x_4000, RK_quasirand_halton_it_4000, linewidth=1.5, color='orange', label=r'SRK - Halton')
plt.plot(x_4000, RK_quasirand_halton_it_4000, linewidth=1.5, color='orange')
plt.scatter(x_4000, RK_quasirand_sobol_it_4000, linewidth=1.5, color='blue', label=r'SRK - Sobol')
plt.plot(x_4000, RK_quasirand_sobol_it_4000, linewidth=1.5, color='blue')
plt.scatter(x_4000, RK_norep_rand_noshuffle_it_4000, linewidth=1.5, color='red', label=r'SRKWOR')
plt.plot(x_4000, RK_norep_rand_noshuffle_it_4000, linewidth=1.5, color='red')
plt.grid()
plt.legend(loc='upper right')
# plt.title(plot_title)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$m$')
plt.ylabel(r'Iterations')

filename_fig = "Variations_RK_quasirand_N_it_4000"

# plt.show()
fig.savefig("plots/seq/pdf/"+filename_fig+".pdf", bbox_inches='tight')
fig.savefig("plots/seq/png/"+filename_fig+".png", bbox_inches='tight')
plt.close()

plot_title = r'Overdetermined systems sampled from $N(\mu,\sigma)$ with $n=10000$'

fig = plt.figure(figsize=(10, 7))
plt.scatter(x_10000, RK_time_10000, linewidth=1.5, color='green', label=r'RK')
plt.plot(x_10000, RK_time_10000, linewidth=1.5, color='green')
plt.scatter(x_10000, RK_quasirand_halton_time_10000, linewidth=1.5, color='orange', label=r'SRK - Halton')
plt.plot(x_10000, RK_quasirand_halton_time_10000, linewidth=1.5, color='orange')
plt.scatter(x_10000, RK_quasirand_sobol_time_10000, linewidth=1.5, color='blue', label=r'SRK - Sobol')
plt.plot(x_10000, RK_quasirand_sobol_time_10000, linewidth=1.5, color='blue')
plt.scatter(x_10000, RK_norep_rand_noshuffle_time_10000, linewidth=1.5, color='red', label=r'SRKWOR')
plt.plot(x_10000, RK_norep_rand_noshuffle_time_10000, linewidth=1.5, color='red')
plt.grid()
plt.legend(loc='upper right')
# plt.title(plot_title)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$m$')
plt.ylabel(r'Total Time (s)')

filename_fig = "Variations_RK_quasirand_N_time_10000"

# plt.show()
fig.savefig("plots/seq/pdf/"+filename_fig+".pdf", bbox_inches='tight')
fig.savefig("plots/seq/png/"+filename_fig+".png", bbox_inches='tight')
plt.close()

plot_title = r'Overdetermined systems sampled from $N(\mu,\sigma)$ with $n=10000$'

fig = plt.figure(figsize=(10, 7))
plt.scatter(x_10000, RK_it_10000, linewidth=1.5, color='green', label=r'RK')
plt.plot(x_10000, RK_it_10000, linewidth=1.5, color='green')
plt.scatter(x_10000, RK_quasirand_halton_it_10000, linewidth=1.5, color='orange', label=r'SRK - Halton')
plt.plot(x_10000, RK_quasirand_halton_it_10000, linewidth=1.5, color='orange')
plt.scatter(x_10000, RK_quasirand_sobol_it_10000, linewidth=1.5, color='blue', label=r'SRK - Sobol')
plt.plot(x_10000, RK_quasirand_sobol_it_10000, linewidth=1.5, color='blue')
plt.scatter(x_10000, RK_norep_rand_noshuffle_it_10000, linewidth=1.5, color='red', label=r'SRKWOR')
plt.plot(x_10000, RK_norep_rand_noshuffle_it_10000, linewidth=1.5, color='red')
plt.grid()
plt.legend(loc='upper right')
# plt.title(plot_title)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$m$')
plt.ylabel(r'Iterations')

filename_fig = "Variations_RK_quasirand_N_it_10000"

# plt.show()
fig.savefig("plots/seq/pdf/"+filename_fig+".pdf", bbox_inches='tight')
fig.savefig("plots/seq/png/"+filename_fig+".png", bbox_inches='tight')
plt.close()

plot_title = r'Overdetermined systems sampled from $N(\mu,\sigma)$ with $n=20000$'

fig = plt.figure(figsize=(10, 7))
plt.scatter(x_20000, RK_time_20000, linewidth=1.5, color='green', label=r'RK')
plt.plot(x_20000, RK_time_20000, linewidth=1.5, color='green')
plt.scatter(x_20000, RK_quasirand_halton_time_20000, linewidth=1.5, color='orange', label=r'SRK - Halton')
plt.plot(x_20000, RK_quasirand_halton_time_20000, linewidth=1.5, color='orange')
plt.scatter(x_20000, RK_quasirand_sobol_time_20000, linewidth=1.5, color='blue', label=r'SRK - Sobol')
plt.plot(x_20000, RK_quasirand_sobol_time_20000, linewidth=1.5, color='blue')
plt.scatter(x_20000, RK_norep_rand_noshuffle_time_20000, linewidth=1.5, color='red', label=r'SRKWOR')
plt.plot(x_20000, RK_norep_rand_noshuffle_time_20000, linewidth=1.5, color='red')
plt.grid()
plt.legend(loc='upper right')
# plt.title(plot_title)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$m$')
plt.ylabel(r'Total Time (s)')

filename_fig = "Variations_RK_quasirand_N_time_20000"

# plt.show()
fig.savefig("plots/seq/pdf/"+filename_fig+".pdf", bbox_inches='tight')
fig.savefig("plots/seq/png/"+filename_fig+".png", bbox_inches='tight')
plt.close()

plot_title = r'Overdetermined systems sampled from $N(\mu,\sigma)$ with $n=20000$'

fig = plt.figure(figsize=(10, 7))
plt.scatter(x_20000, RK_it_20000, linewidth=1.5, color='green', label=r'RK')
plt.plot(x_20000, RK_it_20000, linewidth=1.5, color='green')
plt.scatter(x_20000, RK_quasirand_halton_it_20000, linewidth=1.5, color='orange', label=r'SRK - Halton')
plt.plot(x_20000, RK_quasirand_halton_it_20000, linewidth=1.5, color='orange')
plt.scatter(x_20000, RK_quasirand_sobol_it_20000, linewidth=1.5, color='blue', label=r'SRK - Sobol')
plt.plot(x_20000, RK_quasirand_sobol_it_20000, linewidth=1.5, color='blue')
plt.scatter(x_20000, RK_norep_rand_noshuffle_it_20000, linewidth=1.5, color='red', label=r'SRKWOR')
plt.plot(x_20000, RK_norep_rand_noshuffle_it_20000, linewidth=1.5, color='red')
plt.grid()
plt.legend(loc='upper right')
# plt.title(plot_title)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$m$')
plt.ylabel(r'Iterations')

filename_fig = "Variations_RK_quasirand_N_it_20000"

# plt.show()
fig.savefig("plots/seq/pdf/"+filename_fig+".pdf", bbox_inches='tight')
fig.savefig("plots/seq/png/"+filename_fig+".png", bbox_inches='tight')
plt.close()