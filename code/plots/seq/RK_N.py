import numpy as np
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import gca
import matplotlib.font_manager
import sys

# python3 plots/seq/RK_N.py

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
it_50 = [it[i] for i in indices]
time_50 = [time[i] for i in indices]
x_50 = [2000, 4000, 20000, 40000, 80000, 160000]

indices = (1, 6, 12, 20, 29, 39)
it_100 = [it[i] for i in indices]
time_100 = [time[i] for i in indices]
x_100 = [2000, 4000, 20000, 40000, 80000, 160000]

indices = (2, 7, 13, 21, 30, 40)
it_200 = [it[i] for i in indices]
time_200 = [time[i] for i in indices]
x_200 = [2000, 4000, 20000, 40000, 80000, 160000]

indices = (3, 8, 14, 22, 31, 41)
it_500 = [it[i] for i in indices]
time_500 = [time[i] for i in indices]
x_500 = [2000, 4000, 20000, 40000, 80000, 160000]

indices = (4, 9, 15, 23, 32, 42)
it_750 = [it[i] for i in indices]
time_750 = [time[i] for i in indices]
x_750 = [2000, 4000, 20000, 40000, 80000, 160000]

indices = (10, 16, 24, 33, 43)
it_1000 = [it[i] for i in indices]
time_1000 = [time[i] for i in indices]
x_1000 = [4000, 20000, 40000, 80000, 160000]

indices = (17, 25, 34, 44)
it_2000 = [it[i] for i in indices]
time_2000 = [time[i] for i in indices]
x_2000 = [20000, 40000, 80000, 160000]

indices = (18, 26, 35, 45)
it_4000 = [it[i] for i in indices]
time_4000 = [time[i] for i in indices]
x_4000 = [20000, 40000, 80000, 160000]

indices = (27, 36, 46)
it_10000 = [it[i] for i in indices]
time_10000 = [time[i] for i in indices]
x_10000 = [40000, 80000, 160000]

indices = (37, 47)
it_20000 = [it[i] for i in indices]
time_20000 = [time[i] for i in indices]
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

plot_title = r'RK - Overdetermined systems sampled from $N(\mu,\sigma)$'

fig = plt.figure(figsize=(10, 7))
plt.scatter(x_50, time_50, linewidth=1.5, color='yellow', label=r'$n = 50$')
plt.plot(x_50, time_50, linewidth=1.5, color='yellow')
plt.scatter(x_100, time_100, linewidth=1.5, color='orange', label=r'$n = 100$')
plt.plot(x_100, time_100, linewidth=1.5, color='orange')
plt.scatter(x_200, time_200, linewidth=1.5, color='red', label=r'$n = 200$')
plt.plot(x_200, time_200, linewidth=1.5, color='red')
plt.scatter(x_500, time_500, linewidth=1.5, color='magenta', label=r'$n = 500$')
plt.plot(x_500, time_500, linewidth=1.5, color='magenta')
plt.scatter(x_750, time_750, linewidth=1.5, color='purple', label=r'$n = 750$')
plt.plot(x_750, time_750, linewidth=1.5, color='purple')
plt.scatter(x_1000, time_1000, linewidth=1.5, color='brown', label=r'$n = 1000$')
plt.plot(x_1000, time_1000, linewidth=1.5, color='brown')
plt.scatter(x_2000, time_2000, linewidth=1.5, color='blue', label=r'$n = 2000$')
plt.plot(x_2000, time_2000, linewidth=1.5, color='blue')
plt.scatter(x_4000, time_4000, linewidth=1.5, color='black', label=r'$n = 4000$')
plt.plot(x_4000, time_4000, linewidth=1.5, color='black')
plt.scatter(x_10000, time_10000, linewidth=1.5, color='grey', label=r'$n = 10000$')
plt.plot(x_10000, time_10000, linewidth=1.5, color='grey')
plt.scatter(x_20000, time_20000, linewidth=1.5, color='cyan', label=r'$n = 20000$')
plt.plot(x_20000, time_20000, linewidth=1.5, color='cyan')
plt.grid()
plt.legend(loc='upper right')
# plt.title(plot_title)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$m$')
plt.ylabel(r'Total Time (s)')

filename_fig = "RK_N_time"

plt.show()
fig.savefig("plots/seq/pdf/"+filename_fig+".pdf", bbox_inches='tight')
fig.savefig("plots/seq/png/"+filename_fig+".png", bbox_inches='tight')

plot_title = r'RK - Overdetermined systems sampled from $N(\mu,\sigma)$'

fig = plt.figure(figsize=(10, 7))
plt.scatter(x_50, it_50, linewidth=1.5, color='yellow', label=r'$n = 50$')
plt.plot(x_50, it_50, linewidth=1.5, color='yellow')
plt.scatter(x_100, it_100, linewidth=1.5, color='orange', label=r'$n = 100$')
plt.plot(x_100, it_100, linewidth=1.5, color='orange')
plt.scatter(x_200, it_200, linewidth=1.5, color='red', label=r'$n = 200$')
plt.plot(x_200, it_200, linewidth=1.5, color='red')
plt.scatter(x_500, it_500, linewidth=1.5, color='magenta', label=r'$n = 500$')
plt.plot(x_500, it_500, linewidth=1.5, color='magenta')
plt.scatter(x_750, it_750, linewidth=1.5, color='purple', label=r'$n = 750$')
plt.plot(x_750, it_750, linewidth=1.5, color='purple')
plt.scatter(x_1000, it_1000, linewidth=1.5, color='brown', label=r'$n = 1000$')
plt.plot(x_1000, it_1000, linewidth=1.5, color='brown')
plt.scatter(x_2000, it_2000, linewidth=1.5, color='blue', label=r'$n = 2000$')
plt.plot(x_2000, it_2000, linewidth=1.5, color='blue')
plt.scatter(x_4000, it_4000, linewidth=1.5, color='black', label=r'$n = 4000$')
plt.plot(x_4000, it_4000, linewidth=1.5, color='black')
plt.scatter(x_10000, it_10000, linewidth=1.5, color='grey', label=r'$n = 10000$')
plt.plot(x_10000, it_10000, linewidth=1.5, color='grey')
plt.scatter(x_20000, it_20000, linewidth=1.5, color='cyan', label=r'$n = 20000$')
plt.plot(x_20000, it_20000, linewidth=1.5, color='cyan')
plt.grid()
plt.legend(loc='upper right')
# plt.title(plot_title)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$m$')
plt.ylabel(r'Iterations')

filename_fig = "RK_N_it"

plt.show()
fig.savefig("plots/seq/pdf/"+filename_fig+".pdf", bbox_inches='tight')
fig.savefig("plots/seq/png/"+filename_fig+".png", bbox_inches='tight')