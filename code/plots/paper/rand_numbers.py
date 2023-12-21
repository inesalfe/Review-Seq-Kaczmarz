import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import gca
import matplotlib.font_manager
import sys

import matplotlib.pylab as pylab
params = {'legend.fontsize': 'xx-large',
         'axes.labelsize': 'xx-large',
         'axes.titlesize': 'xx-large',
         'xtick.labelsize': 'xx-large',
         'ytick.labelsize': 'xx-large'}
pylab.rcParams.update(params)

# python3 plots/paper/rand_numbers.py

filename = "outputs/seq/rand_rows.txt";

with open(filename) as f:
	lines = f.read().splitlines()

file_size = len(lines)

rand = []
for i in range(file_size):
	rand.append(int(lines[i].split()[0]))

filename = "outputs/seq/sobol_rows.txt";

with open(filename) as f:
	lines = f.read().splitlines()

file_size = len(lines)

sobol = []
for i in range(file_size):
	sobol.append(int(lines[i].split()[0]))

filename = "outputs/seq/halton_rows.txt";

with open(filename) as f:
	lines = f.read().splitlines()

file_size = len(lines)

halton = []
for i in range(file_size):
	halton.append(int(lines[i].split()[0]))

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

fig = plt.figure(figsize=(5, 2))

plt.hlines(1, 1, 1000, linewidth=2, color='red')
plt.eventplot(rand, linewidth=1, linelengths=0.5, orientation='horizontal', colors='blue')
plt.eventplot([1,1000], linewidth=2.5, linelengths=1, orientation='horizontal', color='black')
plt.axis('off')
plt.xlim([0, 1001])
plt.ylim([0, 2])

filename_fig = "rand_numbers"

plt.show()
fig.savefig("plots/paper/pdf/"+filename_fig+".pdf", bbox_inches='tight')
fig.savefig("plots/paper/png/"+filename_fig+".png", bbox_inches='tight')

fig = plt.figure(figsize=(5, 2))

plt.hlines(1, 1, 1000, linewidth=2, color='red')
plt.eventplot(sobol, linewidth=1, linelengths=0.5, orientation='horizontal', colors='blue')
plt.eventplot([1,1000], linewidth=2.5, linelengths=1, orientation='horizontal', color='black')
plt.axis('off')
plt.xlim([0, 1001])
plt.ylim([0, 2])

filename_fig = "sobol_numbers"

plt.show()
fig.savefig("plots/paper/pdf/"+filename_fig+".pdf", bbox_inches='tight')
fig.savefig("plots/paper/png/"+filename_fig+".png", bbox_inches='tight')

fig = plt.figure(figsize=(5, 2))

plt.hlines(1, 1, 1000, linewidth=2, color='red')
plt.eventplot(halton, linewidth=1, linelengths=0.5, orientation='horizontal', colors='blue')
plt.eventplot([1,1000], linewidth=2.5, linelengths=1, orientation='horizontal', color='black')
plt.axis('off')
plt.xlim([0, 1001])
plt.ylim([0, 2])

filename_fig = "halton_numbers"

plt.show()
fig.savefig("plots/paper/pdf/"+filename_fig+".pdf", bbox_inches='tight')
fig.savefig("plots/paper/png/"+filename_fig+".png", bbox_inches='tight')