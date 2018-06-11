# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import os

import utils

folder = "./results/convergE3/"
#folder = ""
results2D_all_n = utils.load_results(folder+"rNs2D")
results1D1_all_n = utils.load_results(folder+"rNs1D1")
results1D2_all_n = utils.load_results(folder+"rNs1D2")

x = []
y = []
for key, value in results2D_all_n.items():
    x.append(key)
    y.append(value[0].freqHz())
    
x1D1 = []
y1D1 = []
for key, value in results1D1_all_n.items():
    x1D1.append(key)
    y1D1.append(value[0].freqHz())
    
x1D2 = []
y1D2 = []
for key, value in results1D2_all_n.items():
    x1D2.append(key)
    y1D2.append(value[0].freqHz()+0.5)
    
#print(x)
#print(y)

os.environ["PATH"] += os.pathsep + '/usr/local/texlive/2017/bin/x86_64-darwin'

plt.rc('text', usetex=True)
   
plt.rc('font', family='serif')

SMALL_SIZE = 24
MEDIUM_SIZE = 28
BIGGER_SIZE = 32

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

plt.plot(x, y, 'o-', linewidth=2.0, markersize=8, markeredgewidth=2, markeredgecolor='r', markerfacecolor='None', label = "General 2D theory")
plt.plot(x1D1, y1D1, 'x--', linewidth=2.0, markersize=8, markeredgewidth=2, markeredgecolor='b', markerfacecolor='None', label = "Mindlin-Reissner theory")
plt.plot(x1D2, y1D2, 'v:', linewidth=2.0, markersize=8, markeredgewidth=2, markeredgecolor='g', markerfacecolor='None', label = "Square 1D theory")
plt.xlabel(r"$N$ - elements along $\alpha_1$")
plt.ylabel(r"$\omega_{min}$, Hz")

plt.legend(loc='best')

plt.title("Convergence")
#plt.title("Збіжність")

plt.grid()
plt.show()
