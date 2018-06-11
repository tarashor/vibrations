# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import os
import platform

import utils

folder = "./results/shear/"
#folder = ""
#results2D_all_n = utils.load_results(folder+"rNs2D")
#results1D1_all_n = utils.load_results(folder+"rNs1D1")
results1D2_all_n = utils.load_results(folder+"rNs1D2")

x1D2 = []
y1D2 = []
for key, value in results1D2_all_n.items():
    x1D2.append(key)
    y1D2.append(value[0].freqHz()+0.5)
    
#print(x)
#print(y)
    
tex_path = '/usr/local/texlive/2017/bin/x86_64-darwin'
if (platform.system() == 'Windows'):
    tex_path = "C:\Program Files\MiKTeX 2.9\miktex/bin/x64"

os.environ["PATH"] += os.pathsep + tex_path

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

#plt.plot(x, y, 'o-', linewidth=2.0, markersize=8, markeredgewidth=2, markeredgecolor='r', markerfacecolor='None', label = "General 2D theory")
#plt.plot(x1D1, y1D1, 'x--', linewidth=2.0, markersize=8, markeredgewidth=2, markeredgecolor='b', markerfacecolor='None', label = "Mindlin-Reissner theory")
plt.plot(x1D2, y1D2, 'v-', linewidth=2.0, markersize=8, markeredgewidth=2, markeredgecolor='r', markerfacecolor='None', label = "Square 1D theory")
plt.xlabel(r"$k$ - shear factor parameter")
plt.ylabel(r"$\omega_{min}$, Hz")

plt.legend(loc='best')

plt.title(r"Shear influence")
#plt.title(r"Convergence $\frac{E}{E_3} = 1$")
#plt.title("Збіжність")

plt.grid()
plt.show()
