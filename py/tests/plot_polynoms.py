import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.patches as patches
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

import os

os.environ["PATH"] += os.pathsep + "C:\Program Files\MiKTeX 2.9\miktex/bin/x64"
#'/usr/local/texlive/2017/bin/x86_64-darwin'

plt.rc('text', usetex=True)
   
plt.rc('font', family='serif')

SMALL_SIZE = 42
MEDIUM_SIZE = 42
BIGGER_SIZE = 42

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title



h = 1

a3 = np.linspace(-h/2, h/2, 200)
p0 = []
p1 = []
p2 = []
for a in a3:
    p0.append(0.5-a/h)
    p1.append(0.5+a/h)
    p2.append(1-(2*a/h)**2)

plt.plot(a3, p0, dashes=[6, 2], color="r", label=r"$p_0$")
plt.plot(a3, p1, dashes=[1, 2, 6, 2], color="b", label=r"$p_1$")
plt.plot(a3, p2, color="g", label=r"$p_2$")
plt.xlabel(r"$\alpha_3$")
#plt.ylabel(r"$\omega_{min}$, Гц", fontsize=20)

#plt.title(r"Залежність $\omega_{min}$ від $g_v$")
# + r"($N={}, M={}$)".format(N, M))
#plt.tick_params(axis='both', which='major', labelsize=20)
#plt.tick_params(axis='both', which='minor', labelsize=16)

#plt.axes().set_aspect('equal', 'datalim')
plt.legend(loc='best')

#plt.ylim(ymin=0)
plt.grid()
plt.show()