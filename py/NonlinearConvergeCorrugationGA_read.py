import os
import platform
import matplotlib.pyplot as plt

import utils


folder = "./results/corrugated/"

y_per_hr = utils.load_results(folder+"y_per_cf")
x_per_hr = utils.load_results(folder+"x_per_cf")

markers = ['o', '*', 'p', 'd', 'x', 'v', 's', '1']
    
plt.figure()
    
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

#plt.plot(y2D, x, 'ro-', linewidth=2.0, markersize=7, markeredgewidth=2, markeredgecolor='r', markerfacecolor='r', label = "Total nonlinearity")
#plt.plot(y1D2O, x, 'gs--', linewidth=2.0, markersize=7, markeredgewidth=2, markeredgecolor='g', markerfacecolor='g', label = "Second order")
#plt.plot(y1D1O, x, 'bx-.', linewidth=2.0, markersize=7, markeredgewidth=2, markeredgecolor='b', markerfacecolor='b', label = "Mindlin-Reisner")

#plt.plot(ya, x, 'mv:', linewidth=2.0, markersize=7, markeredgewidth=2, markeredgecolor='m', markerfacecolor='m', label = "Analytical")

i = 0
for key, value in y_per_hr.items():
    print('===========')
    print(key)
    print(y_per_hr[key])
    print(x_per_hr[key])
    m=markers[i]
    plt.plot(y_per_hr[key], x_per_hr[key], marker=m, markersize=7, label = r"$g_A = {}$".format(key))
    i += 1

plt.xlabel(r"$\frac{\omega_{NL}}{\omega_{L}}$")
plt.ylabel(r"$\frac{w_{max}}{h}$")

plt.xlim(xmin=0)
plt.legend(loc='best')

plt.grid()
plt.show()
    
    

