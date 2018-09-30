#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 30 12:31:04 2018

@author: tarasgoriachko
"""

import matplotlib.pyplot as plt

g_a=[0, 0.015, 0.03, 0.06, 0.1, 0.2, 0.25, 0.3]
w_min=[101, 143, 217, 377, 622, 675, 552, 461]

plt.figure()
plt.plot(g_a, w_min, 'o-', linewidth=3.0, markersize=7, markeredgewidth=2, markeredgecolor='r', markerfacecolor='None')
plt.xlabel(r"$g_A$, м", fontsize=20)
plt.ylabel(r"$\omega_{min}$, Гц", fontsize=20)

#plt.title(r"Залежність $\omega_{min}$ від $g_А$")
    # + r"($N={}, M={}$)".format(N, M))
plt.tick_params(axis='both', which='major', labelsize=20)
plt.tick_params(axis='both', which='minor', labelsize=16)

plt.ylim(ymin=0)
plt.grid()
plt.show()