import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.patches as patches
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from fem.general2D import matrices2D
from matplotlib import animation
import os
import platform

plot_x1_elements = 400
plot_x2_elements = 20

def init():
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


def plot_strain(result, x1_start, x1_end, x2_start, x2_end, time, strain_index):

    dx1 = (x1_end - x1_start) / plot_x1_elements
    dx2 = (x2_end - x2_start) / plot_x2_elements

    x1 = set()
    x2 = set()
    
    plt.figure()

    v = np.zeros((plot_x1_elements + 1, plot_x2_elements + 1))

    for i in range(plot_x1_elements + 1):
        for j in range(plot_x2_elements + 1):
            x1c = x1_start + i * dx1
            x2c = x2_start + j * dx2
            v[i, j] = result.get_strain(x1c, 0, x2c, time)[strain_index]
            x1.add(x1c)
            x2.add(x2c)

    x1 = sorted(x1)
    x2 = sorted(x2)

    (X1, X2) = np.meshgrid(x1, x2)
    surf = plt.contourf(X1, X2, v.T, cmap=cm.rainbow)
    e_i, e_j = matrices.get_index_conv(strain_index)
    plt.title("Деформації " + r"$E_{{{}{}}}$".format(e_i + 1, e_j + 1))
    plt.colorbar(surf)
    plt.show()
    
    
def plot_strain_2(result, N, M, x1_start, x1_end, x2_start, x2_end, time, strain_index):

    x1 = set()
    x2 = set()
    
    plt.figure()

    v_nl = np.zeros((N + 1, M + 1))
    v = np.zeros((N + 1, M + 1))

    for node in result.mesh.nodes:
        x1c = node.x1
        x2c = node.x2
        j = node.index // (N+1)
        i = node.index % (N+1)
        e_nl = result.get_strain_nl(x1c, 0, x2c, time)
        v_nl[i, j] = e_nl[strain_index]
        e = result.get_strain(x1c, 0, x2c, time)
        v[i, j] = e[strain_index]
        x1.add(x1c)
        x2.add(x2c)

    x1 = sorted(x1)
    x2 = sorted(x2)

    (X1, X2) = np.meshgrid(x1, x2)
    
    plt.subplot(2, 1, 1)
    
    e_i, e_j = matrices.get_index_conv(strain_index)
    plt.title("Деформації " + r"$E_{{{}{}}}$".format(e_i + 1, e_j+1))
    
    surf = plt.contourf(X1, X2, v.T, cmap=cm.rainbow)
    plt.colorbar(surf)
    plt.ylabel("Лінійні")
    
    plt.subplot(2, 1, 2)
    surf = plt.contourf(X1, X2, v_nl.T, cmap=cm.rainbow)
    plt.colorbar(surf)
    plt.ylabel("Нелінійні")
    
    plt.show()


def plot_point_in_time(result, x1, x2, x3, time_points):
    u1 = []
    u2 = []
    u3 = []
    
    for time in time_points:
        u = result.get_displacement_and_deriv(x1, x2, x3, time)
        u1.append(u[0])
        u2.append(u[4])
        u3.append(u[8])


    plt.plot(time_points, u1, "r", label="u1")
    plt.plot(time_points, u3, "b", label="u3")
    plt.title("Переміщення точки")
    # plt.title(r"Форма панелі $L={}, h={}, K={}, g_A={}, g_v={}$".format(x1_end - x1_start, x2_end - x2_start, result.geometry.curvature, result.geometry.corrugation_amplitude, result.geometry.corrugation_frequency))
    plt.axes().set_aspect('equal', 'datalim')
    plt.legend(loc='best')
    plt.xlabel(r"$час$, sec", fontsize=12)
    plt.ylabel(r"$переміщення$, м", fontsize=12)
    plt.grid()
    plt.show()


def plot_init_and_deformed_geometry(result, x1_start, x1_end, x3_start, x3_end, time):

    dx1 = (x1_end - x1_start) / plot_x1_elements

    X_init = []
    Y_init = []
    X_deformed = []
    Y_deformed = []

    x2 = 0
    x3 = x3_end
    for i in range(plot_x1_elements + 1):
        x1 = x1_start + i * dx1
        u = result.get_displacement_and_deriv(x1, x2, x3, time)
        u1 = u[0]
        u2 = u[4]
        u3 = u[8]

        x, y, z = result.geometry.to_cartesian_coordinates(x1, x2, x3)

#        print("==x1 = {}, x2 = {}".format(x1,x2))
#        print("x = {}, y = {}".format(x,y))

        X_init.append(x)
        Y_init.append(z)

        R1x, R1y, R1z = result.geometry.R1(x1, x2, x3)
        R2x, R2y, R2z = result.geometry.R2(x1, x2, x3)
        R3x, R3y, R3z = result.geometry.R3(x1, x2, x3)
        
        x = x + u1*R1x + u2*R2x + u3*R3x
        z = z + u1*R1z + u2*R2z + u3*R3z

        X_deformed.append(x)
        Y_deformed.append(z)

    x3 = x3_start
    for i in range(plot_x1_elements + 1):
        x1 = x1_end - i * dx1
        u = result.get_displacement_and_deriv(x1, x2, x3, time)
        u1 = u[0]
        u2 = u[4]
        u3 = u[8]

        x, y, z = result.geometry.to_cartesian_coordinates(x1, x2, x3)

        X_init.append(x)
        Y_init.append(z)

        R1x, R1y, R1z = result.geometry.R1(x1, x2, x3)
        R2x, R2y, R2z = result.geometry.R2(x1, x2, x3)
        R3x, R3y, R3z = result.geometry.R3(x1, x2, x3)
        
        x = x + u1*R1x + u2*R2x + u3*R3x
        z = z + u1*R1z + u2*R2z + u3*R3z

        X_deformed.append(x)
        Y_deformed.append(z)

    X_init.append(X_init[0])
    Y_init.append(Y_init[0])
    X_deformed.append(X_deformed[0])
    Y_deformed.append(Y_deformed[0])

    plt.plot(X_init, Y_init, "r", label="init configuration")
    plt.plot(X_deformed, Y_deformed, "b", label="current configuration")
    plt.title("Displacement related to minimal natural frequency")
    # plt.title(r"Форма панелі $L={}, h={}, K={}, g_A={}, g_v={}$".format(x1_end - x1_start, x2_end - x2_start, result.geometry.curvature, result.geometry.corrugation_amplitude, result.geometry.corrugation_frequency))
    plt.axes().set_aspect('equal', 'datalim')
    plt.legend(loc='best')
    plt.xlabel(r"$x_1$, m", fontsize=12)
    plt.ylabel(r"$x_3$, m", fontsize=12)
    plt.grid()
    plt.show()
    
def plot_init_and_deformed_geometry_in_cartesian(result, x1_start, x1_end, x3_start, x3_end, time, to_cartesian_coordinates=None):
    
    alphas_toward = np.vstack((np.linspace(x1_start, x1_end, num=plot_x1_elements), np.linspace(x3_end,x3_end,plot_x1_elements)))
    alphas_backward = np.vstack((np.linspace(x1_end, x1_start, num=plot_x1_elements), np.linspace(x3_start,x3_start,plot_x1_elements)))
    
    alphas = np.concatenate((alphas_toward,alphas_backward), axis=1)
    
    rows, cols = alphas.shape
    
    X_init = []
    Y_init = []
    X_deformed = []
    Y_deformed = []

    x2 = 0
    
    for i in range(cols):
        x1 = alphas[0, i]
        x3 = alphas[1, i]
        
        u = result.get_displacement_and_deriv(x1, x2, x3, time)
        u1 = u[0]
        u2 = u[4]
        u3 = u[8]

        x=x1
        y=x2
        z=x3
        if (to_cartesian_coordinates != None):
            x,y,z=to_cartesian_coordinates(x1,x2,x3)
            
        xu=x1 + u1
        yu=x3 + u3
        zu=x3 + u3
        if (to_cartesian_coordinates != None):
            xu,yu,zu=to_cartesian_coordinates(x1 + u1,x2 + u2,x3 + u3)
            
        X_init.append(x)
        Y_init.append(z)
        X_deformed.append(xu)
        Y_deformed.append(zu)

    X_init.append(X_init[0])
    Y_init.append(Y_init[0])
    X_deformed.append(X_deformed[0])
    Y_deformed.append(Y_deformed[0])

    plt.plot(X_init, Y_init, "r", label="init configuration")
    plt.plot(X_deformed, Y_deformed, "b", label="current configuration")
    plt.title("Displacement related to minimal natural frequency")
    # plt.title(r"Форма панелі $L={}, h={}, K={}, g_A={}, g_v={}$".format(x1_end - x1_start, x2_end - x2_start, result.geometry.curvature, result.geometry.corrugation_amplitude, result.geometry.corrugation_frequency))
    plt.axes().set_aspect('equal', 'datalim')
    plt.legend(loc='best')
    plt.xlabel(r"$x_1$, m", fontsize=12)
    plt.ylabel(r"$x_3$, m", fontsize=12)
    plt.grid()
    plt.show()
    

    
def plot_mesh(mesh, L, h):
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111, aspect='equal')
    for fe in mesh.elements:        
        ax2.add_patch(
            patches.Rectangle(
                (fe.bottom_left.x1, fe.bottom_left.x2),
                fe.width(),
                fe.height(),
                fill=False      # remove background
            )
        )
     
    for n in mesh.nodes:        
        plt.text(n.x1, n.x2, "{}".format(n.index))
    
            
    x_eps = L*0.1
    y_eps = h*0.1
            
    plt.xlim([0-x_eps, L+x_eps])
    plt.ylim([-h/2 - x_eps, h/2 + x_eps])
    fig2.show()

def plot_deformed_mesh(result, L, h):
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111, aspect='equal')
    for fe in result.mesh.elements:
        points = []
        u1 = result.u1[fe.top_left.index]
        u3 = result.u3[fe.top_left.index]
        points.append((fe.top_left.x1 + u1, fe.top_left.x2 + u3))
        u1 = result.u1[fe.top_right.index]
        u3 = result.u3[fe.top_right.index]
        points.append((fe.top_right.x1 + u1, fe.top_right.x2 + u3))
        u1 = result.u1[fe.bottom_right.index]
        u3 = result.u3[fe.bottom_right.index]
        points.append((fe.bottom_right.x1 + u1, fe.bottom_right.x2 + u3))
        u1 = result.u1[fe.bottom_left.index]
        u3 = result.u3[fe.bottom_left.index]
        points.append((fe.bottom_left.x1 + u1, fe.bottom_left.x2 + u3))
        ax2.add_patch(
            patches.Polygon(points, closed=True, fill=False)
        )
     
    for n in result.mesh.nodes:     
        u1 = result.u1[n.index]
        u3 = result.u3[n.index]
        plt.text(n.x1 + u1, n.x2 + u3, "{}".format(n.index))
    
            
    x_eps = L*0.1
    y_eps = h*0.1
            
    plt.xlim([0-x_eps, L+x_eps])
    plt.ylim([-h/2 - x_eps, h/2 + x_eps])
    
    plt.grid()
    fig2.show()


def plot_init_geometry(geometry, x1_start, x1_end, x3_start, x3_end, time):

    dx1 = (x1_end - x1_start) / plot_x1_elements

    X_init = []
    Y_init = []

    x2 = 0
    x3 = x3_end
    for i in range(plot_x1_elements + 1):
        x1 = x1_start + i * dx1

        x, y, z = geometry.to_cartesian_coordinates(x1, x2, x3)

        X_init.append(x)
        Y_init.append(z)

    x3 = x3_start
    for i in range(plot_x1_elements + 1):
        x1 = x1_end - i * dx1

        x, y, z = geometry.to_cartesian_coordinates(x1, x2, x3)

        X_init.append(x)
        Y_init.append(z)

    X_init.append(X_init[0])
    Y_init.append(Y_init[0])

    plt.plot(X_init, Y_init, "r", label="init configuration")

#    geometry_title = str(geometry)
#    plot_title = r"Форма панелі $L={}, h={}$".format(x1_end - x1_start, x3_end - x3_start)
#    if (len(geometry_title) > 0):
#        plot_title = r"Форма панелі $L={}, h={}, {}$".format(x1_end - x1_start, x3_end - x3_start, geometry_title)

#    plt.title(plot_title)
    plt.axes().set_aspect('equal', 'datalim')
    plt.legend(loc='best')
    plt.xlabel(r"$x_1$, m", fontsize=12)
    plt.ylabel(r"$x_3$, m", fontsize=12)
    plt.grid()
    plt.show()
    
def plot_init_geometry_2(x1_start, x1_end, x3_start, x3_end, to_cartesian_coordinates):

    dx1 = (x1_end - x1_start) / plot_x1_elements

    X_init = []
    Y_init = []

    x2 = 0
    x3 = x3_end
    for i in range(plot_x1_elements + 1):
        x1 = x1_start + i * dx1

        x, y, z = to_cartesian_coordinates(x1, x2, x3)

        X_init.append(x)
        Y_init.append(z)

    x3 = x3_start
    for i in range(plot_x1_elements + 1):
        x1 = x1_end - i * dx1

        x, y, z =to_cartesian_coordinates(x1, x2, x3)

        X_init.append(x)
        Y_init.append(z)

    X_init.append(X_init[0])
    Y_init.append(Y_init[0])

    plt.plot(X_init, Y_init, "r", label="init configuration")

    plt.axes().set_aspect('equal', 'datalim')
    plt.legend(loc='best')
    plt.xlabel(r"$x_1$, m", fontsize=12)
    plt.ylabel(r"$x_3$, m", fontsize=12)
    plt.grid()
    plt.show()
    plt.savefig('plot1.png')


def plot_freq_from_corrugated_freq(g_v, w_min, N, M):
    plt.plot(g_v, w_min, 'o-', linewidth=3.0, markersize=7, markeredgewidth=1, markeredgecolor='r', markerfacecolor='None')
    plt.xlabel(r"$g_v$", fontsize=20)
    plt.ylabel(r"$\omega_{min}$, Гц", fontsize=20)

    plt.title(r"Залежність $\omega_{min}$ від $g_v$")
    # + r"($N={}, M={}$)".format(N, M))
    plt.tick_params(axis='both', which='major', labelsize=20)
    plt.tick_params(axis='both', which='minor', labelsize=16)

    plt.ylim(ymin=0)
    plt.grid()
    plt.show()
    
def plot_vibrations(time, u_linear, u_nonlinear, u_nonlinear2):
    init()
    plt.plot(time, u_linear, 'go-', linewidth=1.0, markersize=5, markeredgewidth=1, markeredgecolor='g', markerfacecolor='None', label = "linear")
    plt.plot(time, u_nonlinear, 'rv:', linewidth=1.0, markersize=5, markeredgewidth=1, markeredgecolor='r', markerfacecolor='None', label = "nonlinear")
    plt.plot(time, u_nonlinear2, 'bx--', linewidth=1.0, markersize=5, markeredgewidth=1, markeredgecolor='b', markerfacecolor='None', label = "nonlinear2")
    plt.xlabel(r"$time, sec$")
    plt.ylabel(r"$u_3, m$")

    # + r"($N={}, M={}$)".format(N, M))
    plt.legend(loc='best')
#    plt.tick_params(axis='both', which='major', labelsize=20)
#    plt.tick_params(axis='both', which='minor', labelsize=16)

    plt.grid()
    plt.show()
    
    
#os.environ["PATH"] += os.pathsep + '/usr/local/texlive/2017/bin/x86_64-darwin'

def plot_vectors(x1_start, x1_end, x3, to_cartesian_coordinates, normal_in_cartesian_coordinates):
    
    dx1 = (x1_end - x1_start) / plot_x1_elements

    x2 = 0
    
    X_init = []
    Y_init = []
    
    for i in range(plot_x1_elements + 1):
        x1 = x1_start + i * dx1

        x, y, z = to_cartesian_coordinates(x1, x2, x3)

        X_init.append(x)
        Y_init.append(z)
        
    plt.plot(X_init, Y_init, "g", dashes=[4, 4])

    
    arrows_elements = 50
    
    dx1 = (x1_end - x1_start) / arrows_elements
    
    start_arrow_x = []
    start_arrow_y = []
    
    end_arrow_x = []
    end_arrow_y = []

   
    for i in range(arrows_elements + 1):
        x1 = x1_start + i * dx1

        x, y, z = to_cartesian_coordinates(x1, x2, x3)

        start_arrow_x.append(x)
        start_arrow_y.append(z)
        
        n_x, n_y, n_z = normal_in_cartesian_coordinates(x1, x2, x3)

        end_arrow_x.append(n_x)
        end_arrow_y.append(n_z)
        
        
    plt.quiver(start_arrow_x, start_arrow_y, end_arrow_x, end_arrow_y, color='r', width=0.002)


#    geometry_title = str(geometry)
#    start_title = "The shape of the panel with the following parameters:\n"
#    L = x1_end - x1_start
#    h = x2_end - x2_start
#    plot_title = start_title + r"$L={}$m, $h={}$m".format(L, h)
#    if (len(geometry_title) > 0):
#        plot_title += r", {}".format(geometry_title)
        
    
    #fig.text(0.5, 0.01, plot_title, ha='center')

    #plt.title(plot_title)
    plt.gca().set_aspect('equal', 'datalim')
#    plt.legend(loc='best')
    plt.xlabel(r"$x_1$,m")
    plt.ylabel(r"$x_3$,m")
    
    #plt.gca().xaxis.set_label_coords(1, -0.02) 
#    plt.gca().yaxis.set_label_coords(-0.005, 1) 
    plt.grid()
    #fig.subplots_adjust(bottom=0.2) 
    #plt.tight_layout()
    plt.show()
    
    #plt.savefig(filename + '.png')
    
def deformed_shell(result, x1_start, x1_end, x3_start, x3_end, time, to_cartesian_coordinates=None):
    
    alphas_toward = np.vstack((np.linspace(x1_start, x1_end, num=plot_x1_elements), np.linspace(x3_end,x3_end,plot_x1_elements)))
    alphas_backward = np.vstack((np.linspace(x1_end, x1_start, num=plot_x1_elements), np.linspace(x3_start,x3_start,plot_x1_elements)))
    
    alphas = np.concatenate((alphas_toward,alphas_backward), axis=1)
    
    rows, cols = alphas.shape
    
    X_deformed = []
    Y_deformed = []

    x2 = 0
    
    for i in range(cols):
        x1 = alphas[0, i]
        x3 = alphas[1, i]
        
        u = result.get_displacement_and_deriv(x1, x2, x3, time)
        u1 = u[0]
        u2 = u[4]
        u3 = u[8]

        x=x1
        y=x2
        z=x3
        if (to_cartesian_coordinates != None):
            x,y,z=to_cartesian_coordinates(x1,x2,x3)
            
        xu=x1 + u1
        yu=x3 + u3
        zu=x3 + u3
        if (to_cartesian_coordinates != None):
            xu,yu,zu=to_cartesian_coordinates(x1 + u1,x2 + u2,x3 + u3)
            
        X_deformed.append(xu)
        Y_deformed.append(zu)

    X_deformed.append(X_deformed[0])
    Y_deformed.append(Y_deformed[0])
    
    return X_deformed, Y_deformed

def plot_animate_in_cartesian(result, x1_start, x1_end, x3_start, x3_end, T, count, to_cartesian_coordinates=None):
    
    
    deltat = T / count
    
    # First set up the figure, the axis, and the plot element we want to animate
    fig = plt.figure()
    ax = plt.axes()
    x, y = deformed_shell(result, x1_start, x1_end, x3_start, x3_end, 0, to_cartesian_coordinates)
    line, = ax.plot(x, y)
    
    # initialization function: plot the background of each frame
    def init():
        print('init')
        x, y = deformed_shell(result, x1_start, x1_end, x3_start, x3_end, 0, to_cartesian_coordinates)
        line.set_data(x, y)
        return line,
    
    # animation function.  This is called sequentially
    def animate(i):
        time = deltat*i
        print(time)
        x, y = deformed_shell(result, x1_start, x1_end, x3_start, x3_end, time, to_cartesian_coordinates)
        line.set_data(x, y)
        return line,
    
    print('count = {}'.format(count))
#    
#    # call the animator.  blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=count, interval=20, blit=True)
    
#    init()
    
    
    plt.show()
