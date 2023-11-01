import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


def phase_diagram(system, ss_x, ss_y, system_class_list, xmin=-10, xmax=10, ymin=-10, ymax=10, phase_diagram_resolution=20, arrow_scale=100, trajectories=False, saveFig=False, filename='', savefigpath=''):
    # Create a grid of points
    x = np.linspace(xmin, xmax, phase_diagram_resolution)
    y = np.linspace(ymin, ymax, phase_diagram_resolution)

    X, Y = np.meshgrid(x, y)

    # Compute the direction at each grid point
    u, v = np.zeros(X.shape), np.zeros(Y.shape)

    NI, NJ = X.shape
    for i in range(NI):
        for j in range(NJ):
            x = X[i, j]
            y = Y[i, j]
            dx, dy = system([x, y])
            u[i,j] = dx
            v[i,j] = dy

    # Use streamplot to visualize the vector field
    plt.streamplot(X, Y, u, v, color='darkseagreen', density=1.5, linewidth=1, arrowsize=1.5, arrowstyle='->')

    # Plot the steady state as a red dot
    for n in range(len(ss_x)):
        print(system_class_list[n])

        if system_class_list[n] == 'simple stable':
            inside_color = 'seagreen'
            edgecolor = 'seagreen'
        elif system_class_list[n] == 'simple unstable':
            inside_color = 'white'
            edgecolor = 'seagreen'
        else:
            print('a')
            inside_color = 'darkslategrey'
            edgecolor = 'darkslategrey'
        plt.scatter(ss_x[n], ss_y[n], color=inside_color, s=80, zorder=3, edgecolors=edgecolor)

    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    plt.xlabel('$U$')
    plt.ylabel('$V$')

# trajectories
    if trajectories==True:
        random_u0 = np.random.uniform(xmin,xmax,10)
        random_v0 = np.random.uniform(ymin,ymax,10)
        for uv0 in zip(random_u0,random_v0) :
            tspan = np.linspace(0, 50, 1000)

            ys = odeint(system, uv0,tspan)
            plt.plot(ys[:,0], ys[:,1], 'b-', color='darkseagreen',zorder=2) # path
            plt.plot([ys[0,0]], [ys[0,1]], 'o', color='slategrey',zorder=2) # start
            # plt.plot([ys[-1,0]], [ys[-1,1]], 's') # end

    plt.scatter(ss_x, ss_y, color='orangered', s=45)  # s sets the size of the dot
    if saveFig==True:
        plt.savefig('%s%s.pdf'%(savefigpath,filename))
        plt.show()
        plt.close()

    else:
        plt.show()
    
