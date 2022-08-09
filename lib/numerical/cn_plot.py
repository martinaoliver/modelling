import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
cmap = cm.Spectral_r
from sklearn import preprocessing

def plot1D(U,morphogen='both', savefig=False,filename='1',pad=0.01):
    U = np.round(U,decimals=3)
    print('agr')
    if morphogen == 0:
        plt.plot(U[0], label='U')
    if morphogen ==1: 
        plt.plot(U[1], label='V')
    if morphogen == 'both': 
        fig,ax = plt.subplots()
        ax.plot(U[0], label='U', color='blue')
        ax.set_ylim(np.amin(U[0])+pad, np.amax(U[0])-pad)
        ax.legend(loc=2) #upper left
        ax.ticklabel_format(useOffset=False)

        ax2=ax.twinx()
        ax2.plot(U[1], label='V', color='red')
        print(np.amax(U[1])) 
        ax2.set_ylim(np.amin(U[1])+pad, np.amax(U[1])-0.9*pad)
        ax2.legend(loc=1)#upper right

        ax.ticklabel_format(useOffset=False)




    plt.ticklabel_format(useOffset=False)
    plt.xlabel('Space')
    plt.ylabel('Time')
    if savefig==True:
        plt.savefig('%s_final.png'%filename)

    plt.show()
    # return fig


def surfpattern(results,grids,growth='linear', rate=0, morphogen = 0,savefig=False,filename='1',logResults=False, normalize=False):
    if normalize == True:
        results = [preprocessing.normalize(array, norm="l1") for array in results]

    results = np.transpose(results[morphogen])
    x_grid = grids[0]
    t_grid = grids[1]
    values = results.reshape(len(t_grid),len(x_grid))
    x, t = np.meshgrid(x_grid, t_grid)
    plt.contourf(x,t,results, cmap=cmap)
    if logResults==True:
        plt.colorbar(label='Concentration (logscale)')
    else:
        plt.colorbar()


    plt.ylabel('Time')
    plt.xlabel('Space')
    if savefig==True:
        plt.savefig('%s_overtime.png'%filename)
    plt.show()
