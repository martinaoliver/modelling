import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
cmap = cm.Spectral_r
cmap=cm.coolwarm
cmap = cm.magma
# from sklearn import preprocessing

def plot1D(U,morphogen='both', savefig=False,filename='',savefigpath='',pad=0.001,round=False, plotPeaks=False, peaks=False, L=1):
    if round==True:
        U = np.round(U,decimals=3)
    
    if morphogen == 0:
        plt.plot(U[0], label='U')
    if morphogen ==1: 
        plt.plot(U[1], label='V')
    if morphogen == 'both': 
        fig,ax = plt.subplots()
        ax.plot(U[0], label='U', color='blue')
        ax.set_ylim(np.amin(U[0])-pad, np.amax(U[0])+pad)
        ax.legend(loc=2) #upper left
        # ax.ticklabel_format(useOffset=False)

        ax2=ax.twinx()
        ax2.plot(U[1], label='V', color='red')
        ax2.set_ylim(np.amin(U[1])-pad, np.amax(U[1])+pad)
        ax2.legend(loc=1) #upper right


        # ax.ticklabel_format(useOffset=False)

        locs, labels = plt.xticks()
        new_labels=locs/20
        plt.xticks(ticks=locs, labels=new_labels)
        plt.xlim(0,len(U[0]))
        if plotPeaks==True:
            ax.plot(peaks[0],U[0][peaks[0]], 'o', color='blue')
            ax2.plot(peaks[1],U[1][peaks[1]], 'o', color='red')


    # plt.ticklabel_format(useOffset=False)
    plt.xlabel('Space')
    plt.ylabel('Concentration')
    if savefig==True:
        plt.savefig('%s%s.jpeg'%(savefigpath,filename))
    else:
        plt.show()


def surfpattern(results,L,dx,J,T,record_every_x_hours=10,growth='linear', rate=0, morphogen = 0,savefig=False,filename='1',logResults=False, normalize=False, cmap=cm.magma, space_crop=None):
    
    
    dx = float(L)/float(J-1)
    x_grid = np.array([j*dx for j in range(J)])
    t_grid = np.arange(0,T,10) 
    
    if normalize == True:
        print('NEEDS NORMALIZATION')
    results = results[morphogen]

    values = results.reshape(len(x_grid),len(t_grid))
    x, t = np.meshgrid(x_grid, t_grid)

    # t,x = np.meshgrid(t_grid, x_grid)
    # plt.contourf(t,x,results, cmap=cmap)
    plt.contourf(x,t,results, cmap=cmap)
    if logResults==True:
        plt.colorbar(label='Concentration (logscale)')
    else:
        plt.colorbar()


    plt.ylabel('Time')
    plt.xlabel('Space')
    if savefig==True:
        plt.savefig('%s_overtime.png'%filename)
    # plt.show()
    
def surfpattern1(results,grids,growth='linear', rate=0, morphogen = 0,savefig=False,filename='1',logResults=False, normalize=False):
    fig,ax = plt.subplots(1,len(morphogen))
    for count,n in enumerate(morphogen):
        if normalize == True:
            # results = [preprocessing.normalize(array, norm="l1") for array in results]
            print('NEEDS NORMALIZATION')

        results = results[n]
        x_grid = grids[0] 
        t_grid = grids[1]
        values = results.reshape(len(x_grid),len(t_grid))
        x, t = np.meshgrid(x_grid, t_grid)

        # t,x = np.meshgrid(t_grid, x_grid)
        # plt.contourf(t,x,results, cmap=cmap)
        ax[count].contourf(x,t,results, cmap=cmap)
        if logResults==True:
            plt.colorbar(label='Concentration (logscale)')
        else:
            plt.colorbar()


        plt.ylabel('Time')
        plt.xlabel('Space')
        if savefig==True:
            plt.savefig('%s_overtime.png'%filename)
        # plt.show()

def surfpattern2(results,grids,growth='linear', rate=0, morphogen = 0,savefig=False,filename='1',logResults=False, normalize=False):
    fig,ax = plt.subplots(1,len(morphogen))
    x_grid = grids[0] 
    t_grid = grids[1]
    x, t = np.meshgrid(x_grid, t_grid)


    for count,n in enumerate(morphogen):
        values = results[n]
        ax[count].contourf(x,t,values, cmap=cmap)



    ax[0].set_ylabel('Time')
    ax[0].set_xlabel('Space')
    plt.show()