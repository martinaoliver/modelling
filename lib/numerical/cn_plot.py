import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm

from matplotlib.colors import ListedColormap
import seaborn as sns
my_cmap = ListedColormap(sns.color_palette("Spectral",256))   
# from sklearn import preprocessing
import seaborn as sns

def plot1D(U,dx=0.05,morphogen='both', savefig=False,filename='',savefigpath='',pad=0.001,round=False, plotPeaks=False, peaks=False, L=1):
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
        new_labels=locs*dx
        plt.xticks(ticks=locs, labels=new_labels)
        plt.xlim(0,len(U[0]))
        if plotPeaks==True:
            ax.plot(peaks[0],U[0][peaks[0]], 'o', color='blue')
            ax2.plot(peaks[1],U[1][peaks[1]], 'o', color='red')


    # plt.ticklabel_format(useOffset=False)
    plt.xlabel('Space')
    plt.ylabel('Concentration')
    if savefig==True:
        plt.savefig('%s%s.pdf'%(savefigpath,filename))
        plt.show()
        plt.close()
    else:
        plt.show()


def surfpattern(results,L,dx,J,T, record_every_x_hours=10,growth='linear', rate=0, morphogen = 0,savefig=False,filename='',savefigpath='',logResults=False, normalize=False, cmap=my_cmap, space_crop=None):
    

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
    plt.contourf(x,t,results, levels=100, cmap=cmap)
    if logResults==True:
        plt.colorbar(label='Concentration (logscale)')
    else:
        plt.colorbar()


    plt.ylabel('Time')
    plt.xlabel('Space')
    if savefig==True:
        plt.savefig('%s%s.pdf'%(savefigpath,filename))
        plt.show()
        plt.close()

    else:
        plt.show()
    
def surfpattern1(results,grids,savefigpath,growth='linear', rate=0, morphogen = 0,savefig=False,filename='1',logResults=False, normalize=False):
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
        contour = ax[count].contourf(x,t,results, levels=100, cmap=cmap)
        fig.colorbar(contour, ax=ax[count], orientation='vertical')

        ax[count].set_ylabel('Time')
        ax[count].set_xlabel('Space')

    if savefig==True:
        plt.savefig('%s%s.pdf'%(savefigpath,filename))
    plt.show()

def surfpattern2(results,grids,savefigpath,growth='linear', rate=0, morphogen = 0,savefig=False,filename='1',logResults=False, normalize=False):
    fig,ax = plt.subplots(1,len(morphogen))
    x_grid = grids[0] 
    t_grid = grids[1]
    x, t = np.meshgrid(x_grid, t_grid)


    for count,n in enumerate(morphogen):
        values = results[n]
        contour = ax[count].contourf(x,t,results, levels=100, cmap=cmap)
        ax[count].set_ylabel('Time')
        ax[count].set_xlabel('Space')

    if savefig==True:
        plt.savefig('%s%s.pdf'%(savefigpath,filename))
    plt.show()




def surfpattern_growth(results,L,dx,J,T, masking=False, record_every_x_hours=10,growth='linear', rate=0, morphogen = 0,savefig=False,filename='',savefigpath='',logResults=False, normalize=False, cmap=my_cmap, space_crop=None):
    
    def create_growth_mask(shape):
        height, width = shape
        middle_point = width // 2
        mask = np.zeros((height, width))

        for t in range(height):
            growth_extent = int(1 + t*width/2/height)
            start = max(middle_point - growth_extent, 0)
            end = min(middle_point + growth_extent, width)
            mask[t, start:end] = 1
            
        return mask

    mask = create_growth_mask(np.shape(results[0]))




    dx = float(L)/float(J-1)
    x_grid = np.array([j*dx for j in range(J)])
    t_grid = np.arange(0,T,10) 
    
    if normalize == True:
        print('NEEDS NORMALIZATION')
    results = results[morphogen]
    if masking == True:
        results = results * mask
        # Create a masked array where zeros are masked
        results = np.ma.masked_where(results == 0, results)


    values = results.reshape(len(x_grid),len(t_grid))
    x, t = np.meshgrid(x_grid, t_grid)

    # t,x = np.meshgrid(t_grid, x_grid)
    # plt.contourf(t,x,results, cmap=cmap)]

    
    print(np.shape(x), np.shape(t), np.shape(results))
    plt.contourf(x,t,results, levels=100, cmap=my_cmap)
    if logResults==True:
        plt.colorbar(label='Concentration (logscale)')
    else:
        plt.colorbar()


    plt.ylabel('Time')
    plt.xlabel('Space')
    if savefig==True:
        plt.savefig('%s%s.pdf'%(savefigpath,filename))
        plt.show()
        plt.close()

    else:
        plt.show()