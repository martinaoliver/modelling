#############
###paths#####
#############
import sys
import os

from importlib_metadata import distribution
pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')
#############

import numpy as np
import matplotlib.pyplot as plt
import pickle 
from tqdm import tqdm


def plotAllFunction(parIDdict, circuit_n, mechanism, filename, start=0, stop=10, modellingpath=modellingpath, saveFig=True,dpi=2000, tqdm_disble=True):
    fulldataset=len(parIDdict)
    parIDdict = dict(sorted(parIDdict.items(), key=lambda item: item[1])) #sort from lower to higher values
    parIDdict = dict(list(parIDdict.items())[:stop]) #trim to the first stop values

    num=len(parIDdict)
    n_col = int(np.sqrt(num))
    n_row = int(np.floor(num/n_col)+1)    # number of rows in the figure of the cluster


    fig = plt.figure(figsize=(n_col/10+12, n_row/10+12))

    for count,parID in tqdm(enumerate(parIDdict.keys()),disable=tqdm_disble):
        ax=plt.subplot(n_row,n_col, count+1)
        U = pickle.load( open(modellingpath + '/growth/out/numerical/%s/%s/data/2Dfinal_%s.pkl'%(circuit_n,mechanism,filename(parID)), 'rb'))
        U=np.round(U,decimals=3)
        pad=0.001
        ax.plot(U[0], label='U', color='green', alpha=0.5)
        ax.set_ylim(np.amin(U[0])-pad, np.amax(U[0])+pad)

        ax.plot(U[1], label='U', color='red', alpha=0.5)
        ax.set_ylim(np.amin(U[1])-pad, np.amax(U[1])+pad)
        ax.set(yticks=[],xticks=[],yticklabels=[],xticklabels=[])
        ax.set_ylabel('%r-%r'%(parID,np.round(parIDdict[parID],decimals=3)),size= 1,c='y', labelpad=0.35)
    if saveFig==False:
        plt.show()
    if saveFig==True:
        if stop==fulldataset:
            plt.savefig(modellingpath + '/growth/out/numerical/%s/%s/largeFig/largeFig_%s.png'%(circuit_n,mechanism,filename('x')),dpi=dpi)
        
        else:
            plt.savefig(modellingpath + '/growth/out/numerical/%s/%s/largeFig/largeFig_%s_%s-%s.png'%(circuit_n,mechanism,filename('x'),start,stop),dpi=dpi)
            print('not full')
            plt.close()
    # plt.savefig(modelling_home + '/3954/numerical_confocal/results/figures/%s/large_images/%s_%s-%s.png'%(shape,filename,start,stop), dpi=2000)
    # plt.savefig(modelling_home + '/3954/numerical_confocal/results/figures/ca/large_images/%s_%s.png'%(filename,details), dpi=2000)
    # plt.savefig(modelling_home + '/3954/numerical_confocal/results/figures/ca/large_images/%s_%s.png'%(filename,details), dpi=2000)
    x='x'
    print(f'Done plotting {filename(x)}')