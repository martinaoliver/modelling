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
from numerical.generalFunctions import round_it


def plotAllFunction(parIDdict, circuit_n, mechanism, filename, round=True,start=0, stop=10, modellingpath=modellingpath, saveFig=True,dpi=2000, tqdm_disable=True,pad=0.01, round_to=3,metric='peakDistVar'):
    fulldataset=len(parIDdict)
    parIDdict = dict(sorted(parIDdict.items(), key=lambda item: item[1])) #sort from lower to higher values
    parIDdict = dict(list(parIDdict.items())[:stop]) #trim to the first stop values

    num=len(parIDdict)
    n_col = int(np.sqrt(num))
    n_row = int(np.floor(num/n_col)+1)    # number of rows in the figure of the cluster


    fig = plt.figure(figsize=(n_col/10+12, n_row/10+12))

    for count,parID in tqdm(enumerate(parIDdict.keys()),disable=tqdm_disable):
        print(parID)

        ax=plt.subplot(n_row,n_col, count+1)
        U = pickle.load( open(modellingpath + '/growth/out/numerical/%s/%s/data/2Dfinal_%s.pkl'%(circuit_n,mechanism,filename(parID)), 'rb'))
        print(np.amax(U[0]), np.amin(U[0]))
        if round==True:
            U = [[round_it(U0x,round_to) for U0x in U0] for U0 in U]
        ax.plot(U[0], label='U', color='green', alpha=0.5)
        ax.set_ylim(np.amin(U[0])-pad, np.amax(U[0])+pad)
        ax.set(yticks=[],xticks=[],yticklabels=[],xticklabels=[])

        ax2=ax.twinx()
        ax2.plot(U[1], label='U', color='red', alpha=0.5)
        ax2.set_ylim(np.amin(U[1])-pad, np.amax(U[1])+pad)
        ax2.set(yticks=[],xticks=[],yticklabels=[],xticklabels=[])
        ax2.set_ylabel('%r-%r'%(parID,np.round(parIDdict[parID],decimals=3)),size= 1,c='y', labelpad=0.35)

        # ax2=ax.twinx()
        # ax2.plot(U[1], label='V', color='red')
        # ax2.set_ylim(np.amin(U[1])-pad, np.amax(U[1])+pad)
        # ax2.legend(loc=1)#upper right


    
    
    if saveFig==False:
        plt.show()
    if saveFig==True:
        if stop==fulldataset:
            plt.savefig(modellingpath + '/growth/out/numerical/%s/%s/largeFig/largeFig_%s_%s.png'%(circuit_n,mechanism,filename('x'),metric),dpi=dpi)
        
        else:
            plt.savefig(modellingpath + '/growth/out/numerical/%s/%s/largeFig/largeFig_%s_%s-%s_%s.png'%(circuit_n,mechanism,filename('x'),start,stop,metric),dpi=dpi)
            print('not full')
            plt.close()
    # plt.savefig(modelling_home + '/3954/numerical_confocal/results/figures/%s/large_images/%s_%s-%s.png'%(shape,filename,start,stop), dpi=2000)
    # plt.savefig(modelling_home + '/3954/numerical_confocal/results/figures/ca/large_images/%s_%s.png'%(filename,details), dpi=2000)
    # plt.savefig(modelling_home + '/3954/numerical_confocal/results/figures/ca/large_images/%s_%s.png'%(filename,details), dpi=2000)
    x='x'
    print(f'Done plotting {filename(x)}')