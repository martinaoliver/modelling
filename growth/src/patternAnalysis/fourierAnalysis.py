#############
###paths#####
#############
import sys
import os

pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')
#############

from numerical.cn_plot import plot1D, surfpattern
from numerical.fourierAnalysisFunctions import psEntropyFunction, plotFourier
from numerical.generalFunctions import round_it
import pickle
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from math import log10 , floor

circuit_n='turinghill'
variant= 0
n_species=2
mechanism='nogrowth'
L=50; x_gridpoints=5; J=L*x_gridpoints;I=J 
T=2000; t_gridpoints = 25; N=T*t_gridpoints #Number of timepoints
filename= lambda parID: '%s_variant%s_%s_ID%s_L%r_J%r_T%r_N%r'%(circuit_n,variant,mechanism,parID,L,J,T,N)

parID_list = pickle.load(open( modellingpath + '/growth/out/numerical/%s/%s/data/parID_list_%s.pkl'%(circuit_n,mechanism,filename('x')), "rb" ) )
start=0
stop=len(parID_list)
parID_list = [int(i) for i in parID_list[start:stop]] #turn string list into integer list
parID_list.sort() #sort from lower to higher values
parIDPsEntropy = {}
# parID_list=[41018,30997,2, 4]
test=False
for count,parID in enumerate(tqdm(parID_list, disable=False)):
    U = pickle.load( open(modellingpath + '/growth/out/numerical/%s/%s/data/2Dfinal_%s.pkl'%(circuit_n,mechanism,filename(parID)), 'rb'))
    if test==True:
        plt.subplot(121)
        plot1D(U)

    round_to = 3
    U = [[round_it(Uxx,round_to) for Uxx in Ux] for Ux in U]
    H = [psEntropyFunction(Ux) for Ux in U]
    meanH = np.mean(H )
    parIDPsEntropy[parID] = meanH
    if test==True:
        print(H)
        plt.subplot(122)
        plot1D(U)
        plt.show()
        plt.subplot(121)
        plotFourier(U[0])
        plt.subplot(122)
        plotFourier(U[1])
        plt.show()
        print('--------------')
        print('--------------')
        print('--------------')



    if count%10000==0:
        pickle.dump(parIDPsEntropy, open( modellingpath + '/growth/out/patternAnalysis/%s/%s/psEntropy/parIDpsEntropyDict_%s_batch%r.pkl'%(circuit_n,mechanism,filename('x'), count), 'wb'))

if test==False:
    pickle.dump(parIDPsEntropy, open( modellingpath + '/growth/out/patternAnalysis/%s/%s/psEntropy/parIDpsEntropyDict_%s.pkl'%(circuit_n,mechanism,filename('x')), 'wb'))
    print('saved')

#add column to lsa_df with Hps and save it to file
if test==False:
    n_param_sets = 100000
    lsa_df= pickle.load( open(modellingpath + '/growth/out/analytical/lsa_dataframes/lsa_df_%s_variant%r_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), "rb"))
    lsa_df = lsa_df.xs(0, level=1)
    for parID in lsa_df.index:
        if parID not in parIDPsEntropy.keys():
            parIDPsEntropy[parID]= 0 

    lsa_df['psEntropy'] = lsa_df.index.to_series().map(parIDPsEntropy)
    pickle.dump(lsa_df , open( modellingpath + '/growth/out/patternAnalysis/%s/%s/psEntropy/psEntropy_df_%s.pkl'%(circuit_n,mechanism,filename('x')), 'wb'))
    # #add column to lsa_df with Hps and save it to file
    # lsa_df= pickle.load( open(modellingpath + '/growth/out/analytical/lsa_dataframes/lsa_df_%s_variant%r_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), "rb"))
    # lsa_df_single = lsa_df.xs(0, level=1)
    print(lsa_df.head())