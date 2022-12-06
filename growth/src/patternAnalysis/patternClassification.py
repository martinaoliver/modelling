#############
###paths#####
#############
import sys
import os

pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
modellingephemeral = '/rds/general/ephemeral/user/mo2016/ephemeral/Documents/modelling'
sys.path.append(modellingpath + '/lib')
#############

from numerical.cn_plot import plot1D, surfpattern
# from numerical.fourierAnalysisFunctions import psEntropyFunction, plotFourier
from numerical.generalFunctions import round_it
from analytical.linear_stability_analysis import detailed_turing_analysis_dict
from randomfunctions import plot_all_dispersion

import pickle
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from math import log10 , floor
from scipy.signal import find_peaks

circuit_n='turinghill'
variant= 0
n_species=2
mechanism='nogrowth'
L=50; dx =1; J = int(L/dx)
T =3000; dt = 0.05; N = int(T/dt)
boundaryCoeff=2;rate=0.01
filename= lambda mechanism, parID: 'circuit%s_variant%s_bc%s_%s_rate%s_ID%s_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundaryCoeff, mechanism,rate,parID,L,J,T,N)
n_param_sets=2000000
lsa_df= pickle.load( open(modellingpath + '/growth/out/analytical/lsa_dataframes/lsa_df_%s_variant%r_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), "rb"))
# single_lsa_df =  lsa_df.xs(0, level=1)
print('done')

def parID_display(parID, L,J,T,round=False,crop=100,ss_n=0,single_lsa_df = lsa_df):
    #data
    U = pickle.load( open(modellingpath + '/growth/out/numerical/%s/%s/simulation/2Dfinal_%s.pkl'%(circuit_n,mechanism,filename(mechanism,parID)), 'rb'))
    #plot
    plot1D(U,round=round)

    plt.subplots(figsize=(10,4))

    #dispersion
    plt.subplot(121)
    parID_dispersion(parID,crop,ss_n)
    #convergence
    plt.subplot(122)

    parID_surfpattern(parID,L,J,T)
    plt.tight_layout()
    plt.show()
def parID_surfpattern(parID,L,J,T,record_every_x_hours = 10):
    #simulation 
    U_record = pickle.load( open(modellingephemeral + '/growth/out/numerical/%s/%s/simulation/2Drecord_%s.pkl'%(circuit_n,mechanism,filename(mechanism,parID)), 'rb'))
    
    #grids
    dx = float(L)/float(J-1)
    x_grid = np.array([j*dx for j in range(J)])
    reduced_t_grid = np.arange(0,T,record_every_x_hours) 

    #plot
    surfpattern(U_record, [x_grid, reduced_t_grid], 'linear',  morphogen=1, rate=0, savefig=False,filename='',logResults=False,normalize=False)
def parID_dispersion(parID,crop, ss_n):
    #dispersion
    par_dict = single_lsa_df.loc[parID].to_dict() #converts a dataframe row into a dictionary outputing a dictionary for a specific parameter set
    out = detailed_turing_analysis_dict(par_dict, circuit_n, n_species)
    plot_all_dispersion(out[-3][ss_n],2, crop=crop)



parID_list = pickle.load(open( modellingpath + '/growth/out/numerical/%s/%s/simulation/parID_list_%s.pkl'%(circuit_n,mechanism,filename(mechanism,'x')), "rb" ) )
start=0
stop=len(parID_list)
# stop=10
parID_list = [i for i in parID_list[start:stop]] #turn string list into integer list
# parID_list.sort() #sort from lower to higher values
patternDict = {}
# parID_list=[41018,30997,2, 4]
test=False
for count,parIDss in enumerate(tqdm(parID_list, disable=False)):
    # print(parID)
    #load records 
    U_final = pickle.load( open(modellingpath + '/growth/out/numerical/%s/%s/simulation/2Dfinal_%s.pkl'%(circuit_n,mechanism,filename(mechanism,parIDss)), 'rb'))
    U_record = pickle.load( open(modellingephemeral + '/growth/out/numerical/%s/%s/simulation/2Drecord_%s.pkl'%(circuit_n,mechanism,filename(mechanism,parIDss)), 'rb'))
    

    # parID_display(parID,L,J,T, crop=100)

    #check if flat with diff function
    diffUfinal = np.round(np.diff(U_final),decimals=3)
    if np.all(diffUfinal==0)==True:
        flat=True
    else:
        flat=False

    #check if convergence with diff function
    diffUrecord = np.round(np.diff(U_record[1][-10:], axis=0),decimals=6)
    
    if np.all(diffUrecord==0) ==True:
        converged=True
    else:
        converged=False


    if flat==True and converged==True:
        pattern='Homogeneous'
    if flat==True and converged==False:
        pattern='Temporal Oscillator'
    if flat==False and converged==True:
        # pattern = 'Stationary periodic wave'
        pattern = 'Stationary spatial wave'
    if flat==False and converged==False:
        pattern = 'Non-Stationary heterogeneity'

    # patternDict[parIDss] = pattern
    parID,ss = parIDss.split('.')
    patternDict[(parID,ss)]=pattern
#     if count%10000==0:
#         pickle.dump(parIDPsEntropy, open( modellingpath + '/growth/out/patternAnalysis/%s/%s/psEntropy/parIDpsEntropyDict_%s_batch%r.pkl'%(circuit_n,mechanism,filename('x'), count), 'wb'))

# if test==False:
#     pickle.dump(parIDPsEntropy, open( modellingpath + '/growth/out/patternAnalysis/%s/%s/psEntropy/parIDpsEntropyDict_%s.pkl'%(circuit_n,mechanism,filename('x')), 'wb'))
#     print('saved')
# test=False
# parIDPsEntropy = pickle.load(open( modellingpath + '/growth/out/patternAnalysis/%s/%s/psEntropy/parIDpsEntropyDict_%s.pkl'%(circuit_n,mechanism,filename('x')), 'rb'))
# #add column to lsa_df with Hps and save it to file
if test==False:
    n_param_sets = 2000000
    lsa_df= pickle.load( open(modellingpath + '/growth/out/analytical/lsa_dataframes/lsa_df_%s_variant%r_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), "rb"))
    lsa_df['pattern'] = np.nan
    for parID_ss,pattern in tqdm(patternDict.items()):
        lsa_df.loc[(int(parID_ss[0]),int(parID_ss[1])), 'pattern']=pattern

    # lsa_df['pattern'] = lsa_df.index #you get the index as tuple
    # # lsa_df['pattern'] = np.nan
    # lsa_df['pattern'] = lsa_df['pattern'].map(patternDict)
    # lsa_df = lsa_df.set_index('pattern', append=True)


    # lsa_df['pattern'] = lsa_df.index.to_series().map(patternDict)
    # print(lsa_df)
    pickle.dump(lsa_df , open( modellingpath + '/growth/out/patternAnalysis/%s/%s/pattern/pattern_df_%s.pkl'%(circuit_n,mechanism,filename(mechanism,'x')), 'wb'))
#     # #add column to lsa_df with Hps and save it to file
#     # lsa_df= pickle.load( open(modellingpath + '/growth/out/analytical/lsa_dataframes/lsa_df_%s_variant%r_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), "rb"))
#     # lsa_df_single = lsa_df.xs(0, level=1)
#     print(lsa_df.head())
    # print(lsa_df)
    print(lsa_df['pattern'].value_counts())
