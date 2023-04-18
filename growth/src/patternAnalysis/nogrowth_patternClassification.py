#%%
#####
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
from randomfunctions import plot_all_dispersion, plot_highest_dispersion

import pickle
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from math import log10 , floor
from scipy.signal import find_peaks

def NormalizeData(data):
    return (data - np.min(data)) / (np.max(data) - np.min(data)+1e-8)

def countPeaks(U, showPlot1D=True):
    peaks = [0,0]
    peaks[0], _ = find_peaks(U[0], prominence=0.1)
    peaks[1], _ = find_peaks(U[1], prominence=0.1)
    if showPlot1D == True:
        plot1D(U,plotPeaks=True, peaks=peaks)

    return peaks
def parID_surfpattern(parIDss,L,J,T,morphogen=0,record_every_x_hours = 10):
    #data 
    U_record = pickle.load( open(modellingephemeral + '/growth/out/numerical/%s/%s/simulation/%s/2Drecord_%s.pkl'%(circuit_n,mechanism,folder,filename(mechanism,parIDss)), 'rb'))    
    #grids
    dx = float(L)/float(J-1)
    x_grid = np.array([j*dx for j in range(J)])
    reduced_t_grid = np.arange(0,T,record_every_x_hours) 

    #plot
    surfpattern(U_record, [x_grid, reduced_t_grid], 'linear',  morphogen=morphogen, rate=0, savefig=False,filename='',logResults=False,normalize=False)
# %%
def patternClassification(U_final, U_record, normalize=True):
    #check if flat
    relRangeFlat = [(np.amax(U) - np.amin(U))/(np.amax(U)+1e-8) for U in U_final]
    if any(i<0.001 for i in relRangeFlat):
        flat=True
    else:
        flat=False

    #check if converged
    relRangeConverged=[0,0]
    for count,Ux_record in enumerate(U_record):
        relRangeConverged[count] = [(np.amax(x) - np.amin(x))/(np.amax(x)+1e-8) for x in np.transpose(Ux_record[-5:])]
    # if np.amax(relRangeConverged[0])>0.001 or np.amax(relRangeConverged[1])>0.001:
    if np.amax(relRangeConverged[0])>0.001 or np.amax(relRangeConverged[1])>0.001:
        converged=False
    else:
        converged=True

    #check if regular
    U_final_norm = [NormalizeData(U) for U in U_final]
    peaks = countPeaks(U_final_norm, showPlot1D=False)

    #calculate distance between peaks in peak0
    std=[0,0]
    for count,singleUpeak in enumerate(peaks):
        if len(singleUpeak)>2:
            peak_dist = [np.linalg.norm(singleUpeak[i]-singleUpeak[i+1]) for i in range(len(singleUpeak)-1)]
            peak_dist = peak_dist/np.sum(peak_dist)
            std[count] = np.std(peak_dist)
        else:
            std[count] = 1
    if std[0]<0.01 and std[1]<0.01:
        regular=True
    else:
        regular=False
        


    if flat==True and converged==True:
        pattern='Homogeneous'
    if flat==True and converged==False:
        pattern='Temporal Oscillator'
    if flat==False and converged==True and regular==True:
        pattern = 'Stationary regular pattern'
    if flat==False and converged==True and regular==False:
        pattern = 'Stationary irregular pattern'

    if flat==False and converged==False and regular==True:
        pattern = 'Non-Stationary regular pattern'
    if flat==False and converged==False and regular==False:
        pattern = 'Non-Stationary irregular pattern'
    # if var[0]>0.1 and var[1]>0.1:





    return pattern, converged, flat, regular




circuit_n='turinghill'
variant= 4
n_species=2
mechanism='nogrowth'
folder = 'turinghill_variant4_nogrowth'

L=15; dx =0.1; J = int(L/dx)
T =40000; dt = 0.02; N = int(T/dt)
boundaryCoeff=1;rate=0.1


# pattern_df = pickle.load(open( modellingpath + '/growth/out/patternAnalysis/%s/%s/pattern/pattern_df_%s.pkl'%(circuit_n,mechanism,filename(mechanism,'x')), 'rb'))

filename= lambda mechanism, parID: 'circuit%s_variant%s_bc%s_%s_rate%s_ID%s_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundaryCoeff, mechanism,rate,parID,L,J,T,N)
n_param_sets=2000000
lsa_df= pickle.load( open(modellingpath + '/growth/out/analytical/lsa_dataframes/lsa_df_%s_variant%r_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), "rb"))
# turing_df= pickle.load( open(modellingpath + '/growth/out/analytical/turing_data/lsa_df_%s_variant%r_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), "rb"))
parID_list = pickle.load(open( modellingephemeral + '/growth/out/numerical/%s/%s/simulation/%s/parID_list_%s.pkl'%(circuit_n,mechanism,folder,filename(mechanism,'x')), "rb" ) )
# parID_list=['1009985.2', '1009876.3']

start=0;stop=len(parID_list);stop=30
parID_list = [i for i in parID_list[start:stop]] #turn string list into integer list

# parID_list.sort() #sort from lower to higher values
patternDict = {}
# parID_list=[41018,30997,2, 4]
test=True

for count,parIDss in enumerate(tqdm(parID_list, disable=False)):
    print(parID)
    #load records 
    U_final = pickle.load( open(modellingephemeral + '/growth/out/numerical/%s/%s/simulation/%s/2Dfinal_%s.pkl'%(circuit_n,mechanism,folder,filename(mechanism,parIDss)), 'rb'))
    U_final = np.round(U_final,decimals=4)
    U_record = pickle.load( open(modellingephemeral + '/growth/out/numerical/%s/%s/simulation/%s/2Drecord_%s.pkl'%(circuit_n,mechanism,folder,filename(mechanism,parIDss)), 'rb'))    
    
    pattern, converged, flat, regular = patternClassification(U_final, U_record)

    parID,ss = parIDss.split('.')
    patternDict[(parID,ss)]=pattern

    if test==True:
        parID_surfpattern(parIDss,L,J,T,morphogen=1)
        plt.show()
        parID_surfpattern(parIDss,L,J,T,morphogen=0)
        plot1D(U_final)
        print(f'Pattern = {pattern}')
        print(f'Converged = {converged}, Flat = {flat}, Regular = {regular}')
        print('-----------------------------------')
# %%
if test==False:
    n_param_sets = 2000000
    lsa_df= pickle.load( open(modellingpath + '/growth/out/analytical/lsa_dataframes/lsa_df_%s_variant%r_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), "rb"))
    lsa_df['pattern'] = np.nan
    for parID_ss,pattern in tqdm(patternDict.items()):
        lsa_df.loc[(int(parID_ss[0]),int(parID_ss[1])), 'pattern']=pattern

    pickle.dump(lsa_df , open( modellingpath + '/growth/out/patternAnalysis/%s/%s/pattern/pattern_df_%s.pkl'%(circuit_n,mechanism,filename(mechanism,'x')), 'wb'))
    print(lsa_df['pattern'].value_counts())

# %%