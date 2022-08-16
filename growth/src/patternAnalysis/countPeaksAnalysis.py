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

from numerical.cn_plot import plot1D, surfpattern
from numerical.countPeaksAnalysisFunctions import countPeaks, varPeakDistFunction

from scipy.signal import find_peaks
import pickle
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

circuit_n='turinghill'
variant= 0
n_species=2
mechanism='edgegrowth1'
L=50; x_gridpoints=5; J=L*x_gridpoints;I=J 
T=2000; t_gridpoints = 30; N=T*t_gridpoints #Number of timepoints
filename= lambda parID: '%s_variant%s_%s_ID%s_L%r_J%r_T%r_N%r'%(circuit_n,variant,mechanism,parID,L,J,T,N)
n_param_sets=100000


parID_list = pickle.load(open( modellingpath + '/growth/out/numerical/%s/%s/data/parID_list_%s.pkl'%(circuit_n,mechanism,filename('x')), "rb" ) )
start=0
stop=len(parID_list)
parID_list = [int(i) for i in parID_list[start:stop]] #turn string list into integer list
parID_list.sort() #sort from lower to higher values
parIDvarDict = {}
# parID_list=[41018, 30997, 18622,28817]
# parID_list=[76,167,1]

test=False
for count,parID in enumerate(tqdm(parID_list, disable=True)):
    U = pickle.load( open(modellingpath + '/growth/out/numerical/%s/%s/data/2Dfinal_%s.pkl'%(circuit_n,mechanism,filename(parID)), 'rb'))

    varPeakDist = varPeakDistFunction(U)
    meanvarPeakDist= np.mean(varPeakDist)
    parIDvarDict[parID] = meanvarPeakDist
    if test==True:
        plot1D(U)
        print(varPeakDist)
    if count%10000==0:
        pickle.dump(parIDvarDict, open( modellingpath + '/growth/out/patternAnalysis/%s/%s/peakDistVar/peakDistVar_%s_batch%r.pkl'%(circuit_n,mechanism,filename('x'), count), 'wb'))

if test==False:
    pickle.dump(parIDvarDict, open( modellingpath + '/growth/out/patternAnalysis/%s/%s/peakDistVar/peakDistVar_%s.pkl'%(circuit_n,mechanism,filename('x')), 'wb'))
    print('saved')

lsa_df= pickle.load( open(modellingpath + '/growth/out/analytical/lsa_dataframes/lsa_df_%s_variant%r_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), "rb"))
lsa_df = lsa_df.xs(0, level=1)
for parID in lsa_df.index:
    if parID not in parIDvarDict.keys():
        parIDvarDict[parID]= 0 

if test==False:
    lsa_df['peakDistVar'] = lsa_df.index.to_series().map(parIDvarDict)
    pickle.dump(lsa_df , open( modellingpath + '/growth/out/patternAnalysis/%s/%s/peakDistVar/peakDistVar_df_%s.pkl'%(circuit_n,mechanism,filename('x')), 'wb'))
    # #add column to lsa_df with Hps and save it to file
    # lsa_df= pickle.load( open(modellingpath + '/growth/out/analytical/lsa_dataframes/lsa_df_%s_variant%r_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), "rb"))
    # lsa_df_single = lsa_df.xs(0, level=1)
    print(lsa_df.head())