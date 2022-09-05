#############
###paths#####
#############
import sys
import os

pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')
#############

from numerical.plotAllFunctions import plotAllFunction

import pickle
import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib import colors
# from mpl_toolkits.axes_grid1 import make_axes_locatable


#############################

circuit_n='turinghill'
variant= 0
n_species=2
mechanism='nogrowth'
L=50; x_gridpoints=5; J=L*x_gridpoints;I=J 
T=5000; t_gridpoints = 30; N=T*t_gridpoints #Number of timepoints
filename= lambda parID: '%s_variant%s_%s_ID%s_L%r_J%r_T%r_N%r'%(circuit_n,variant,mechanism,parID,L,J,T,N)
  
metric = 'psEntropy'
metricDict = pickle.load( open( modellingpath + '/growth/out/patternAnalysis/%s/%s/%s/parID%sDict_%s.pkl'%(circuit_n,mechanism,metric,metric,filename('x')), 'rb'))
# metricDict = pickle.load( open( modellingpath + '/growth/out/patternAnalysis/%s/%s/%s/parID%sDict_%s_batch60000.pkl'%(circuit_n,mechanism,metric,metric,filename('x')), 'rb'))
if metric=='psEntropy':
    patternsDict = {parID:psEntropy for parID,psEntropy in metricDict.items() if psEntropy<5.5174528964647065}
elif metric=='peakDistVar':
    patternsDict = {parID:peakDistVar for parID,peakDistVar in metricDict.items() if peakDistVar!=0}
sortedDict = dict(sorted(metricDict.items(), key=lambda item: item[1])) #important to sort out dictionary by metric values
# del_list = [ 10718,11190, 10393, 10177,100,1000 ,101,10486,10749,10873, 10,10667, 10558,10811] 
# [sortedDict.pop(parID, None) for parID in del_list]
print(patternsDict)
# print(sortedDict)
start = 0
stop = len(patternsDict) 
# stop=10
plotAllFunction(patternsDict, circuit_n, mechanism, filename, start=start, stop=stop,round=False, tqdm_disable=False, saveFig=False,pad=0, metric=metric)


# df = pickle.load( open( modellingpath + '/growth/out/patternAnalysis/%s/%s/%s/%s_df_%s.pkl'%(circuit_n,mechanism,metric,metric,filename('x')), 'rb'))
# print(df['system_class'].value_counts())

    # pickle.dump(lsa_df , open( modellingpath + '/growth/out/patternAnalysis/%s/%s/%s/%s_df_%s.pkl'%(circuit_n,mechanism,metric,metric,filename('x')), 'wb'))
