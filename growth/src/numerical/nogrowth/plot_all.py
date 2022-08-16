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
T=2000; t_gridpoints = 25; N=T*t_gridpoints #Number of timepoints
filename= lambda parID: '%s_variant%s_%s_ID%s_L%r_J%r_T%r_N%r'%(circuit_n,variant,mechanism,parID,L,J,T,N)
  
metric = 'peakDistVar'
peakDistVarDict = pickle.load( open( modellingpath + '/growth/out/patternAnalysis/%s/%s/%s/%s_%s.pkl'%(circuit_n,mechanism,metric,metric,filename('x')), 'rb'))
underoneDict = {parID:peakDistVar for parID,peakDistVar in peakDistVarDict.items() if peakDistVar<1}

start = 0
stop = len(underoneDict) 
stop=1000
plotAllFunction(underoneDict, circuit_n, mechanism, filename, start=start, stop=stop, tqdm_disable=False, saveFig=True,pad=0)

