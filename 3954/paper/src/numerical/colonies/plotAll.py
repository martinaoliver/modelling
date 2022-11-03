#############
###paths#####
#############
import sys
import os

pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')
#############

from numerical.plotAllFunctions import plotAllFunctionColonies

import pickle
import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib import colors
# from mpl_toolkits.axes_grid1 import make_axes_locatable


#############################

circuit_n=2
variant= 0#'48257gaussian0.1nsr'
# variant= '48257gaussian0.21nsr'
n_species=6
shape='ca'
folder = 'circuit2variant0_1M'
# folder = 'circuit2variant0_instabilities'
# folder='circuit2variant48257gaussian0.21nsr'
boundarycoeff = 1.7
p_division=0.5;seed=1

L=8; dx =0.05; J = int(L/dx)
T =125; dt = 0.05; N = int(T/dt)
# T =1; dt = 0.05; N = int(T/dt)
x_gridpoints=int(1/dx)

filename= lambda parID: 'circuit%r_variant%s_bc%s_%s_ID%s_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundarycoeff, shape,parID,L,J,T,N)
data_path = modellingpath + '/3954/paper/out/numerical/colonies/simulation/%s/'%(folder)
parID_list = pickle.load( open(data_path + '/parID_list_%s.pkl'%(filename('x')), "rb" ) )

start = 0
stop = len(parID_list) 
# stop=10
plotAllFunctionColonies(parID_list, circuit_n, shape, filename, L,x_gridpoints,folder=folder,start=start, stop=stop, tqdm_disable=False, saveFig=True)


# df = pickle.load( open( modellingpath + '/growth/out/patternAnalysis/%s/%s/%s/%s_df_%s.pkl'%(circuit_n,mechanism,metric,metric,filename('x')), 'rb'))
# print(df['system_class'].value_counts())

    # pickle.dump(lsa_df , open( modellingpath + '/growth/out/patternAnalysis/%s/%s/%s/%s_df_%s.pkl'%(circuit_n,mechanism,metric,metric,filename('x')), 'wb'))
