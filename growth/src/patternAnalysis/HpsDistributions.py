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
import pickle
import matplotlib.pyplot as plt
circuit_n='turinghill'
variant= 0
n_species=2
mechanism='nogrowth'
L=50; x_gridpoints=5; J=L*x_gridpoints;I=J 
T=2000; t_gridpoints = 25; N=T*t_gridpoints #Number of timepoints
filename= lambda parID: '%s_variant%s_%s_ID%s_L%r_J%r_T%r_N%r'%(circuit_n,variant,mechanism,parID,L,J,T,N)

parIDHpsDict = pickle.load(open( modellingpath + '/growth/out/patternAnalysis/%s/%s/parIDHpsDict%s_batch%r.pkl'%(circuit_n,mechanism,filename('x'),30000), 'rb'))
plt.hist(parIDHpsDict.values(),log=True, bins=100)
plt.show()