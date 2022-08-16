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
import seaborn as sns

circuit_n='turinghill'
variant= 0
n_species=2
mechanism='edgegrowth1'
n_param_sets = 100000
L=50; x_gridpoints=5; J=L*x_gridpoints;I=J 
T=2000; t_gridpoints = 30; N=T*t_gridpoints #Number of timepoints
filename= lambda parID: '%s_variant%s_%s_ID%s_L%r_J%r_T%r_N%r'%(circuit_n,variant,mechanism,parID,L,J,T,N)

# lsa_df= pickle.load( open(modellingpath + '/growth/out/patternAnalysis/peakDistVar/peakDistVar_df_%s_variant%r.pkl'%(circuit_n,variant,n_param_sets), "rb"))
lsa_df = pickle.load(open( modellingpath + '/growth/out/patternAnalysis/%s/%s/peakDistVar/peakDistVar_df_%s.pkl'%(circuit_n,mechanism,filename('x')), 'wb'))


# parIDHpsDict = pickle.load(open( modellingpath + '/growth/out/patternAnalysis/%s/%s/parIDHpsDict%s_batch%r.pkl'%(circuit_n,mechanism,filename('x'),30000), 'rb'))


# lsa_df = lsa_df.xs(0, level=1)
# for parID in lsa_df.index:
#     if parID not in parIDHpsDict.keys():
#         parIDHpsDict[parID]= 0 

#add column to lsa_df with Hps and save it to file
lsa_df['Hps'] = lsa_df.index.to_series().map(parIDHpsDict)


boring_states= ['simple stable', 'simple unstable']  
interesting_states = ['turing I, turing II', 'turing I hopf', 'turing I oscillatory', 'turing II hopf','hopf']  
boring_df = lsa_df.loc[lsa_df['system_class'].isin(boring_states)]
interesting_df = lsa_df.loc[lsa_df['system_class'].isin(interesting_states)]



sns.histplot(data=boring_df, x='Hps', log_scale=(False,True), bins=100, hue='system_class', palette='Set1', alpha=0.01)
sns.histplot(data=interesting_df, x='Hps', log_scale=(False,True), bins=100, hue='system_class', palette='Set1', alpha=0.5)
plt.show()

print(lsa_df['system_class'].value_counts())
