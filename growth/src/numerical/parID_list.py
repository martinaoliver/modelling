#############
###paths#####
#############
from fileinput import filename
import sys
import os

from importlib_metadata import distribution
pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')
#############



import glob
import pickle


circuit_n='turinghill'
mechanism='nogrowth'
variant=2

L=50; x_gridpoints=5; J=L*x_gridpoints;I=J 
T=5000; t_gridpoints = 30; N=T*t_gridpoints #Number of timepoints

filename= lambda parID: '%s_variant%s_%s_ID%s_L%r_J%r_T%r_N%r'%(circuit_n,variant,mechanism,parID,L,J,T,N)
datafile= modellingpath + '/growth/out/numerical/%s/%s/data/2Dfinal_%s.pkl'%(circuit_n,mechanism,filename('*'))
print(modellingpath + '/growth/out/numerical/%s/%s/data/2Dfinal_%s.pkl'%(circuit_n,mechanism,filename('*')))
files = glob.glob(datafile)
print(len(files))
parID_list = []
for f in files:
    f0 = f.rpartition('ID')[2]
    f1 = f0.rpartition('_L')[0]
    parID_list.append(f1)
print(parID_list[:10])
print(len(parID_list))

pickle.dump( parID_list, open( modellingpath + '/growth/out/numerical/%s/%s/data/parID_list_%s.pkl'%(circuit_n,mechanism,filename('x')), "wb" ) )
print('------')
print(parID_list.count('1'))
