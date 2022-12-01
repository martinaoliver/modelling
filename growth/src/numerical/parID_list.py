#############
###paths#####
#############
from fileinput import filename
import sys
import os

pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')
#############



import glob
import pickle


circuit_n='turinghill'
mechanism='edgegrowth2'
mechanism='nogrowth'
variant=0

L=500; dx =1; J = int(L/dx)
T =3000; dt = 0.05; N = int(T/dt)
boundaryCoeff=2;rate=0.1


L=50; dx =1; J = int(L/dx)
T =3000; dt = 0.05; N = int(T/dt)
boundaryCoeff=2;rate=0.01

filename= lambda mechanism, parID: 'circuit%s_variant%s_bc%s_%s_rate%s_ID%s_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundaryCoeff, mechanism,rate,parID,L,J,T,N)
print(filename(mechanism, 1))
datafile= modellingpath + '/growth/out/numerical/%s/%s/simulation/2Dfinal_%s.pkl'%(circuit_n,mechanism,filename(mechanism, '*'))
print(modellingpath + '/growth/out/numerical/%s/%s/simulation/2Dfinal_%s.pkl'%(circuit_n,mechanism,filename(mechanism, '*')))
files = glob.glob(datafile)
print(len(files))
parID_list = []
for f in files:
    f0 = f.rpartition('ID')[2]
    f1 = f0.rpartition('_L')[0]
    parID_list.append(f1)
print(parID_list[:10])
print(len(parID_list))
print(modellingpath + '/growth/out/numerical/%s/%s/simulation/parID_list_%s.pkl'%(circuit_n,mechanism,filename(mechanism, 'x')))

pickle.dump( parID_list, open( modellingpath + '/growth/out/numerical/%s/%s/simulation/parID_list_%s.pkl'%(circuit_n,mechanism,filename(mechanism,'x')), "wb" ) )
print('------')
print(parID_list.count('1.0'))
