#############
###paths#####
#############
from fileinput import filename
import sys
import os

pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
modellingephemeral = '/rds/general/user/mo2016/ephemeral/Documents/modelling'
sys.path.append(modellingpath + '/lib')
#############



import glob
import pickle


circuit_n='turinghill'
mechanism='edgegrowth2'
# mechanism='openboundary'
# mechanism='nogrowth'
variant=9
folder = f'turinghill_variant{variant}'
#solver parameters
L=50; dx =0.1; J = int(L/dx)
T =2000; dt = 0.02; N = int(T/dt)

boundaryCoeff=2;rate=L/T
suggesteddt = float(dx*dx*2)



filename= lambda mechanism, parID: 'circuit%s_variant%s_bc%s_%s_rate%s_ID%s_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundaryCoeff, mechanism,rate,parID,L,J,T,N)
print(filename(mechanism, 1))
datafile= modellingephemeral + '/growth/out/numerical/%s/simulation/%s/2Dfinal_%s.pkl'%(mechanism,folder,filename(mechanism, '*'))
print(modellingephemeral + '/growth/out/numerical/%s/simulation/%s/2Dfinal_%s.pkl'%(mechanism,folder,filename(mechanism, '*')))
files = glob.glob(datafile)
print(len(files))
parID_list = []
for f in files:
    f0 = f.rpartition('ID')[2]
    f1 = f0.rpartition('_L')[0]
    parID_list.append(f1)
print(parID_list[:10])
print(len(parID_list))
print( modellingpath + '/growth/out/numerical/%s/simulation/%s/parID_list_%s.pkl'%(mechanism,folder,filename(mechanism, 'x')))

pickle.dump( parID_list, open( modellingephemeral + '/growth/out/numerical/%s/simulation/%s/parID_list_%s.pkl'%(mechanism,folder,filename(mechanism, 'x')), "wb" ) )
print( modellingephemeral + '/growth/out/numerical/%s/simulation/%s/parID_list_%s.pkl'%(mechanism,folder,filename(mechanism, 'x')))
print('------')
print(parID_list.count('47.0'))

