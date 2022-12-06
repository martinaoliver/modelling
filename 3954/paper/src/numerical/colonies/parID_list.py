#############
###paths#####
#############
import sys
import os

pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')
#############



import glob
import pickle


circuit_n=2
shape='ca'
variant=1
# variant='48257gaussian0.21nsr'

boundarycoeff = 1.7
p_division=0.5;seed=1
folder = 'circuit2variant1_turing'
# folder = 'circuit2variant48257gaussian0.21nsr'
L=8; dx =0.05; J = int(L/dx)
T =125; dt = 0.05; N = int(T/dt)


L=4; dx =0.05; J = int(L/dx)
T =65; dt = 0.05; N = int(T/dt)
boundarycoeff = 1.7
p_division=0.5;seed=1
divisionTimeHours = 1
# T =1; dt = 0.05; N = int(T/dt)
#------
L=4; dx =0.025; J = int(L/dx)
T =65; dt = 0.0025; N = int(T/dt)
boundarycoeff = 1.7
divisionTimeHours=0.5
p_division=0.5;seed=1

filename= lambda parID: 'circuit%r_variant%s_bc%s_%s_ID%s_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundarycoeff, shape,parID,L,J,T,N)
print(f'filename: {filename(1)}')
datafile= modellingpath + '/3954/paper/out/numerical/colonies/simulation/%s/2Dfinal_%s.pkl'%(folder,filename('*'))
print(datafile)
files = glob.glob(datafile)
print(len(files))
parID_list = []
for f in files:
    f0 = f.rpartition('ID')[2]
    f1 = f0.rpartition('_L')[0]
    parID_list.append(f1)
print(parID_list[:10])
print(len(parID_list))

pickle.dump( parID_list, open( modellingpath + '/3954/paper/out/numerical/colonies/simulation/%s/parID_list_%s.pkl'%(folder,filename('x')), "wb" ) )
print('------')
print(parID_list.count('1'))
