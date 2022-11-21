#############
###paths#####
#############
import os
import sys

from importlib_metadata import distribution

pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')

import pickle
import time

import matplotlib.pyplot as plt
import numpy as np
#############
###Imports#####
#############
from numerical.adi_ca_function_openclosed_nodilution import \
    adi_ca_openclosed_nodilution
from numerical.adi_ca_function_openclosed_nodilution_preMask import \
    adi_ca_openclosed_nodilution_preMask
from numerical.adi_ca_function_openclosed_nodilution_preMask_numba import \
    adi_ca_openclosed_nodilution_preMask as \
    adi_ca_openclosed_nodilution_preMask_numba
from numerical.adi_square_function import adi
from numerical.plotting_numerical import *

# %matplotlib inline
#############
###execution parameters#####
#############
# %matplotlib inline
shape = 'ca'
# parID = int(sys.argv[1])
parID = 524126
parID = 546780
parID = 964698
# parID = 824507
circuit_n=2
variant= 1
n_species=6

folder = 'circuit2variant1_turing'
nsamples =  1000000
save_figure = False
tqdm_disable = False #disable tqdm
# boundarycoeff = float(sys.argv[6])


# open parameter dictionaries
df= pickle.load( open(modellingpath + '/3954/paper/input/lhs_parameterfiles/df_circuit%s_variant%s_%rparametersets.pkl'%(circuit_n,variant,nsamples), "rb" ) )
# instabilities_df= pickle.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/instabilities_dataframes/instability_df_circuit%s_variant%r_%rparametersets.pkl'%(circuit_n,variant,nsamples), "rb" ) )
# df= pickle.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/turing_dataframes/turing_df_circuit%s_variant%r_%rparametersets.pkl'%(circuit_n,variant,nsamples), "rb" ) )


#solver parameters
L=8; dx =0.05; J = int(L/dx)
T =125; dt = 0.05; N = int(T/dt)
divisionTimeHours=1
# L=8; dx =0.05; J = int(L/dx)
# T =125; dt = 0.05; N = int(T/dt)
boundarycoeff = 1.7
p_division=0.5;seed=1

cell_matrix_record = pickle.load( open(modellingpath + "/3954/paper/out/numerical/masks/caMask_seed%s_pdivision%s_L%s_J%s_T%s_N%s.pkl"%(seed,p_division,L,J,T,N), "rb" ) )
daughterToMotherDictList = pickle.load( open(modellingpath + "/3954/paper/out/numerical/masks/caMemory_seed%s_pdivision%s_L%s_J%s_T%s_N%s.pkl"%(seed,p_division,L,J,T,N), "rb" ) )
# T =1; dt = 0.05; N = int(T/dt)
# 0.5,0.02, 0.005
filename= lambda parID: 'circuit%r_variant%s_bc%s_%s_ID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundarycoeff, shape,parID,L,J,T,N)


# for parID in range(10):
par_dict = df.loc[parID].to_dict()

D = np.zeros(n_species)
D[:2] = [par_dict['DA'],par_dict['DB'] ]

print(par_dict)

# U_record,U_final =  adi_ca_openclosed_nodilution_preMask(par_dict,L,dx,J,T,dt,N, circuit_n, n_species,D,cell_matrix_record, daughterToMotherDictList,tqdm_disable=False, p_division=0.5,stochasticity=0, seed=1,growth='Slow', boundarycoeff=boundarycoeff)
# U_record,U_final =  adi(par_dict,L,L,J,J,T,N, circuit_n, n_species,D,tqdm_disable=False,stochasticity=0, steadystates=0)
# get the start time
st = time.time()
U_record,U_final =  adi_ca_openclosed_nodilution_preMask_numba(par_dict,L,dx,J,T,dt,N, circuit_n, n_species,D,cell_matrix_record, daughterToMotherDictList,tqdm_disable=False,divisionTimeHours=divisionTimeHours, stochasticity=0, seed=1, boundarycoeff=boundarycoeff)
elapsed_time = time.time() - st
print('Execution time numba:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
plt.imshow(U_final[-1])
plt.show()
print(U_final[0,12,12])



pickle.dump(U_final, open(modellingpath + '/3954/paper/out/numerical/colonies/simulation/%s/2Dfinal_%s.pkl'%(folder,filename(parID)), "wb" ) )
pickle.dump(U_record, open(modellingpath + '/3954/paper/out/numerical/colonies/simulation/%s/2Drecord_%s.pkl'%(folder,filename(parID)), 'wb'))

print(modellingpath + '/3954/paper/out/numerical/colonies/simulation/%s/2Drecord_%s.pkl'%(folder,filename(parID)))
print('saved')

















