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
parID = 383859
circuit_n=2
variant= 3
n_species=6


nsamples =  1000000
save_figure = False
tqdm_disable = False #disable tqdm
# boundarycoeff = float(sys.argv[6])
boundarycoeff = 1.7
p_division=0.5;seed=1

# open parameter dictionaries
df= pickle.load( open(modellingpath + '/3954/paper/input/lhs_parameterfiles/df_circuit%s_variant%s_%rparametersets.pkl'%(circuit_n,variant,nsamples), "rb" ) )
# instabilities_df= pickle.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/instabilities_dataframes/instability_df_circuit%s_variant%r_%rparametersets.pkl'%(circuit_n,variant,nsamples), "rb" ) )


#solver parameters
L=8; dx =0.05; J = int(L/dx)
T =125; dt = 0.05; N = int(T/dt)

# L=8; dx =0.5; J = int(L/dx)
# T =125; dt = 0.5; N = int(T/dt)


cell_matrix_record = pickle.load( open(modellingpath + "/3954/paper/out/numerical/masks/caMask_seed%s_pdivision%s_L%s_J%s_T%s_N%s.pkl"%(seed,p_division,L,J,T,N), "rb" ) )
daughterToMotherDictList = pickle.load( open(modellingpath + "/3954/paper/out/numerical/masks/caMemory_seed%s_pdivision%s_L%s_J%s_T%s_N%s.pkl"%(seed,p_division,L,J,T,N), "rb" ) )
# T =50; dt = 0.1; N = int(T/dt)

filename= lambda parID: 'circuit%r_variant%s_bc%s_%s_ID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundarycoeff, shape,parID,L,J,T,N)


# for parID in range(10):
par_dict = df.iloc[parID].to_dict()

D = np.zeros(n_species)
D[:2] = [par_dict['DA'],par_dict['DB'] ]

print(par_dict)

# U_record,U_final =  adi_ca_openclosed_nodilution_preMask(par_dict,L,dx,J,T,dt,N, circuit_n, n_species,D,cell_matrix_record, daughterToMotherDictList,tqdm_disable=False, p_division=0.5,stochasticity=0, seed=1,growth='Slow', boundarycoeff=boundarycoeff)
# U_record,U_final =  adi(par_dict,L,L,J,J,T,N, circuit_n, n_species,D,tqdm_disable=False,stochasticity=0, steadystates=0)
# get the start time
st = time.time()
U_record,U_final =  adi_ca_openclosed_nodilution_preMask_numba(par_dict,L,dx,J,T,dt,N, circuit_n, n_species,D,cell_matrix_record, daughterToMotherDictList,tqdm_disable=False, p_division=0.5,stochasticity=0, seed=1,growth='Slow', boundarycoeff=boundarycoeff)
elapsed_time = time.time() - st
print('Execution time numba:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
plt.imshow(U_final[0])
plt.show()
print(U_final[0,12,12])

st = time.time()
U_record,U_final =  adi_ca_openclosed_nodilution_preMask(par_dict,L,dx,J,T,dt,N, circuit_n, n_species,D,cell_matrix_record, daughterToMotherDictList,tqdm_disable=False, p_division=0.5,stochasticity=0, seed=1,growth='Slow', boundarycoeff=boundarycoeff)
elapsed_time = time.time() - st
print('Execution time no numba:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))




print(U_final[0,12,12])
plt.imshow(U_final[0])
plt.show()

# plt.imshow(U_final[0])
# plt.show()
# pickle.dump(U_final, open(modellingpath + '/3954/paper/out/numerical/colonies/simulation/2Dfinal_%s.pkl'%filename(parID), "wb" ) )




# savefigpath = modellingpath + '/3954/paper/out/numerical/colonies/figures/'
# plot_redgreen_contrast(U_final,mm=L)





























