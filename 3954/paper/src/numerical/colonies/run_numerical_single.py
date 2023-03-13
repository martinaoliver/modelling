#%%
#############
###paths#####
#############
import os
import sys


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
from numerical.cn_plot import *
from colonyMaskCreation import *
#%%
# # %matplotlib inline
#############
###execution parameters#####
#############
# %matplotlib inline
shape = 'ca'
circuit_n=14;variant='2nd';n_species=6
# Specifiy number of parameter sets in parameterset file to be loaded
n_param_sets = 10
# balance = 'balanced'
folder = 'circuit14variant2ndBalancedTuring'
nsamples =  1000000
save_figure = False
tqdm_disable = False #disable tqdm
# boundarycoeff = float(sys.argv[6])


# open parameter dictionaries


# df= pickle.load( open(modellingpath + '/3954/paper/input/balanced_parameterfiles/df_circuit%r_variant%s_%rparametersets_balanced.pkl'%(circuit_n,variant,nsamples), "rb" ) )
# df= pickle.load( open(modellingpath + '/3954/paper/input/fitted_parameterfiles/df_circuit%r_variant%s_%rparametersets.pkl'%(circuit_n,variant,nsamples), "rb" ) )
# instabilities_df= pickle.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/instabilities_dataframes/instability_df_circuit%s_variant%s_%rparametersets.pkl'%(circuit_n,variant,nsamples), "rb" ) )
with open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/turing_dataframes/turing_df_circuit%s_variant%s_%rparametersets_balanced.pkl'%(circuit_n,variant,nsamples), "rb" ) as f:
    df = pickle.load(f)
print(df)
#solver parameters
# specify dimensions of system
L=4; dx =0.01; J = int(L/dx)
T =50; dt = 0.005; N = int(T/dt)
boundarycoeff = 1

divisionTimeHours=0.5
p_division=0.3;seed=1
shape = 'ca'
x_gridpoints=int(1/dx)
# divisionTimeHours=0.5
# p_division=0.22;seed=1
# L=int(sys.argv[1]); dx =float(sys.argv[2]); J = int(L/dx)
# T =int(sys.argv[3]); dt = float(sys.argv[4]); N = int(T/dt)

try:
    cell_matrix_record = pickle.load( open(modellingpath + "/3954/paper/out/numerical/masks/caMask_seed%s_pdivision%s_L%s_J%s_T%s_N%s.pkl"%(seed,p_division,L,J,T,N), "rb" ) )
    daughterToMotherDictList = pickle.load( open(modellingpath + "/3954/paper/out/numerical/masks/caMemory_seed%s_pdivision%s_L%s_J%s_T%s_N%s.pkl"%(seed,p_division,L,J,T,N), "rb" ) )
except:
    #file does not exist
    FileNotFoundError
    maskFunction(L=L,dx=dx, T=T, dt=dt, divisionTimeHours=divisionTimeHours, p_division=p_division, plot1D=True, plotScatter=True)
    cell_matrix_record = pickle.load( open(modellingpath + "/3954/paper/out/numerical/masks/caMask_seed%s_pdivision%s_L%s_J%s_T%s_N%s.pkl"%(seed,p_division,L,J,T,N), "rb" ) )
    daughterToMotherDictList = pickle.load( open(modellingpath + "/3954/paper/out/numerical/masks/caMemory_seed%s_pdivision%s_L%s_J%s_T%s_N%s.pkl"%(seed,p_division,L,J,T,N), "rb" ) )

    
    # T =1; dt = 0.05; N = int(T/dt)
# 0.5,0.02, 0.005
# T =2; dt = 0.05; N = int(T/dt)
filename= lambda parID: 'circuit%r_variant%s_bc%s_%s_ID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundarycoeff, shape,parID,L,J,T,N)
#%%

parID=316386
print('parID = ' + str(parID))
par_dict = df.loc[parID].to_dict()
D = np.zeros(n_species)
Dr = float(par_dict['Dr'])
D[:2] = [1,Dr ]
degDiv = 1
par_dict['muASV'] =par_dict['muASV']/degDiv
par_dict['muLVA'] = par_dict['muLVA'] /degDiv
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



pickle.dump(U_final, open(modellingpath + '/3954/paper/out/numerical/colonies/simulation/%s/2Dfinal_%s_degDiv%r.pkl'%(folder,filename(parID), degDiv), "wb" ) )
pickle.dump(U_record, open(modellingpath + '/3954/paper/out/numerical/colonies/simulation/%s/2Drecord_%s_degDiv%r.pkl'%(folder,filename(parID), degDiv), 'wb'))

print(modellingpath + '/3954/paper/out/numerical/colonies/simulation/%s/2Drecord_%s_degDiv%s.pkl'%(folder,filename(parID), degDiv))
print('saved')
# %%
rgb = plot_redgreen_contrast(U_final,L,parID=parID,scale_factor=x_gridpoints,save_figure=False)

# %%
plt.imshow(U_final[-1])
plt.colorbar()
# %%
