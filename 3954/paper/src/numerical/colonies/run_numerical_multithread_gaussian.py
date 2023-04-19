#############
###paths#####
#############
import sys
import os

pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
modellingephemeral = '/rds/general/ephemeral/user/mo2016/ephemeral/Documents/modelling'
sys.path.append(modellingpath + '/lib')
#############

from numerical.adi_ca_function_openclosed_nodilution_preMask_numba import \
    adi_ca_openclosed_nodilution_preMask as \
    adi_ca_openclosed_nodilution_preMask_numba
from numerical.adi_square_function import adi
from numerical.plotting_numerical import *

from numerical.cn_plot import plot1D, surfpattern


import pickle
from datetime import date
import pandas as pd
import numpy as np
import time
import multiprocessing
import matplotlib.pyplot as plt
'''
====================================================
    Code
====================================================
'''
print('f')
# Set number of threads to 1 if no valid number provided
if len(sys.argv) > 1:
    Number_of_Threads = int(sys.argv[1])
else:
    Number_of_Threads = 1
print('Number of Threads set to ', Number_of_Threads)

# Specify name of circuit and variant investigated
# circuit_n=14;variant='fitted1';n_species=6
circuit_n=14;variant=int(sys.argv[2]);n_species=6;nsr=float(sys.argv[3])

# Specifiy number of parameter sets in parameterset file to be loaded
# balance = 'balanced'
# folder = 'circuit14variantfitted1'
# folder = 'circuit14variant2ndBalancedTuring'
folder = f'circuit14variant{variant}'
modelArgs = [circuit_n,variant,n_species,folder]

# Specifiy number of parameter sets in parameterset file to be loaded
# nsamples = 2000000
nsamples = 2000
# specify dimensions of system

L=20; dx =0.1; J = int(L/dx)
T =50; dt = 0.02; N = int(T/dt)
boundarycoeff = 1
divisionTimeHours=0.5
p_division=1;seed=1



systemArgs = [L, dx, J, T, dt, N, boundarycoeff, p_division, seed, divisionTimeHours]




# Specify date today
date = date.today().strftime('%m_%d_%Y')
# cell_matrix_record = pickle.load( open(modellingpath + "/3954/paper/out/numerical/masks/caMask_seed%s_pdivision%s_L%s_J%s_T%s_N%s.pkl"%(seed,p_division,L,J,T,N), "rb" ) )
# daughterToMotherDictList = pickle.load( open(modellingpath + "/3954/paper/out/numerical/masks/caMemory_seed%s_pdivision%s_L%s_J%s_T%s_N%s.pkl"%(seed,p_division,L,J,T,N), "rb" ) )
with open(modellingpath + "/3954/paper/out/numerical/masks/caMask_seed%s_pdivision%s_L%s_J%s_T%s_N%s.pkl"%(seed,p_division,L,J,T,N), "rb" ) as f:
    cell_matrix_record = pickle.load(f)
with open(modellingpath + "/3954/paper/out/numerical/masks/caMemory_seed%s_pdivision%s_L%s_J%s_T%s_N%s.pkl"%(seed,p_division,L,J,T,N), "rb" ) as f:
    daughterToMotherDictList = pickle.load(f)


def numerical_check(df, circuit_n,modelArgs=modelArgs, systemArgs=systemArgs,cell_matrix_record = cell_matrix_record,daughterToMotherDictList=daughterToMotherDictList, variant = variant, n_species=n_species, folder=folder):
    # L=8; dx =0.02; J = int(L/dx)
    # T =125; dt = 0.05; N = int(T/dt)
    # boundarycoeff = 1.7
    # p_division=0.7;seed=1
    # divisionTimeHours = 0.5
    circuit_n, variant, n_species, folder = modelArgs
    L, dx, J, T, dt, N, boundarycoeff, p_division, seed, divisionTimeHours = systemArgs
    df_index = np.unique(df.index.get_level_values(0))

    if int(sys.argv[1]) == 1:
        T =1; dt =0.5; N = int(T/dt)
        tqdm_disable = False
    else:
        tqdm_disable = True
    print(systemArgs)
    shape = 'ca'

    for parID in df_index:
        
        print('parID = ' + str(parID))
        par_dict = df.loc[parID].to_dict()
        D = np.zeros(n_species)
        Dr = float(par_dict['Dr'])
        D[:2] = [1,Dr ]
        # par_dict['muASV'] =par_dict['muASV']/degDiv
        # par_dict['muLVA'] = par_dict['muLVA'] /degDiv
        # steadystates=par_dict['ss_list']

        filename= lambda parID: 'circuit%r_variant%snsr%s_bc%s_%s_ID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant,nsr,boundarycoeff, shape,parID,L,J,T,N)
        savefig=False
        savefigpath = modellingpath + '/3954/paper/out/numerical/colonies/figures/%s/'%(folder)

        try:
            U_record,U_final =  adi_ca_openclosed_nodilution_preMask_numba(par_dict,L,dx,J,T,dt,N, circuit_n, n_species,D,cell_matrix_record, daughterToMotherDictList,tqdm_disable=tqdm_disable,divisionTimeHours=divisionTimeHours, stochasticity=0, seed=1, boundarycoeff=boundarycoeff)

            # pickle.dump(U_final, open(modellingpath + '/3954/paper/out/numerical/colonies/simulation/%s/2Dfinal_%s.pkl'%(folder,filename), "wb" ) )
            # pickle.dump(U_record, open(modellingpath + '/3954/paper/out/numerical/colonies/simulation/%s/2Drecord_%s.pkl'%(folder,filename), 'wb'))

            
            save = True
            if save == True:
                        
                with open(modellingephemeral + '/3954/paper/out/numerical/colonies/simulation/%s/2Dfinal_%s.pkl'%(folder,filename(parID)), "wb" ) as f:
                    pickle.dump(U_final, f)
                # with open(modellingephemeral + '/3954/paper/out/numerical/colonies/simulation/%s/2Drecord_%s.pkl'%(folder,filename), "wb" ) as f:
                with open(modellingephemeral + '/3954/paper/out/numerical/colonies/simulation/%s/2Drecord_%s.pkl'%(folder,filename(parID)), "wb" ) as f:
                    pickle.dump(U_record, f)
                
            del U_record
            del U_final
            # print(np.shape(U_record))
            if savefig==True:
                plot1D(U_final, savefig=True,filename=filename, savefigpath=savefigpath)
                plt.close()
        except ValueError:
            print('!!!!!!!!!!!!!')
            print(f'ValueError --> unstable solution in {parID}')
            print('!!!!!!!!!!!!!')
            print()

            pass

# Runs if the module is the main program
# if __name__ == '__main__':

start_time = time.perf_counter()

#parameters
# L=int(sys.argv[2]); x_gridpoints = int(sys.argv[3])
# T =int(sys.argv[4]); t_gridpoints = int(sys.argv[5]) 


# Load dataframe of parameter sets
# df= pickle.load( open(modellingpath + '/3954/paper/input/gaussian_parameterfiles/df_circuit%s_variant%s_%rparametersets.pkl'%(circuit_n,variant,nsamples), "rb" ) )
# df= pickle.load( open(modellingpath + '/3954/paper/input/lhs_parameterfiles/df_circuit%s_variant%s_%rparametersets.pkl'%(circuit_n,variant,nsamples), "rb" ) )
# instabilities_df= pickle.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/instabilities_dataframes/instability_df_circuit%s_variant%r_%rparametersets.pkl'%(circuit_n,variant,nsamples), "rb" ) )
# instabilities_df= pickle.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/instabilities_dataframes/instability_df_circuit%s_variant%r_%rparametersets.pkl'%(circuit_n,variant,nsamples), "rb" ) )
# with open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/instabilities_dataframes/instability_df_circuit%s_variant%s_%rparametersets.pkl'%(circuit_n,variant,nsamples), "rb" ) as f:
# with open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/turing_dataframes/turing_df_circuit%s_variant%s_%sparametersets_balanced.pkl'%(circuit_n,variant,nsamples), "rb" ) as f:
#     df = pickle.load(f)
with open(modellingpath + '/3954/paper/input/gaussian_parameterfiles/df_circuit%r_variant%sgaussian%snsr_%rparametersets.pkl'%(circuit_n,variant,nsr,nsamples), "rb" ) as f:
    df = pickle.load(f)


# turing_df= pickle.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/turing_dataframes/turing_df_circuit%s_variant%s_%rparametersets.pkl'%(circuit_n,variant,nsamples), "rb" ) )
total_params=len(df)
# total_params=10
print(total_params)
print('loaded')
batch_size = int(total_params/Number_of_Threads) + 1
instabilities_df = df.iloc[0:total_params]
# df = df.iloc[85:total_params]
print(df.head())
batch_indices = list(range(0, len(df), batch_size))
# Create a pool of workers
pool = multiprocessing.Pool(Number_of_Threads)

# Run lsa_check function in parallel across different threads
pool_output = []
for start_batch_index in batch_indices:

    print('main' + str(start_batch_index))
    df_batch = df.iloc[start_batch_index:start_batch_index+batch_size]
    pool_output.append(pool.apply_async(numerical_check, args=(df_batch, circuit_n)))

# Close the parallel processing job
pool.close()
pool.join()
print('Run finished')
# Report time taken
finish_time = time.perf_counter()
time_taken = finish_time-start_time
print("Time taken: %d s" %time_taken)

for count,start_batch_index in enumerate(batch_indices):
    print('error' + str(start_batch_index))
    pool_output[count].get()
