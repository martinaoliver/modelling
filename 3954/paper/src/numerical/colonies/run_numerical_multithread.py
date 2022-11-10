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

from numerical.adi_ca_function_openclosed_nodilution_preMask import adi_ca_openclosed_nodilution_preMask
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
circuit_n=2
# variant= '48257gaussian0.21nsr'
variant=3
n_species=6
folder = 'circuit2variant3_turing'
# folder = 'circuit2variant48257gaussian0.21nsr'
# Specifiy number of parameter sets in parameterset file to be loaded
# nsamples = 2000
nsamples = 1000000


# Specify date today
date = date.today().strftime('%m_%d_%Y')



def numerical_check(df,circuit_n, variant = variant, n_species=n_species, folder=folder, test=False, nsamples=nsamples, date=date):
    L=8; dx =0.05; J = int(L/dx)
    T =125; dt = 0.05; N = int(T/dt)
    boundarycoeff = 1.7
    p_division=0.5;seed=1
    df_index = np.unique(df.index.get_level_values(0))

    cell_matrix_record = pickle.load( open(modellingpath + "/3954/paper/out/numerical/masks/caMask_seed%s_pdivision%s_L%s_J%s_T%s_N%s.pkl"%(seed,p_division,L,J,T,N), "rb" ) )
    daughterToMotherDictList = pickle.load( open(modellingpath + "/3954/paper/out/numerical/masks/caMemory_seed%s_pdivision%s_L%s_J%s_T%s_N%s.pkl"%(seed,p_division,L,J,T,N), "rb" ) )
    if test==True:
        T =1; dt = 0.05; N = int(T/dt)
    D = np.zeros(n_species)


    for parID in df_index:
        print('parID = ' + str(parID))
        shape = 'ca'


        par_dict = df.loc[parID].to_dict()
        print(par_dict)
        D[:2] = [par_dict['DA'],par_dict['DB'] ]

        # steadystates=par_dict['ss_list']

        filename = 'circuit%r_variant%s_bc%s_%s_ID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundarycoeff, shape,parID,L,J,T,N)
        savefig=False
        savefigpath = modellingpath + '/3954/paper/out/numerical/colonies/figures/%s/'%(folder)

        try:
            U_record,U_final =  adi_ca_openclosed_nodilution_preMask(par_dict,L,dx,J,T,dt,N, circuit_n, n_species,D,cell_matrix_record, daughterToMotherDictList,tqdm_disable=True, p_division=p_division,stochasticity=0, seed=seed,growth='Slow', boundarycoeff=boundarycoeff)
            pickle.dump(U_final, open(modellingpath + '/3954/paper/out/numerical/colonies/simulation/%s/2Dfinal_%s.pkl'%(folder,filename), "wb" ) )
            # pickle.dump(U_record, open(modellingpath + '/growth/out/numerical/%s/%s/data/2Drecord_%s.pkl'%(circuit_n,mechanism,filename), 'wb'))
            # print(np.shape(U_record))
            if savefig==True:
                plot1D(U_final, savefig=True,filename=filename, savefigpath=savefigpath)
                plt.close()
        except ValueError:
            print('!!!!!!!!!!!!!')
            print('ValueError --> unstable solution')
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
turing_df= pickle.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/turing_dataframes/turing_df_circuit%s_variant%s_%rparametersets.pkl'%(circuit_n,variant,nsamples), "rb" ) )
df = turing_df
total_params=len(df)
# total_params=10
print(total_params)
print('loaded')
batch_size = int(total_params/Number_of_Threads) + 1
instabilities_df = df.iloc[0:total_params]
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

for count,start_batch_index in enumerate(batch_indices):
    print('error' + str(start_batch_index))
    pool_output[count].get()
# Report time taken
finish_time = time.perf_counter()
time_taken = finish_time-start_time
print("Time taken: %d s" %time_taken)
