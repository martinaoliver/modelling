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

from numerical.edgegrowth2 import cn_edgegrowth2
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
circuit_n='turinghill'
variant= 0
n_species=2

# Specifiy number of parameter sets in parameterset file to be loaded
n_param_sets = 2000000



# Specify date today
date = date.today().strftime('%m_%d_%Y')
# Specify size of batches in which to complete computations
# Does not need to be a factor of number of parameter sets
total_params=10



def numerical_check(df,x_gridpoints, t_gridpoints,T,L,circuit_n=circuit_n, variant = variant, n_species=n_species):

    df_index = np.unique(df.index.get_level_values(0))
    for parID in df_index:
        print('parID = ' + str(parID))
        mechanism = 'nogrowth'


        par_dict = df.loc[parID].to_dict()


        J = L *x_gridpoints  # number of equally spaced gridpoints in space domain (larger J means more spatial precision(tends towards continuum solution) )

        N = T * t_gridpoints
        steadystates=par_dict['ss_list']

        filename = '%s_variant%s_%s_ID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant,mechanism,parID,L,J,T,N)
        savefig=False
        savefigpath = modellingpath + '/growth/out/numerical/%s/%s/fig/'%(circuit_n,mechanism)

        try:
            U_final,U_record, U0, x_grid, reduced_t_grid, cellMatrix= cn_edgegrowth2(par_dict,L,J,T,N, circuit_n, tqdm_disable=True)            
            pickle.dump(U_final, open(modellingpath + '/growth/out/numerical/%s/%s/data/2Dfinal_%s.pkl'%(circuit_n,mechanism,filename), 'wb'))
            pickle.dump(U_record, open(modellingpath + '/growth/out/numerical/%s/%s/data/2Drecord_%s.pkl'%(circuit_n,mechanism,filename), 'wb'))
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

L=50; x_gridpoints=5
T=5000; t_gridpoints = 30


# L=2; x_gridpoints=2
# T=100; t_gridpoints = 2

# Load dataframe of parameter sets
# multiple_df= pickle.load( open(modellingpath + "/growth/out/analytical/lsa_dataframes/lsa_df_%s_variant%r_%rparametersets.pkl"%(circuit_n,variant,n_param_sets), "rb"))
# multiple_df= pickle.load( open(modellingpath + "/growth/input/parameterfiles/df_%s_variant%r_%rparametersets.pkl"%(circuit_n,variant,n_param_sets), "rb"))
# df= pickle.load( open(modellingpath + "/growth/input/parameterfiles/df_%s_variant%r_%rparametersets.pkl"%(circuit_n,variant,n_param_sets), "rb"))
# df = multiple_df.xs(0, level=1)
df= pickle.load( open(modellingpath + '/growth/out/analytical/instability/instability_df_%s_variant%r_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), "rb"))
total_params=len(df)
print(total_params)
print('loaded')
batch_size = int(total_params/Number_of_Threads) + 1
df = df.iloc[0:total_params]
print(df.head())
batch_indices = list(range(0, len(df), batch_size))
# Create a pool of workers
pool = multiprocessing.Pool(Number_of_Threads)

# Run lsa_check function in parallel across different threads
pool_output = []
for start_batch_index in batch_indices:

    print('main' + str(start_batch_index))
    df_batch = df.iloc[start_batch_index:start_batch_index+batch_size]
    pool_output.append(pool.apply_async(numerical_check, args=(df_batch,x_gridpoints, t_gridpoints,T,L)))

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