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

from analytical.linear_stability_analysis import big_turing_analysis_df, detailed_turing_analysis_dict


import pickle
from datetime import date
import pandas as pd
import numpy as np
import time
import multiprocessing


'''
====================================================
    Code
====================================================
'''

# Set number of threads to 1 if no valid number provided
if len(sys.argv) > 1:
    Number_of_Threads = int(sys.argv[1])
else:
    Number_of_Threads = 1
print('Number of Threads set to ', Number_of_Threads)
# Number_of_Threads=48
# Specify name of circuit and variant investigated

circuit_n='circuit14'
# variant='fitted7'
variant='2nd'
n_species=6
# balance='notBalanced'
balance=str(sys.argv[2])
# Specifiy number of parameter sets in parameterset file to be loaded
n_samples =1000000
n_analysed_param_sets = 1000000
print(n_analysed_param_sets)
# Specify date today
date = date.today().strftime('%m_%d_%Y')

# Specify size of batches in which to complete computations
# Does not need to be a factor of number of parameter sets
# batch_n_samples =10000000
# batch_size = 2
# batch_size =5
# batch_size =int(n_samples/Number_of_Threads)
# print(f'batch_size: {batch_size}')


# Define work to be done per batch of parameter sets
def lsa_check(start_batch_index,n_samples,df,circuit_n=circuit_n, variant=variant, n_species=n_species):
    print('pool' + str(start_batch_index))
    output_df = big_turing_analysis_df(df,circuit_n,variant, n_samples, n_species,print_parID=False, saveInstability = False)
    print('calculated')
    
    pickle.dump(output_df, open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/all_dataframes/lsa_df_%s_variant%s_%rparametersets_%s_batch%r.pkl'%(circuit_n,variant,n_samples,balance, start_batch_index), 'wb'))
    # pickle.dump(output_df, open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/all_dataframes/lsa_df_%s_variant%s_%rparametersets_%s_batch%r.pkl'%(circuit_n,variant,n_samples,balance, start_batch_index), 'wb'))
    print('saved')
# Runs if the module is the main program
# if __name__ == '__main__':
print('start_time')
start_time = time.perf_counter()
start_parameter = int(0)
# Load dataframe of parameter sets
# df= pickle.load( open('../lhs_parameterfiles/df_circuit%r_variant%s_%rparametersets.pkl'%(circuit_n,variant,n_samples), "rb" ) )
df= pickle.load( open(modellingpath + "/3954/paper/input/balanced_parameterfiles/df_%s_variant%s_%rparametersets_%s.pkl"%(circuit_n,variant,n_samples, balance), "rb"))

# df= pickle.load( open(modellingpath + "/3954/paper/input/fitted_parameterfiles/df_%s_variant%s_%rparametersets_%s.pkl"%(circuit_n,variant,n_samples,balance), "rb"))
print(df)
print('df_loaded')

# n_analysed_param_sets = len(df)
# print(n_analysed_param_sets)
df = df.iloc[start_parameter:start_parameter+n_analysed_param_sets]

batch_size =int(n_analysed_param_sets/Number_of_Threads)
batch_indices = list(range(0+start_parameter, n_analysed_param_sets + start_parameter, batch_size))

# Create a pool of workers
pool = multiprocessing.Pool(Number_of_Threads)

# Define jobs as different batches of parameter sets
# Run lsa_check function in parallel across different threads
pool_output = []
print('start_loop')
for start_batch_index in batch_indices:

    print('main' + str(start_batch_index))
    df_batch = df.iloc[start_batch_index:start_batch_index+batch_size]
    # df_batch = df.iloc[start_batch_index:start_batch_index+2]

    pool_output.append(pool.apply_async(lsa_check, args=(start_batch_index, n_samples,df_batch)))
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

my_data = {}

# Load all batch dataframes
for start_batch_index in batch_indices:
    try:
        my_data[start_batch_index] = pickle.load(open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/all_dataframes/lsa_df_%s_variant%s_%rparametersets_%s_batch%r.pkl'%(circuit_n,variant,n_samples,balance,start_batch_index), "rb" ) )
        print(start_batch_index)
    except FileNotFoundError:
        print(start_batch_index, 'not Found')
    # my_data[start_batch_index] = pickle.load(open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/all_dataframes/lsa_df_%s_variant%s_%rparametersets_%s_batch%r.pkl'%(circuit_n,variant,n_samples,balance,start_batch_index), "rb" ) )
    # my_data[start_batch_index] = pickle.load(open('../results/output_dataframes/lsa_df_circuit%r_variant%s_%rparametersets_batch%r_rbslibrary0.pkl'%(circuit_n,variant,n_samples,start_batch_index), "rb" ) )
# Join all batch results to large results dataframe
results_df = pd.concat(my_data.values(), ignore_index=False)

# Pickle and save results dataframe
tupled_index =  [tuple(l) for l in results_df.index]
multi_index = pd.MultiIndex.from_tuples(tupled_index)
results_df = results_df.set_index(multi_index)
# pickle.dump(results_df, open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/all_dataframes/lsa_df_%s_variant%s_%rparametersets_balancedSemiBalanced.pkl'%(circuit_n,variant,n_samples), 'wb'))
pickle.dump(results_df, open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/all_dataframes/lsa_df_%s_variant%s_%rparametersets_%s.pkl'%(circuit_n,variant,n_samples,balance), 'wb'))
