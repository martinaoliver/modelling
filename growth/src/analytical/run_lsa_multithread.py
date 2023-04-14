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
circuit_n='turinghill'
variant= 4
n_species=2
# Specifiy number of parameter sets in parameterset file to be loaded
df_lenght = 2000000
n_param_sets = 2000000
# df_lenght = 10
# n_param_sets = 10


# Specify date today
date = date.today().strftime('%m_%d_%Y')

# Specify size of batches in which to complete computations
# Does not need to be a factor of number of parameter sets
# batch_size = 20000
batch_size = 41666
# batch_size = 10
print(batch_size)


# Define work to be done per batch of parameter sets
def lsa_check(start_batch_index,n_param_sets,df,circuit_n=circuit_n, variant=variant, n_species=n_species):
    print('pool' + str(start_batch_index))
    output_df = big_turing_analysis_df(df,circuit_n,n_species,print_parID=False)
    print('calculated')
    pickle.dump(output_df, open(modellingpath + '/growth/out/analytical/lsa_dataframes/lsa_df_%s_variant%r_%rparametersets_batch%r.pkl'%(circuit_n,variant,n_param_sets,start_batch_index), 'wb'))
    print('saved')
# Runs if the module is the main program
# if __name__ == '__main__':
print('start_time')
start_time = time.perf_counter()
start_parameter = int(0)
# Load dataframe of parameter sets
print('df_%s_variant%r_%rparametersets.pkl'%(circuit_n,variant,n_param_sets))
# df= pickle.load( open('../parameterfiles/df_circuit%r_variant%r_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), "rb" ) )
df= pickle.load( open(modellingpath + "/growth/input/parameterfiles/df_%s_variant%r_%rparametersets.pkl"%(circuit_n,variant,n_param_sets), "rb"))
# df = df.iloc[:10]
# df= pickle.load( open("../parameterfiles/df_circuit2_variant1_1954parametersets_rbslibrary0.pkl", "rb"))
batch_indices = list(range(0+start_parameter, len(df) + start_parameter, batch_size))
# batch_indices = list(range(0+start_parameter, 10 + start_parameter, batch_size))

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

    pool_output.append(pool.apply_async(lsa_check, args=(start_batch_index, n_param_sets,df_batch)))
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
    my_data[start_batch_index] = pickle.load(open(modellingpath + '/growth/out/analytical/lsa_dataframes/lsa_df_%s_variant%r_%rparametersets_batch%r.pkl'%(circuit_n,variant,n_param_sets,start_batch_index), "rb" ) )
    # my_data[start_batch_index] = pickle.load(open('../results/output_dataframes/lsa_df_circuit%r_variant%r_%rparametersets_batch%r_rbslibrary0.pkl'%(circuit_n,variant,n_param_sets,start_batch_index), "rb" ) )
# Join all batch results to large results dataframe
results_df = pd.concat(my_data.values(), ignore_index=False)

# Pickle and save results dataframe
tupled_index =  [tuple(l) for l in results_df.index]
multi_index = pd.MultiIndex.from_tuples(tupled_index)
results_df = results_df.set_index(multi_index)
pickle.dump(results_df, open(modellingpath + '/growth/out/analytical/lsa_dataframes/lsa_df_%s_variant%r_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), 'wb'))
# print(results_df)