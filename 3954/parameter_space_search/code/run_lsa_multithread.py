#IMPORTS#
#############################
import sys
import os

from joblib import parallel_backend
pwd = os.getcwd()
root = pwd.rpartition("mo2016")[0] + pwd.rpartition("mo2016")[1] #/Volumes/mo2016/ or '/Users/mo2016/' or '/rds/general/mo2016/'
print(root)
if root == '/Users/mo2016':
    modelling_ephemeral = '/Volumes/mo2016/ephemeral/Documents/modelling'
    modelling_home = '/Volumes/mo2016/home/Documents/modelling'
    modelling_local = root + '/Documents/modelling'
    import matplotlib as mpl
    mpl.use('tkagg')

if root == '/Volumes/mo2016' or root=='/rds/general/user/mo2016': #'/rds/general' or root=='/Volumes':
        modelling_ephemeral = root + '/ephemeral/Documents/modelling'
        modelling_home = root  + '/home/Documents/modelling'
        modelling_local = modelling_home

modulepath = modelling_local + '/3954/modules'

sys.path.append(modulepath)

from linear_stability_analysis import *


# Change directory, useful for testing on local machine, as script needs to find other files
#path = "C:\\Users\\anna_\\OneDrive - Imperial College London\\Computing project\\Code\\Global_param_search\\UPDATED_GLobalLSA"
#os.chdir(path)

import sys
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
circuit_n=2
variant= 9

# Specifiy number of parameter sets in parameterset file to be loaded
df_lenght = 1000000
n_param_sets = 1000000
# n_param_sets = 10000


# Specify date today
date = date.today().strftime('%m_%d_%Y')

# Specify size of batches in which to complete computations
# Does not need to be a factor of number of parameter sets
# batch_size = 20000
batch_size = 20833
# batch_size = 2

# Get starting parameter set indices for batches


# Define work to be done per batch of parameter sets
def lsa_check(start_batch_index,n_param_sets,df,circuit_n=2, variant=9, n_species=6):
    print('pool' + str(start_batch_index))
    output_df = big_turing_analysis_df(df,circuit_n,n_species,print_parID=False)
    print('calculated')
    pickle.dump(output_df, open('../results/output_dataframes/lsa_df_circuit%r_variant%r_%rparametersets_batch%r.pkl'%(circuit_n,variant,n_param_sets,start_batch_index), 'wb'))
    # pickle.dump(output_df, open('../results/output_dataframes/lsa_df_circuit%r_variant%r_%rparametersets_batch%r_rbslibrary0.pkl'%(circuit_n,variant,n_param_sets,start_batch_index), 'wb'))
    print('saved')
# Runs if the module is the main program
# if __name__ == '__main__':

start_time = time.perf_counter()
start_parameter = int(0)
# Load dataframe of parameter sets
df= pickle.load( open('../parameterfiles/df_circuit%r_variant%r_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), "rb" ) )
# df= pickle.load( open("../parameterfiles/df_circuit2_variant1_1954parametersets_rbslibrary0.pkl", "rb"))
batch_indices = list(range(0+start_parameter, len(df) + start_parameter, batch_size))
# batch_indices = list(range(0+start_parameter, 10 + start_parameter, batch_size))

# Create a pool of workers
pool = multiprocessing.Pool(Number_of_Threads)

# Define jobs as different batches of parameter sets
# Run lsa_check function in parallel across different threads
pool_output = []
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
    my_data[start_batch_index] = pickle.load(open('../results/output_dataframes/lsa_df_circuit%r_variant%r_%rparametersets_batch%r.pkl'%(circuit_n,variant,n_param_sets,start_batch_index), "rb" ) )
    # my_data[start_batch_index] = pickle.load(open('../results/output_dataframes/lsa_df_circuit%r_variant%r_%rparametersets_batch%r_rbslibrary0.pkl'%(circuit_n,variant,n_param_sets,start_batch_index), "rb" ) )
# Join all batch results to large results dataframe
results_df = pd.concat(my_data.values(), ignore_index=False)

# Pickle and save results dataframe
tupled_index =  [tuple(l) for l in results_df.index]
multi_index = pd.MultiIndex.from_tuples(tupled_index)
results_df = results_df.set_index(multi_index)
pickle.dump(results_df, open('../results/output_dataframes/lsa_df_circuit%r_variant%r_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), 'wb'))
# pickle.dump(results_df, open('../results/output_dataframes/lsa_df_circuit%r_variant%r_%rparametersets_rbslibrary0_concat.pkl'%(circuit_n,variant,n_param_sets), 'wb'))
