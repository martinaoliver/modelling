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
    number_of_cpus = int(sys.argv[1])
else:
    Number_of_cpus = 1
print('Number of Threads set to ', number_of_cpus)
number_of_threads=100000
# Number_of_Threads=48
# Specify name of circuit and variant investigated
circuit_n='turinghill'
variant= int(sys.argv[2])
n_species=2
# Specifiy number of parameter sets to analyse
n_param_sets = 25000000 #so you can find 1000 hopf if there is 8 hopf in 2 million 
batch_size =int(n_param_sets/number_of_threads)
print(batch_size)


# Specify date today
date = date.today().strftime('%m_%d_%Y')





# Define work to be done per batch of parameter sets
def lsa_check(thread_number,n_param_sets=batch_size,circuit_n=circuit_n, variant=variant, n_species=n_species):
    #call function to create df(seed)
    df_batch = create_df_function(thread_number, n_param_sets)

    #calculate df lsa
    output_df = big_turing_analysis_df(df_batch,circuit_n,n_species,print_parID=False)  

   
    #query how many system class in db. 

    if system_class_x_n < 1000:
        #  add system_class_x_n results onto db 
            #add general robustness results onto db 
            #  pickle.dump('dict with results')
        search_finished=False
    
    else:
        search_finished = True


    return search_finished
# Runs if the module is the main program
# if __name__ == '__main__':
print('start_time')
start_time = time.perf_counter()


# Create a pool of workers
pool = multiprocessing.Pool(number_of_cpus)

# Define jobs as different batches of parameter sets
# Run lsa_check function in parallel across different threads
pool_output = []
print('start_loop')
search_finished=False
for thread_number in range(number_of_threads):
    while search_finished==False:
        print('main' + str(thread_number))
        search_finished = pool_output.append(pool.apply_async(lsa_check, args=(thread_number)))
        


# Close the parallel processing job
pool.close()
pool.join()
print('Run finished')

for count,start_batch_index in enumerate(number_of_threads):
    print('error' + str(start_batch_index))
    pool_output[count].get()
# Report time taken
finish_time = time.perf_counter()
time_taken = finish_time-start_time
time_taken_hours = time_taken/3600
print(f'Time taken: {time_taken} s, {time_taken_hours} h')

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