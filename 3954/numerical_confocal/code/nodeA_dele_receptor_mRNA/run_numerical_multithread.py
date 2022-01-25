'''
====================================================
    Imports
====================================================
'''

import sys
modelling_hpc = '/rds/general/user/mo2016/home/Documents/modelling'
modulepath = modelling_hpc + '/3954/modules/new_CN'
modelling_ephemeral = '/rds/general/user/mo2016/ephemeral/Documents/modelling'

sys.path.append(modulepath)
from adi_v1_openclosed_ca_function import *
from plotting_numerical import *

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

# Specify name of circuit and variant investigated
circuit_n=11
variant= 8
n_species = 7
# Specifiy number of parameter sets in parameterset file to be loaded
n_param_sets = 30000


# Specify date today
date = date.today().strftime('%m_%d_%Y')
# Specify size of batches in which to complete computations
# Does not need to be a factor of number of parameter sets
# batch_size = 1000

# Get starting parameter set indices for batches


# Define work to be done per batch of parameter sets


def numerical_check(start_batch_index,n_param_sets,df,x_gridpoints,t_gridpoints,T,L,circuit_n=11, variant = 8, n_species=7,p_division=0.5,seed=1):
    save_figure = True
    tqdm_disable = True #disable tqdm

    df_index = np.unique(df.index.get_level_values(0))
    for parID in df_index:
        print('parID = ' + str(parID))
        mechanism = 'nodeA_dele_receptor_mRNA'
        boundary_coef = 1 #1 is open boundary and 0 is closed boundary
        shape = 'ca'
        growth = False

        # x_gridpoints = int(sys.argv[3])
        # T =int(sys.argv[4])
        par_dict = df.loc[parID].to_dict()


        par_dict['bmB']=1
        par_dict['bmE']=1
        par_dict['bmF']=1

        d_H = par_dict['d_H']
        D = np.zeros(n_species)
        D[1]=d_H

        J = L *x_gridpoints  # number of equally spaced gridpoints in space domain (larger J means more spatial precision(tends towards continuum solution) )
        N = T * t_gridpoints
        initial_condition = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001]

        filename = 'circuit%r_variant%s_%s_%sID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant, shape,mechanism,parID,L,J,T,N)

        L_x = L
        L_y = L
        I = J
        try:
            # Run 2D simulation
            U_record,U_final = adi_ca(par_dict,L_x,L_y,J,I,T,N, circuit_n,n_species, D,tqdm_disable=tqdm_disable)#,p_division=p_division,seed=seed)
            savefig_path = modelling_ephemeral + '/3954/numerical_confocal/results/figures/1M_colony_ca/2D/nodeA_dele_receptor_mRNA'
            plot_redgreen_contrast(U_final,L_x,mechanism,shape,filename,savefig_path,parID=parID,scale_factor=x_gridpoints,save_figure=save_figure)
            # rgb_timeseries = redgreen_contrast_timeseries(records)
            # show_rgbvideo(rgb_timeseries)
            if save_figure ==True:
                pickle.dump(U_final, open(modelling_ephemeral + '/3954/numerical_confocal/results/simulation/1M_colony_ca/2D/nodeA_dele_receptor_mRNA/2Dfinal_%s.pkl'%filename, 'wb'))
                pickle.dump(U_record,open(modelling_ephemeral + '/3954/numerical_confocal/results/simulation/1M_colony_ca/2D/nodeA_dele_receptor_mRNA/2Dtimeseries_%s.pkl'%filename, 'wb'))

            # else:
            #     plt.show()
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
t_gridpoints = int(sys.argv[2])
x_gridpoints = int(sys.argv[3])
T =int(sys.argv[4])
L =int(sys.argv[5])


# Load dataframe of parameter sets
df= pickle.load( open(modelling_hpc + '/3954/parameter_space_search/parameterfiles/df_circuit%r_variant%s_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), "rb" ) )
start = int(sys.argv[6])
end = int(sys.argv[7])
total_params = int(end-start)
print('loaded')
batch_size = int(total_params/Number_of_Threads) + 1

df = df.iloc[start:end]
batch_indices = list(range(0, len(df), batch_size))

# Create a pool of workers
pool = multiprocessing.Pool(Number_of_Threads)

# Define jobs as different batches of parameter sets
# Run lsa_check function in parallel across different threads
pool_output = []
for start_batch_index in batch_indices:

    print('main' + str(start_batch_index))
    df_batch = df.iloc[start_batch_index:start_batch_index+batch_size]
    # df_batch = df.iloc[start_batch_index:start_batch_index+2]
    pool_output.append(pool.apply_async(numerical_check, args=(start_batch_index, n_param_sets,df_batch,x_gridpoints,t_gridpoints,T,L)))
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
