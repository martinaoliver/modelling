#############################
#IMPORTS#
#############################
import sys
import os
pwd = os.getcwd()
root = pwd.rpartition("mo2016")[0] + pwd.rpartition("mo2016")[1] #/Volumes/mo2016/ or '/Users/mo2016/' or '/rds/general/mo2016/'
print(root)
if root == '/Users/mo2016':
    modelling_ephemeral = '/Volumes/mo2016/ephemeral/Documents/modelling'
    modelling_home = '/Volumes/mo2016/home/Documents/modelling'
    modelling_local = root + '/Documents/modelling'
    import matplotlib as mpl
    print("im here")
    mpl.use('tkagg')

if root == '/Volumes/mo2016' or root=='/rds/general/user/mo2016': #'/rds/general' or root=='/Volumes':
        modelling_ephemeral = root + '/ephemeral/Documents/modelling'
        modelling_home = root  + '/home/Documents/modelling'
        modelling_local = modelling_home

modulepath = modelling_local + '/3954/modules/new_CN'

sys.path.append(modulepath)


from adi_ca_function import *
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
circuit_n=2
variant= '5716gaussian'
folder = '5716gaussian'
n_species = 6
# Specifiy number of parameter sets in parameterset file to be loaded
n_param_sets = 30000


# Specify date today
date = date.today().strftime('%m_%d_%Y')
# Specify size of batches in which to complete computations
# Does not need to be a factor of number of parameter sets
# batch_size = 1000

# Get starting parameter set indices for batches


# Define work to be done per batch of parameter sets


def numerical_check(start_batch_index,n_param_sets,df,x_gridpoints, t_gridpoints,T,L,circuit_n=2, variant = variant,folder=folder, n_species=6,p_division=0.5,seed=1):
    save_figure = True
    tqdm_disable = True #disable tqdm

    df_index = np.unique(df.index.get_level_values(0))
    for parID in df_index:
        print('parID = ' + str(parID))
        mechanism = 'fullcircuit'
        shape = 'ca'


        par_dict = df.loc[parID].to_dict()

        d_A = par_dict['d_A']
        d_B = par_dict['d_B']
        D = np.zeros(n_species)
        D[0]=d_A
        D[1]=d_B

        J = L *x_gridpoints  # number of equally spaced gridpoints in space domain (larger J means more spatial precision(tends towards continuum solution) )

        N = T * t_gridpoints
        initial_condition = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001]

        filename = 'circuit%r_variant%s_%s_%sID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant, shape,mechanism,parID,L,J,T,N)


     # Define 2D numerical parameters
        L_x = L
        L_y = L
        I = J
        try:

            # Run 2D simulation
            # U_record,U_final = adi_ca(par_dict,L_x,L_y,J,I,T,N, circuit_n,n_species, D,tqdm_disable=tqdm_disable)#,p_division=p_division,seed=seed)
            U_record,U_final = adi_ca(par_dict,L_x,L_y,J,I,T,N, circuit_n,n_species,D, tqdm_disable=tqdm_disable)#,p_division=p_division,seed=seed)
            savefig_path = modelling_ephemeral + '/3954/numerical_confocal/results/figures/ca/2D/full_circuit/%s'%(folder)
            # plot_2D_final_concentration(U_final,L_x,J,filename,savefig_path,n_species=n_species,save_figure=True)
            plot_redgreen_contrast(U_final,L_x,mechanism,shape,filename,savefig_path,parID=parID,scale_factor=x_gridpoints,save_figure=save_figure)
            # rgb_timeseries = redgreen_contrast_timeseries(records)
            # show_rgbvideo(rgb_timeseries)
            if save_figure ==True:
                pickle.dump(U_final, open(modelling_ephemeral + '/3954/numerical_confocal/results/simulation/ca/2D/full_circuit/%s/2Dfinal_%s.pkl'%(folder,filename), 'wb'))
                pickle.dump(U_record,open(modelling_ephemeral + '/3954/numerical_confocal/results/simulation/ca/2D/full_circuit/%s/2Dtimeseries_%s.pkl'%(folder,filename), 'wb'))

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
L=int(sys.argv[2]); x_gridpoints = int(sys.argv[3])
T =int(sys.argv[4]); t_gridpoints = int(sys.argv[5]) 

# t_gridpoints = int(sys.argv[2])
# x_gridpoints = int(sys.argv[3])
# T =int(sys.argv[4])
# L =int(sys.argv[5])


# Load dataframe of parameter sets
# df= pickle.load( open(modelling_home + '/3954/parameter_space_search/parameterfiles/df_circuit%r_variant%s_%rparametersets.pkl'%(2,variant,n_param_sets), "rb" ) )
df= pickle.load( open(modelling_home + '/3954/parameter_space_search/parameterfiles/5716gaussian/df_circuit%r_variant%s_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), "rb" ) )
start = int(sys.argv[6])
end = int(sys.argv[7])
total_params = int(end-start)
print('loaded')
batch_size = int(total_params/Number_of_Threads) + 1
df = df.iloc[start:end]
batch_indices = list(range(0, len(df), batch_size))

# Create a pool of workers
pool = multiprocessing.Pool(Number_of_Threads)

# Run lsa_check function in parallel across different threads
pool_output = []
for start_batch_index in batch_indices:

    print('main' + str(start_batch_index))
    df_batch = df.iloc[start_batch_index:start_batch_index+batch_size]

    # df_batch = df.iloc[start_batch_index:start_batch_index+2]
    pool_output.append(pool.apply_async(numerical_check, args=(start_batch_index, n_param_sets,df_batch,x_gridpoints, t_gridpoints,T,L)))

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
