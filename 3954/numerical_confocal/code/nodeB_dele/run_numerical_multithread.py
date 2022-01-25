'''
====================================================
    Imports
====================================================
'''

modulepath = '/Users/mo2016/Documents/modelling/3954/modules'
import sys
sys.path.append(modulepath)
modelling_hpc = '/Volumes/mo2016/home/Documents/modelling'
modelling_ephemeral = '/Volumes/mo2016/ephemeral/Documents/modelling'
while True:
    try:
        from numerical_solvers_variableboundary import *
        break
    except ImportError:
        modelling_hpc = '/rds/general/user/mo2016/home/Documents/modelling'
        modulepath = modelling_hpc + '/3954/modules'
        modelling_ephemeral = '/rds/general/user/mo2016/ephemeral/Documents/modelling'
        sys.path.append(modulepath)

from numerical_solvers_variableboundary import *
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

# Specify name of circuit and variant investigated
circuit_n=9
variant= 0
n_species = 5
# Specifiy number of parameter sets in parameterset file to be loaded
n_param_sets = 1000000


# Specify date today
date = date.today().strftime('%m_%d_%Y')

# Specify size of batches in which to complete computations
# Does not need to be a factor of number of parameter sets
# batch_size = 1000

# Get starting parameter set indices for batches


# Define work to be done per batch of parameter sets


def numerical_check(start_batch_index,n_param_sets,df, dimension,x_gridpoints,T,L,circuit_n=9, variant=0, n_species=5,p_division=0.41,seed=1):
    save_figure = True
    tqdm_disable = True #disable tqdm

    df_index = np.unique(df.index.get_level_values(0))
    for parID in df_index:
        print('parID = ' + str(parID))
        mechanism = 'nodeBdele'
        boundary_coef = 1 #1 is open boundary and 0 is closed boundary
        shape = 'ca'
        growth = True
        # x_gridpoints = int(sys.argv[3])
        # T =int(sys.argv[4])
        par_dict = df.loc[parID].to_dict()
        par_dict['ba']=1
        par_dict['bc']=1
        par_dict['bd']=1
        par_dict['be']=1
        par_dict['bf']=1
        J = L *x_gridpoints  # number of equally spaced gridpoints in space domain (larger J means more spatial precision(tends towards continuum solution) )
        t_gridpoints = t_gridpoints_stability(L, J, T)  # number of equally spaced gridpoints in domain (larger N means more temporal precision (tends towards continuum solution) )

        N = T * t_gridpoints
        initial_condition = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001]

        filename = 'circuit%r_variant%r_boundary%r_%s_%sID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundary_coef, shape,mechanism,parID,L,J,T,N)

        if dimension == '1D':
            try:
                records, final_concentration, grids = crank_nicolson_ca(par_dict, initial_condition, L, J, T, N, circuit_n,n_species=n_species,boundary_coef=boundary_coef,growth = growth,tqdm_disable=tqdm_disable)
                plot_1D_final_concentration(final_concentration,grids,mechanism,shape,filename,n_species, parID,save_figure=save_figure,path=modelling_home)
                # plot_redgreen_contrast(final_concentration,L,mechanism,shape,filename,modelling_ephemeral,dimension=dimension,save_figure=False)
                if save_figure ==True:
                    pickle.dump(final_concentration, open(modelling_ephemeral + '/3954/numerical_confocal/results/simulation/1M_colony_ca/1D/1Dfinal_%s.pkl'%filename, 'wb'))
                    pickle.dump(records,open(modelling_ephemeral + '/3954/numerical_confocal/results/simulation/1M_colony_ca/1D/1Dtimeseries_%s.pkl'%filename, 'wb'))
                    timeseries_1D_rgb(records, L, J,T,filename,modelling_home,save_figure=save_figure)
            except ValueError:
                print('!!!!!!!!!!!!!')
                print('ValueError --> unstable solution')
                print('!!!!!!!!!!!!!')
                print()

                pass

            # plot_1D_final_concentration(final_concentration, grids,mechanism,shape,filename,parID,save_figure=save_figure,path=modelling_home)
            # plot_redgreen_contrast(final_concentration,grids,mechanism,shape,filename,parID,modelling_home,dimension=dimension)
            # plt.show()
        if dimension == '2D':
         # Define 2D numerical parameters
            L_x = L
            L_y = L
            I = J
            try:

                # Run 2D simulation

                records,final_concentration,grids = adi_ca(par_dict,initial_condition,L_x,L_y,J,I,T,N, circuit_n, boundary_coef=boundary_coef,tqdm_disable=tqdm_disable,n_species=n_species,p_division=p_division,seed=seed)
                # plot_2D_final_concentration(final_concentration,grids,filename,n_species=n_species)
                nodeB_dele_path = modelling_ephemeral + '/3954/numerical_confocal/results/figures/1M_colony_ca/2D/nodeB_dele'
                plot_redgreen_contrast(final_concentration,L,mechanism,shape,filename,nodeB_dele_path ,parID=parID,dimension=dimension,scale_factor=x_gridpoints,save_figure=save_figure)
                # rgb_timeseries = redgreen_contrast_timeseries(records)
                # show_rgbvideo(rgb_timeseries)
                if save_figure ==True:
                    pickle.dump(final_concentration, open(modelling_ephemeral + '/3954/numerical_confocal/results/simulation/1M_colony_ca/2D/nodeB_dele/2Dfinal_%s.pkl'%filename, 'wb'))
                    pickle.dump(records,open(modelling_ephemeral + '/3954/numerical_confocal/results/simulation/1M_colony_ca/2D/nodeB_dele/2Dtimeseries_%s.pkl'%filename, 'wb'))

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
dimension = str(sys.argv[2])
x_gridpoints = int(sys.argv[3])
T =int(sys.argv[4])
L =int(sys.argv[5])

# Load dataframe of parameter sets
df= pickle.load( open(modelling_hpc + '/3954/parameter_space_search/parameterfiles/df_circuit%r_variant%r_%rparametersets.pkl'%(2,0,n_param_sets), "rb" ) )
start = int(sys.argv[6])
end = int(sys.argv[7])
total_params = int(end-start)

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
    pool_output.append(pool.apply_async(numerical_check, args=(start_batch_index, n_param_sets,df_batch,dimension,x_gridpoints,T,L)))
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
