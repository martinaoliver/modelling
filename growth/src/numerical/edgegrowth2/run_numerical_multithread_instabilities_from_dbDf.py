#%%
#############
###paths#####
#############
import sys
import os

pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
# modellingephemeral = '/rds/general/ephemeral/user/mo2016/ephemeral/Documents/modelling'
modellingephemeral = '/rds/general/user/mo2016/ephemeral/Documents/modelling'
sys.path.append(modellingpath + '/lib')
#############

from numerical.cn_edgegrowth2_numba import cn_edgegrowth2 as cn_edgegrowth2_numba
from numerical.cn_nogrowth import cn_nogrowth

from numerical.cn_plot import plot1D, surfpattern
from database.databaseFunctions import insert_simulationOutput_to_sql, general_query


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
# '''
# Set number of threads to 1 if no valid number provided
if len(sys.argv) > 1:
    Number_of_Threads = int(sys.argv[1])
else:
    Number_of_Threads = 1
print('Number of Threads set to ', Number_of_Threads)


#df with only instabilities of 8 and 9
circuit_n='turinghill'
variant = int(sys.argv[2])
n_samples=1000000

# df= pickle.load( open(modellingpath + f'/growth/out/analytical/instability/multiinstability_df_circuit{circuit_n}_variant{variant}_combinedparametersets.pkl', 'rb'))
df= pickle.load( open(modellingpath + f'/growth/out/analytical/lsa_dataframes/multiinstability_lsa_df_circuitturinghill_variant{variant}_combinedparametersets.pkl','rb'))

query = f'''select mp."parID", so."ssID"  from simulation_output so
join model_param mp on mp.model_param_id = so.model_param_id
where simulation_param_uuid='132323a4-3f93-4287-aca9-d18e84848e37'
and mp.variant='{variant}'
and mp.n_samples={n_samples};'''
simulated_parID_ss = general_query(query)

df = df.drop(simulated_parID_ss[0])

#%%





# Specify date today
date = date.today().strftime('%m_%d_%Y')



def numerical_check(df,a):


    # bigger field
    # L=500; dx =1; J = int(L/dx)
    # T =3000; dt = 0.05; N = int(T/dt)
    # boundaryCoeff=2;rate=0.1
    
    #smaller field
    # L=50; dx =1; J = int(L/dx)
    # T =3000; dt = 0.05; N = int(T/dt)
    # boundaryCoeff=2;rate=0.01

    # smaller time and smaller dt 
    if Number_of_Threads == 1:
        test=True
        print('test')
    else:
        test=False


    # #solver parameters
    # L=50; dx =0.1; J = int(L/dx)
    # T =2000; dt = 0.02; N = int(T/dt)
    # rate=L/T
    # suggesteddt = float(dx*dx*2)

    #solver parameters
    L=25; dx =0.05; J = int(L/dx)
    T =2000; dt = 0.005; N = int(T/dt)
    rate=L/T
    suggesteddt = float(dx*dx*2)

    if test == True:
        T =10; dt = 0.1; N = int(T/dt)
        tqdm_disable = False
        rate=L/T

    else:
        tqdm_disable = True
    # df_index = np.unique(df.index.get_level_values(0))
    for parID,ssID in df.index:
        parIDssID = f'{parID}.{ssID}'
        print('parID = ' + str(parIDssID))
        par_dict = df.loc[(parID,ssID)].to_dict()


        folder=  f'{circuit_n}_variant{variant}'

        model_param_dict = {'parID':parID, 'circuit_n':circuit_n,'variant':variant, 'n_samples':n_samples}


        simulateGrowth=True
        simulateNoGrowth=True
        simulateOpenBoundary=True
        try:

                
            if simulateNoGrowth == True:
                mechanism = 'nogrowth'
                boundaryCoeff=1
                U_final_1D,U_record_1D, U0, x_grid, reduced_t_grid= cn_nogrowth(par_dict,L,J,T,N, circuit_n,boundaryCoeff=boundaryCoeff, tqdm_disable=tqdm_disable)
                filename= lambda mechanism, parID: 'circuit%s_variant%s_bc%s_%s_rate%s_ID%s_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundaryCoeff, mechanism,rate,parID,L,J,T,N)
                simulation_param_dict = {'L':L, 'dx':dx, 'J':J, 'T':T, 'dt':dt, 'N':N, 'boundaryCoeff':boundaryCoeff, 'mechanism':mechanism, 'growth_rate': rate}
                query = insert_simulationOutput_to_sql(simulation_param_dict, model_param_dict,U_final_1D,U_record_1D, ssID,dimensions='1D',allow_update=True)
                with open(modellingephemeral + f'/growth/out/numerical/{mechanism}/simulation/{folder}/2Dfinal_{filename(mechanism,parIDssID)}.pkl', 'wb') as f:
                    pickle.dump(U_final_1D, f)
                with open(modellingephemeral + f'/growth/out/numerical/{mechanism}/simulation/{folder}/2Drecord_{filename(mechanism,parIDssID)}.pkl', 'wb') as f:
                    pickle.dump(U_record_1D, f)


            if simulateOpenBoundary == True:
                mechanism = 'openboundary'
                boundaryCoeff=2
                U_final_1D,U_record_1D, U0, x_grid, reduced_t_grid= cn_nogrowth(par_dict,L,J,T,N, circuit_n,boundaryCoeff=boundaryCoeff, tqdm_disable=tqdm_disable)
                filename= lambda mechanism, parID: 'circuit%s_variant%s_bc%s_%s_rate%s_ID%s_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundaryCoeff, mechanism,rate,parID,L,J,T,N)
                simulation_param_dict = {'L':L, 'dx':dx, 'J':J, 'T':T, 'dt':dt, 'N':N, 'boundaryCoeff':boundaryCoeff, 'mechanism':mechanism, 'growth_rate': rate}
                query = insert_simulationOutput_to_sql(simulation_param_dict, model_param_dict,U_final_1D,U_record_1D,ssID, dimensions='1D',allow_update=True)
                with open(modellingephemeral + f'/growth/out/numerical/{mechanism}/simulation/{folder}/2Dfinal_{filename(mechanism,parIDssID)}.pkl', 'wb') as f:
                    pickle.dump(U_final_1D, f)
                with open(modellingephemeral + f'/growth/out/numerical/{mechanism}/simulation/{folder}/2Drecord_{filename(mechanism,parIDssID)}.pkl', 'wb') as f:
                    pickle.dump(U_record_1D, f)



            if simulateGrowth == True:
                mechanism = 'edgegrowth2'
                boundaryCoeff=2

                U_final_1D,U_record_1D, U0, x_grid, reduced_t_grid, cellMatrix= cn_edgegrowth2_numba(par_dict,L,J,T,N, circuit_n, rate=rate, boundaryCoeff=boundaryCoeff, tqdm_disable=tqdm_disable)
                simulation_param_dict = {'L':L, 'dx':dx, 'J':J, 'T':T, 'dt':dt, 'N':N, 'boundaryCoeff':boundaryCoeff, 'mechanism':mechanism, 'growth_rate': rate}
                filename= lambda mechanism, parID: 'circuit%s_variant%s_bc%s_%s_rate%s_ID%s_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundaryCoeff, mechanism,rate,parID,L,J,T,N)
                query = insert_simulationOutput_to_sql(simulation_param_dict, model_param_dict,U_final_1D,U_record_1D,ssID, dimensions='1D',allow_update=True)
                with open(modellingephemeral + f'/growth/out/numerical/{mechanism}/simulation/{folder}/2Dfinal_{filename(mechanism,parIDssID)}.pkl', 'wb') as f:
                    pickle.dump(U_final_1D, f)
                with open(modellingephemeral + f'/growth/out/numerical/{mechanism}/simulation/{folder}/2Drecord_{filename(mechanism,parIDssID)}.pkl', 'wb') as f:
                    pickle.dump(U_record_1D, f)


        except ValueError:
            print('!!!!!!!!!!!!!')
            print('ValueError --> unstable solution')
            print('!!!!!!!!!!!!!')
            print()

            pass



start_time = time.perf_counter()


df.index.names = ['parID','ssID']
total_params=len(df)
print(df)
print(f'len(df) = {len(df)}')
print('loadedd')

print(f'total params {total_params}')
batch_size = int(total_params/Number_of_Threads) + 1
print(df.head())
batch_indices = list(range(0, total_params, batch_size))
# Create a pool of workers
pool = multiprocessing.Pool(Number_of_Threads)

# Run lsa_check function in parallel across different threads
pool_output = []
for start_batch_index in batch_indices:
    print('main' + str(start_batch_index))
    df_batch = df.iloc[start_batch_index:start_batch_index+batch_size]
    pool_output.append(pool.apply_async(numerical_check, args=(df_batch,'a')))

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

