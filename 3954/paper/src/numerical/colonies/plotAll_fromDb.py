#%%
#############
###paths#####
#############
import sys
import os

pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
modellingephemeral = '/rds/general/ephemeral/user/mo2016/ephemeral/Documents/modelling'
sys.path.append(modellingpath + '/lib')
#############

# from numerical.plotAllFunctions import plotAllFunctionColonies_fromDb

import pickle
import numpy as np
import psycopg2
import matplotlib.pyplot as plt
from numerical.plotting_numerical import *

from tqdm import tqdm
#############################

# # slow
# L=20; dx =0.1; J = int(L/dx)
# T =100; dt = 0.02; N = int(T/dt)
# boundaryCoeff = 1
# divisionTimeHours=0.5
# p_division=0.38;seed=1


# # medium
# L=20; dx =0.1; J = int(L/dx)
# T =50; dt = 0.02; N = int(T/dt)
# boundaryCoeff = 1
# divisionTimeHours=0.5
# p_division=1;seed=1



# fast
L=20; dx =0.1; J = int(L/dx)
T =25; dt = 0.02; N = int(T/dt)
boundaryCoeff = 1
divisionTimeHours=0.2
p_division=0.7;seed=1
x_gridpoints=int(1/dx)


shape='ca'
simulation_param_dict = {'L':L, 'dx':dx, 'J':J, 'T':T, 'dt':dt, 'N':N, 
            'boundaryCoeff':boundaryCoeff, 
            'shape':'ca', 'p_division': p_division, 'seed':seed}




credentials=f"postgresql://moliver:moliver@ld-rendres07.bc.ic.ac.uk/moliver"



ssID=0
circuit_n='14' #circuit_n='circuit14'
variant='2nd' #variant='fitted7'
balance='Balanced'
Kce=100
n_samples = 1000000 #n_samples = 13700000
folder = 'circuit14variant2ndBalancedKce100'
model_param_dict = {'circuit_n':circuit_n,'variant':variant, 'n_samples':n_samples, 'balance':balance}


#%%

def find_db_id_from_dict(table_name, param_dict):
    with psycopg2.connect(credentials) as conn:
        with conn.cursor() as cursor:
            # query = f'SELECT {0}_id FROM {table_name} WHERE'
            query = "SELECT {0}_id FROM {0} WHERE ".format(table_name)
            conditions = []
            values = []
            for key, value in param_dict.items():
                conditions.append('"{0}" = %s'.format(key))
                values.append(value)
            cursor.execute(query + ' AND '.join(conditions), values)
            return cursor.fetchall()
        
def getUFinal1D_fromDb(model_param_dict, sim_param_dict):
    template = """
    
    with 
      params_ss1 as (select p."model_param_id" from model_param p, analytical_output a where p.model_param_id = a.model_param_id
                      AND a.ss_n = 1
                      AND {0}),
    sim_params as (select "simulation_param_id" from simulation_param s where {1}  )
        
    select o."U_final_1D", o."model_param_id"
        from simulation_output o,
            params_ss1,
            sim_params,  model_param p
        where o.model_param_id = params_ss1.model_param_id
        and o.simulation_param_id = sim_params.simulation_param_id
        limit 51
        ;
        """
    param_conditions = " AND ".join([ f'p."{k}" = %s' for k in model_param_dict.keys()])
    sim_conditions = " AND ".join([ f's."{k}" = %s' for k in sim_param_dict.keys() ])
    params = list(model_param_dict.values()) + list(sim_param_dict.values())
    query = template.format(param_conditions, sim_conditions)
    # with psycopg2.connect(credentials) as conn:
    #     with conn.cursor() as cursor:
    # query = f'SELECT {0}_id FROM {table_name} WHERE'
    conn=psycopg2.connect(credentials)
    cursor = conn.cursor()
    cursor.execute(query, params)

    return cursor
            



cursor = getUFinal1D_fromDb(model_param_dict, simulation_param_dict)
# U_final = np.array(cursor.fetchmany(2)[0][0], dtype=float)


# plot_redgreen_contrast(U_final, 10,filename=None, path=None, parID=0, scale_factor=10, save_figure=False)

# def plotAllFunctionColonies_fromDb(parID_list, circuit_n, shape, filename, L,x_gridpoints,start=0, stop=10,folder=None, modellingpath=modellingpath, saveFig=True,dpi=2000, tqdm_disable=True, print_parID=False):

def plotAllFunctionColonies_fromDb(cursor,L, x_gridpoints, len_dataset=2):
    cursor_output=cursor.fetchmany(len_dataset)
    




    n_col = int(np.sqrt(len_dataset))
    n_row = int(np.floor(len_dataset/n_col)+1)    # number of rows in the figure of the cluster


    fig = plt.figure(figsize=(n_col/10+12, n_row/10+12))

    for count,[U_final,model_param_id] in tqdm(enumerate(cursor_output)):
        U_final = np.array(U_final, dtype=float)
        print(count,model_param_id)
        print(np.sum(U_final))
    #     if print_parID == True:
    #         print(parID)

        ax=plt.subplot(n_row,n_col, count+1)
        rgb = plot_redgreen_contrast(U_final,L,parID=model_param_id,scale_factor=x_gridpoints,save_figure='LargeImage')

        ax.set(yticks=[],xticks=[],yticklabels=[],xticklabels=[])
        # ax.imshow(rgb.astype('uint8'), origin= 'lower', norm=colors.LogNorm())
        ax.imshow(rgb.astype('uint8'), origin= 'lower')
        ax.set_ylabel(model_param_id,size= 10,c='y', labelpad=0.35)
    plt.show()


    
    # if saveFig==False:
    #     plt.show()
    # if saveFig==True:
    #     if stop==len_fullDataset and start==0:
    #         plt.savefig(modellingpath + '/3954/paper/out/numerical/colonies/largeFigs/%s/largeFig_%s.png'%(folder,filename('x')),dpi=dpi)
    #         plt.savefig(modellingpath + '/3954/paper/out/numerical/colonies/largeFigs/%s/largeFig_%s.pdf'%(folder,filename('x')),dpi=dpi)
        
    #     else:
    #         plt.savefig(modellingpath + '/3954/paper/out/numerical/colonies/largeFigs/%s/largeFig_%s_%s-%s.png'%(folder,filename('x'),start,stop),dpi=dpi)
    #         plt.savefig(modellingpath + '/3954/paper/out/numerical/colonies/largeFigs/%s/largeFig_%s_%s-%s.pdf'%(folder,filename('x'),start,stop),dpi=dpi)
    #         print('not full')
    #         plt.close()
    # # plt.savefig(modelling_home + '/3954/numerical_confocal/results/figures/%s/large_images/%s_%s-%s.png'%(shape,filename,start,stop), dpi=2000)
    # # plt.savefig(modelling_home + '/3954/numerical_confocal/results/figures/ca/large_images/%s_%s.png'%(filename,details), dpi=2000)
    # # plt.savefig(modelling_home + '/3954/numerical_confocal/results/figures/ca/large_images/%s_%s.png'%(filename,details), dpi=2000)
    # x='x'
    # print(f'Done plotting {filename(x)}')


plotAllFunctionColonies_fromDb(cursor, L, x_gridpoints, len_dataset=4)










# %%
