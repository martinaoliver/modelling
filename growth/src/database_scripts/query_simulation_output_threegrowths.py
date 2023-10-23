#%%
import sys
import os
pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')
from database.databaseFunctions import *
import pickle
#############
#############
###paths#####
#############

import sys
import os
import pickle
import psycopg2
import matplotlib.pyplot as plt

from numerical.cn_plot import plot1D, surfpattern
from database.databaseFunctions import *

#%%




L=50; dx =0.1; J = int(L/dx)
T =2000; dt = 0.02; N = int(T/dt)
rate=L/T
suggesteddt = float(dx*dx*2)

parID = 8248914
circuit_n='turinghill'
variant= '9'
n_samples=1000000
ssID = 0
folder = f'{circuit_n}_variant{variant}'
model_param_dict = {'parID':parID, 'circuit_n':circuit_n,'variant':variant, 'n_samples':n_samples}
#%%
for mechanism, boundaryCoeff in zip(['nogrowth', 'openboundary', 'edgegrowth2'], [1,2,2]):
    simulation_param_dict = {'L':L, 'dx':dx, 'J':J, 'T':T, 'dt':dt, 'N':N, 'boundaryCoeff':boundaryCoeff, 'mechanism':mechanism, 'growth_rate': rate}

    U_final = query_simulationOutput_single_from_sql(simulation_param_dict,model_param_dict,'U_final_1D', ssID=0)
    plot1D(U_final, savefig=False,filename='')

    U_record = query_simulationOutput_single_from_sql(simulation_param_dict,model_param_dict,'U_record_1D', ssID=0)
    surfpattern(U_record,L,dx,J,T, savefig=False,filename='')
    plt.show()
    plt.close('all')


#%%
lsa_df = pickle.load(open( modellingpath + f'/growth/out/analytical/lsa_dataframes/lsa_df_circuitturinghill_variant9_{n_samples}parametersets.pkl', "rb" ) )
# lsa_df = pickle.load(open( modellingpath + f'/growth/out/analytical/instability/instability_df_circuitturinghill_variant8-9_combinedparametersets.pkl', "rb" ) )
lsa_df.loc[parID]

#%%
#linear stability analysis of parID
from analytical.linear_stability_analysis import big_turing_analysis_df, detailed_turing_analysis_dict
from randomfunctions import plot_all_dispersion

# par_dict = lsa_df.loc[parID].iloc[1].to_dict() #converts a dataframe row into a dictionary outputing a dictionary for a specific parameter set
par_dict = lsa_df.loc[parID].to_dict() #converts a dataframe row into a dictionary outputing a dictionary for a specific parameter set
#Run analysis on 1M parameter sets
out = detailed_turing_analysis_dict(par_dict, circuit_n, 2)
plot_all_dispersion(out[4][0],2, crop=20)



    # %%
