#%%
#####
###paths#####
#############
import sys
import os

pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
modellingephemeral = '/rds/general/ephemeral/user/mo2016/ephemeral/Documents/modelling'
sys.path.append(modellingpath + '/lib')
#############

from numerical.cn_plot import plot1D, surfpattern
# from numerical.fourierAnalysisFunctions import psEntropyFunction, plotFourier
from numerical.generalFunctions import round_it
from analytical.linear_stability_analysis import detailed_turing_analysis_dict
from randomfunctions import plot_all_dispersion, plot_highest_dispersion
from pattern_classification.pattern_1D_nogrowth_classification_noRegularity import patternClassification_nogrowth_noRegularity, countPeaks
from pattern_classification.pattern_analysis_functions import *
from database.databaseFunctions import *
import pickle
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm








#%%
# L=25; dx =0.05; J = int(L/dx)
# T =2000; dt = 0.005; N = int(T/dt)

#solver parameters
L=100; dx =0.2; J = int(L/dx)
T =18000; dt = 0.05; N = int(T/dt)


# L=50; dx =0.1; J = int(L/dx)
# T =5000; dt = 0.02; N = int(T/dt)

x_grid = np.array([j*dx for j in range(J)])

boundaryCoeff=1;rate=L/T
suggesteddt = float(dx*dx*2)
mechanism = 'nogrowth'
simulation_param_dict = {'L':L, 'dx':dx, 'J':J, 'T':T, 'dt':dt, 'N':N, 
            'boundaryCoeff':boundaryCoeff, 
            'mechanism':mechanism, 'growth_rate': rate}

parID = 'x'
circuit_n='turinghill'
# variant= int(sys.argv[1])
# n_samples=int(sys.argv[2])

# variant= 0

folder = lambda variant: f'{circuit_n}_variant{variant}'
filename= lambda parID,variant: 'circuit%s_variant%s_bc%s_%s_rate%s_ID%s_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundaryCoeff, mechanism,rate,parID,L,J,T,N)
# %%


data_path = lambda variant: modellingephemeral + f'/growth/out/numerical/{mechanism}/simulation/{folder(variant)}'




query = f'''select mp."parID", so."ssID", mp."variant", mp."n_samples"  from simulation_output so
inner join model_param mp on so.model_param_id = mp.model_param_id
inner join analytical_output ao on (ao.model_param_id,ao."ssID") = (so.model_param_id, so."ssID")

-- where ao.system_class not in ('turing I', 'turing II', 'turing I hopf', 'turing I oscillatory', 'turing II hopf', 'turing semi-hopf')
-- where ao.system_class in ('hopf')
and so.simulation_param_uuid='e9956efa-c4f1-439c-8146-5d009bd107d8'
'''
parIDssID = general_query(query)


 #%%
for parID,ssID,variant,n_samples in tqdm(parIDssID[0]):
    #model param dict
    model_param_dict = {'parID':parID, 'circuit_n':circuit_n,'variant':variant, 'n_samples':n_samples}
    parIDdotssID =f'{parID}.{ssID}'
    print(model_param_dict)
    #load simulations
    U_final_1D = pickle.load( open(data_path(variant) + '/2Dfinal_%s.pkl'%(filename(parIDdotssID,variant)), 'rb'))
    # U_final = np.round(U_final,decimals=4)
    U_record_1D = pickle.load( open(data_path(variant) + '/2Drecord_%s.pkl'%(filename(parIDdotssID,variant)), 'rb'))
    # peaks = countPeaks(U_final, showPlot1D=False)
    # U_final_1D = query_simulationOutput_multiple_from_sql(simulation_param_dict,model_param_dict,'U_final_1D', ssID=0)
    # U_record_1D = query_simulationOutput_multiple_from_sql(simulation_param_dict,model_param_dict,'U_record_1D', ssID=0)

#show simulations
    plot=False
    if plot==True:
        # plot1D(U_final_1D, plotPeaks=False)
        # plt.show()
        surfpattern(U_record_1D,L,dx,J,T, 'linear',  morphogen=0, rate=0, savefig=False,filename='',logResults=False,normalize=False)
        plt.show()

    #classify simulations
    pattern_class, converged, flat = patternClassification_nogrowth_noRegularity(U_record_1D)
    print( pattern_class, converged, flat)

    # insert classification into psql 
    insert_patternClassOutput_to_sql(simulation_param_dict,model_param_dict,ssID,pattern_class, 'pattern_class_nogrowth',allow_update=True)

    # numerical_wavelength = find_wavelenght(U_final, x_grid,showplot1D=False)
    # print('wvl', numerical_wavelength)

    # convergence_time = find_convergence(U_record)
    # print('time', convergence_time)

    # insert_wavelength_convergence_to_sql(simulation_param_dict,model_param_dict,ssID, numerical_wavelength, convergence_time)

    
    


    








# %%
