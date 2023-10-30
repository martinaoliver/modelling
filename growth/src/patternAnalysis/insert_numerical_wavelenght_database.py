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
from pattern_classification.pattern_analysis_functions import find_wavelenght
from pattern_classification.pattern_1D_nogrowth_classification import patternClassification
#%%




#solver parameters
L=25; dx =0.05; J = int(L/dx)
T =2000; dt = 0.005; N = int(T/dt)
x_grid = np.array([j*dx for j in range(J)])

boundaryCoeff=1;rate=L/T
suggesteddt = float(dx*dx*2)
mechanism = 'nogrowth'
simulation_param_dict = {'L':L, 'dx':dx, 'J':J, 'T':T, 'dt':dt, 'N':N, 
            'boundaryCoeff':boundaryCoeff, 
            'mechanism':mechanism, 'growth_rate': rate}
print(simulation_param_dict.keys(), simulation_param_dict.values())



circuit_n='turinghill'
variant= '11'
n_samples=1000000

from scipy.signal import find_peaks







def find_wavelenght(U,x_grid,showplot1D=True):
    peaks = [0, 0]
    prominence=0.05

    peaks[0], _ = find_peaks(U[0], prominence=prominence)
    peaks[1], _ = find_peaks(U[1], prominence=prominence)

    # Calculate the wavelength
    wavelength_x = np.mean(np.diff(x_grid[peaks[0]]))
    wavelength_y = np.mean(np.diff(x_grid[peaks[1]]))
    avg_wavelength = np.mean([wavelength_x, wavelength_y])
    # Plot the 1D signal and peaks
    if showplot1D:
        plot1D(U, peaks=peaks)

    return  avg_wavelength


query = f'''select mp."parID", so."ssID"  from simulation_output so
inner join model_param mp on mp.model_param_id = so.model_param_id
inner join analytical_output ao on (ao.model_param_id,ao."ssID") = (so.model_param_id, so."ssID")
where ao.system_class = 'turing I oscillatory'
and simulation_param_uuid='6952d306-f619-4af1-963c-aa28acb132df'
and mp.variant='{variant}'
and mp.n_samples={n_samples}
and ss_n=1;'''
simulated_parID_ss = general_query(query)[0]
for parID,ssID in simulated_parID_ss[:20]:
    print(parID,ssID)

    model_param_dict = {'parID':parID, 'circuit_n':circuit_n,'variant':variant, 'n_samples':n_samples}
    U_final_1D = query_simulationOutput_single_from_sql(simulation_param_dict,model_param_dict,'U_final_1D', ssID=ssID)
    # U_final_1D = query_simulationOutput_multiple_from_sql(simulation_param_dict,model_param_dict,'U_record_1D', ssID=0,fetch=2)
    plot1D(U_final_1D, L=L)
    plt.show()
    # todo: measure wavelenght
    # todo: add wavelenght to pattern class
    numerical_wavelenght = find_wavelenght(U_final_1D, x_grid,showplot1D=False)
    print('wvl',numerical_wavelenght)




    U_record_1D = query_simulationOutput_single_from_sql(simulation_param_dict,model_param_dict,'U_record_1D', ssID=ssID)


    # %%
    def find_convergence(U_record):
        #check if converged
        relRangeConverged=[0,0]
        for time in np.arange(200, 2,-1):
            for count,Ux_record in enumerate(U_record):
                relRangeConverged[count] = [(np.amax(x) - np.amin(x))/(np.amax(x)+1e-8) for x in np.transpose(Ux_record[time:time+3])]
            # if np.amax(relRangeConverged[0])>0.001 or np.amax(relRangeConverged[1])>0.001:
            if np.amax(relRangeConverged[0])>0.05 or np.amax(relRangeConverged[1])>0.05:
                converged=False
                return time*10
            else:
                converged=True


    surfpattern(U_record_1D,L,dx,J,T)
    plt.show()
    time = find_convergence(U_record_1D)
    print('time', time)

    print('----------------')
# %%
