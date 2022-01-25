%config Completer.use_jedi = False

import pickle
import numpy as np
import sys
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib.animation as animation
modelling_hpc = '/Volumes/mo2016/home/Documents/modelling'

modulepath = '/Users/mo2016/Documents/modelling/6eq/modules'
import sys
sys.path.append(modulepath)

from class_circuit_eq import circuit2_eq
from scipy.integrate import solve_ivp,odeint
from lhs import *

#LOAD FILES
circuit_n=2
variant=0
n_species=6
parametersets_n = 1000
boundary_coef=1
shape='growing_colony'
mechanism = 'general'
parID=504
save_figure=True
L,T,J,N = [8,100,80,19900]

#load simulation file
parametersets_n=1000
# general_df = pickle.load(open(modelling_hpc + '/6eq/parameter_space_search/results/output_dataframes/lsa_df_circuit%r_variant%r_%rparametersets.pkl'%(circuit_n,variant,parametersets_n), "rb"))
general_df = pickle.load(open(modelling_hpc + '/6eq/parameter_space_search/parameterfiles/df_circuit%r_variant%r_%rparametersets.pkl'%(circuit_n,variant,parametersets_n), "rb"))
results_df = pickle.load(open(modelling_hpc + '/6eq/parameter_space_search/results/output_dataframes/lsa_df_circuit%r_variant%r_%rparametersets.pkl'%(circuit_n,variant,parametersets_n), "rb"))
par_dict = general_df.loc[parID].to_dict()


results_df.loc[parID]

#LOAD EQUATIONS AND FIXED POINTS
ss_list = results_df.loc[parID]['ss_list'].values
ss_list[0]
eq = circuit2_eq(par_dict)
def model(y0,t):
    A,B,C,D,E,F = y0
    return eq.dAdt_f(A,B,C,D,E,F),eq.dBdt_f(A,B,C,D,E,F),eq.dCdt_f(A,B,C,D,E,F),eq.dDdt_f(A,B,C,D,E,F),eq.dEdt_f(A,B,C,D,E,F),eq.dFdt_f(A,B,C,D,E,F)


#ODE SOLVER
y0 = np.zeros(n_species)
t_final=100000
t = np.linspace(0,t_final,t_final)
y = odeint(model, y0,t)

def plot_ODE(t,y,crop_start=0,crop_end=-1):
    line1,line2,line3,line4,line5,line6 = plt.plot(t[crop_start:crop_end],y[crop_start:crop_end])
    plt.scatter(np.full(6,crop_end),ss_list[0],c='k',s=10)
    plt.legend((line1,line2,line3,line4,line5,line6),('A','B','C','D','E','F'))
    plt.yscale('log')
    plt.xlabel('Time')
    plt.ylabel('Molecule concentration (log)')
    plt.show()
plot_ODE(t,y,crop_end=100)

def plot_phaseplane(t,y,crop=0,log=False):
    plt.plot(y[crop:,4],y[crop:,3])
    plt.xlabel('mCherry')
    plt.ylabel('GFP')
    if log==True:
        plt.xscale('log')
        plt.yscale('log')
    plt.show()
plot_phaseplane(t,y,crop=0)


#LHS
t_final=100000
t = np.linspace(0,t_final,t_final)
n_initial_conditions=500
initial_conditions = lhs_initial_conditions(n_initial_conditions,n_species)
ode_matrix = np.ones(shape=(n_initial_conditions,t_final,n_species))

for count,y0 in enumerate(initial_conditions):
    y = odeint(model, y0,t)
    ode_matrix[count]=y

np.shape(ode_matrix[:,-1,0])
plt.scatter(ode_matrix[:,-1,0],ode_matrix[:,-1,1])
plt.scatter(ode_matrix[:,-1,2],ode_matrix[:,-1,3])
plt.scatter(ode_matrix[:,-1,4],ode_matrix[:,-1,5])
plot_ODE(t,ode_matrix[1])
