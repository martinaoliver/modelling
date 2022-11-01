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
###Imports#####
#############
from numerical.adi_ca_function_openclosed_nodilution import adi_ca_openclosed_nodilution
from numerical.adi_ca_function_openclosed_nodilution_preMask import adi_ca_openclosed_nodilution_preMask
from numerical.plotting_numerical import *
from numerical.adi_square_function import adi

import numpy as np
import pickle
import matplotlib.pyplot as plt
# %matplotlib inline
#############
###execution parameters#####
#############
mechanism = 'fullcircuit'
shape = 'ca'
# parID = int(sys.argv[1])
parID = 1

circuit_n=12
variant=0

# folder = 'fullcircuit/1M_turingI'#'fullcircuit/1M'#'fullcircuit/1M_turingI'
n_species = 9

parametersets_n = 10 #1000000
save_figure = False
tqdm_disable = False #disable tqdm
# boundarycoeff = float(sys.argv[6])
boundarycoeff = 1.5
seed=1;p_division=1;boundaryCoeff=1.5

# open parameter dictionaries
lsa_df= pickle.load( open(modellingpath + '/3954/paper/input/parameterfiles/df_circuit%s_variant%s_%rparametersets.pkl'%(circuit_n,variant,parametersets_n), "rb" ) )
par_dict = lsa_df.loc[parID].to_dict()
# par_dict = general_df.iloc[parID].to_dict()

# d_A = par_dict['dU']
# d_B = par_dict['dV']

D = np.zeros(n_species)
# D[0]=d_A
# D[1]=d_B
D[:2] = [par_dict['DU'],par_dict['DV'] ]
print(par_dict)
# par_dict['mulva'] = par_dict['mulva'] + np.log(2)*p_division


#solver parameters
# L_x=int(sys.argv[2]); x_gridpoints = int(sys.argv[3]); L=L_x; J = L*x_gridpoints;  L_y=L_x; I=J
# T =int(sys.argv[4]); t_gridpoints = int(sys.argv[5]) ; N = T*t_gridpoints


L=5; dx =0.2; J = int(L/dx)
T =10; dt = 0.1; N = int(T/dt)

# suggested_tgridpoints = x_gridpoints**2


filename = 'circuit%r_variant%s_bc%s_%s_%sID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundarycoeff, shape,mechanism,parID,L,J,T,N)
U_record,U_final = adi_ca_openclosed_nodilution_preMask(par_dict,L_x,L_y,J,I,T,N, circuit_n,n_species,D, seed=seed, p_division=p_division, tqdm_disable=tqdm_disable, growth='Fast', boundarycoeff=boundarycoeff)#,p_division=p_division,seed=seed)
# def adi(par_dict,L_x,L_y,J,I,T,N, circuit_n, n_species,D,tqdm_disable=False,stochasticity=0, steadystates=0):

model_args = [par_dict, circuit_n, n_species, D]
solver_args=[L,dx,J,T,dt,N]
system_args = [p_division, seed, boundaryCoeff]
U_record,U_final =  adi_ca_openclosed_nodilution(model_args, solver_args,system_args)

# par_dict,L,dx,J,T,dt,N, circuit_n, n_species,D,tqdm_disable=False, p_division=0.5,stochasticity=0, seed=1,growth='Slow', boundarycoeff=1.5)

plt.imshow(U_final[-1])
plt.show()
plt.imshow(U_final[0])
plt.show()

































