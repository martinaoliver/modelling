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
# %matplotlib inline
shape = 'ca'
# parID = int(sys.argv[1])
parID = 1

circuit_n=2
variant=0

# folder = 'fullcircuit/1M_turingI'#'fullcircuit/1M'#'fullcircuit/1M_turingI'
n_species = 6

nsamples =  10
save_figure = False
tqdm_disable = False #disable tqdm
# boundarycoeff = float(sys.argv[6])
boundarycoeff = 1.7
p_division=0.5;seed=1

# open parameter dictionaries
df= pickle.load( open(modellingpath + '/3954/paper/input/parameterfiles/df_circuit%s_variant%s_%rparametersets.pkl'%(circuit_n,variant,nsamples), "rb" ) )
# instabilities_df= pickle.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/instabilities_dataframes/instability_df_circuit%s_variant%r_%rparametersets.pkl'%(circuit_n,variant,nsamples), "rb" ) )


#solver parameters
L=8; dx =0.05; J = int(L/dx)
T =125; dt = 0.05; N = int(T/dt)




cell_matrix_record = pickle.load( open(modellingpath + "/3954/paper/out/numerical/masks/caMask_seed%s_pdivision%s_L%s_J%s_T%s_N%s.pkl"%(seed,p_division,L,J,T,N), "rb" ) )
daughterToMotherDictList = pickle.load( open(modellingpath + "/3954/paper/out/numerical/masks/caMemory_seed%s_pdivision%s_L%s_J%s_T%s_N%s.pkl"%(seed,p_division,L,J,T,N), "rb" ) )
T =1; dt = 0.05; N = int(T/dt)
filename= lambda parID: 'circuit%r_variant%s_bc%s_%s_ID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundarycoeff, shape,parID,L,J,T,N)


for parID in range(10):
    par_dict = df.iloc[parID].to_dict()

    D = np.zeros(n_species)
    D[:2] = [par_dict['DA'],par_dict['DB'] ]
    
    print(par_dict)

    U_record,U_final =  adi_ca_openclosed_nodilution_preMask(par_dict,L,dx,J,T,dt,N, circuit_n, n_species,D,cell_matrix_record, daughterToMotherDictList,tqdm_disable=False, p_division=0.5,stochasticity=0, seed=1,growth='Slow', boundarycoeff=boundarycoeff)


    plt.imshow(U_final[0])
    plt.show()
    pickle.dump(U_final, open(modellingpath + '/3954/paper/out/numerical/colonies/simulation/2Dfinal_%s.pkl'%filename(parID), "wb" ) )




    savefigpath = modellingpath + '/3954/paper/out/numerical/colonies/figures/'
    plot_redgreen_contrast(U_final,mm=L)





























