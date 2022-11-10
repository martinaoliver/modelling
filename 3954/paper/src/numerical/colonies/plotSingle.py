#############
###paths#####
#############
import os
import sys

from importlib_metadata import distribution

pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')

import pickle
import time

import matplotlib.pyplot as plt
import numpy as np
#############
###Imports#####
#############
from numerical.plotting_numerical import *

#############
###execution parameters#####
#############
# %matplotlib inline
shape = 'ca'

parID = 383859
parID = 13974
circuit_n=2
variant= 3
n_species=6
boundarycoeff = 1.7
p_division=0.5;seed=1
folder = 'circuit2variant3_turing'

save_figure = False
nsamples=1000000


# open parameter dictionaries
# df= pickle.load( open(modellingpath + '/3954/paper/input/lhs_parameterfiles/df_circuit%s_variant%s_%rparametersets.pkl'%(circuit_n,variant,nsamples), "rb" ) )
df= pickle.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/turing_dataframes/turing_df_circuit%s_variant%r_%rparametersets.pkl'%(circuit_n,variant,nsamples), "rb" ) )
par_dict = df.loc[parID].to_dict()
print(par_dict)

#solver parameters
L=8; dx =0.05; J = int(L/dx)
T =125; dt = 0.05; N = int(T/dt)

T =1; dt = 0.05; N = int(T/dt)




filename= lambda parID: 'circuit%r_variant%s_bc%s_%s_ID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundarycoeff, shape,parID,L,J,T,N)




#load image
U_final = pickle.load( open(modellingpath + '/3954/paper/out/numerical/colonies/simulation/%s/2Dfinal_%s.pkl'%(folder,filename(parID)), 'rb'))
U_record = pickle.load( open(modellingpath + '/3954/paper/out/numerical/colonies/simulation/%s/2Drecord_%s.pkl'%(folder,filename(parID)), 'rb'))
plt.imshow(U_final[-1])
plt.colorbar()
plt.show()
plt.imshow(U_final[-2])
plt.colorbar()
plt.show()

# savefig_path = modellingpath + '/3954/paper/out/numerical/colonies/figures/%s'%(folder)
# rgb = plot_redgreen_contrast(U_final,L,path = savefig_path,filename=filename(parID),parID=parID,scale_factor=int(1/dx),save_figure=True)















