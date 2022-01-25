# %config Completer.use_jedi = False

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

while True:
    try:
        from numerical_solvers_variableboundary import *
        break
    except ImportError:
        modelling_hpc = '/rds/general/user/mo2016/home/Documents/modelling'
        modulepath = modelling_hpc + '/6eq/modules'
        sys.path.append(modulepath)
from numerical_solvers_variableboundary import *
from linear_stability_analysis import detailed_turing_analysis_dict
from randomfunctions import plot_highest_dispersion, plot_all_dispersion
# #system parameters
circuit_n=2
variant=0
parametersets_n = 1000
boundary_coef=1
shape='growing_colony'
mechanism = 'general'
parID=int(sys.argv[1])
save_figure=True
L,T,J,N = [8,100,80,19900]

#load simulation file
parametersets_n=1000
# general_df = pickle.load(open(modelling_hpc + '/6eq/parameter_space_search/results/output_dataframes/lsa_df_circuit%r_variant%r_%rparametersets.pkl'%(circuit_n,variant,parametersets_n), "rb"))
general_df = pickle.load(open(modelling_hpc + '/6eq/parameter_space_search/parameterfiles/df_circuit%r_variant%r_%rparametersets.pkl'%(circuit_n,variant,parametersets_n), "rb"))
par_dict = general_df.loc[parID].to_dict()

steadystatelist, number_steadystates, maxeig_list, pattern_class_list, eigenvalues_list = detailed_turing_analysis_dict(par_dict,2,6)
# print(np.shape(eigenvalues_list[0]))
# # print(eigenvalues_list[0])
plot_highest_dispersion(eigenvalues_list[0])
# plot_all_dispersion(eigenvalues_list[0],crop=10)
