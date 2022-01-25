# %config Completer.use_jedi = False

import pickle
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib.animation as animation
import sys
import os
pwd = os.getcwd()
root = pwd.split("home", 1)[0]
modelling_home = root + 'home/Documents/modelling'
modelling_ephemeral = root + 'ephemeral/Documents/modelling'
modulepath = modelling_home + '/3954/modules'
sys.path.append(modulepath)

from numerical_solvers_variableboundary import *

#system parameters
circuit_n=4
variant=0
parametersets_n = 1000
boundary_coef=0
shape='growing_colony'
mechanism = 'lost_plasmid'
parID=int(sys.argv[1])
# parID=26
# save_figure=False
L,T,J,N = [8,100,64,12600]

#load simulation file
filename = 'timeseries_circuit%r_variant%r_boundary%r_%s_%sID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundary_coef, shape,mechanism,parID,L,J,T,N)
simulation = pickle.load( open( modelling_ephemeral + '/3954/numerical_confocal/results/simulation/1D_df_%s/%s.pkl'%(mechanism,filename), "rb" ) )

show_rgbvideo(simulation)
# timeseries_1D_rgb(simulation, L, J,T,filename,modelling_hpc,save_figure=save_figure)
