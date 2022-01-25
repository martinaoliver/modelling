# %config Completer.use_jedi = False

import pickle
import numpy as np
import matplotlib as mpl
mpl.use('tkagg')
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
circuit_n=8
variant=6
# parametersets_n = 1000
boundary_coef=1
shape='ca'
mechanism = 'nodeAdele'
parID=int(sys.argv[1])
# parID=26
save_figure=False
L,T,J,N = [8,120,80,23880]
# L,T,J,N = [8,300,80,59700]
#load simulation file
filename = 'timeseries_circuit%r_variant%s_boundary%r_%s_%sID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundary_coef, shape,mechanism,parID,L,J,T,N)
simulation = pickle.load( open( modelling_ephemeral + '/3954/numerical_confocal/results/simulation/1M_colony_%s/2D/nodeA_dele/2D%s.pkl'%(shape,filename), "rb" ) )
print(np.shape(simulation))


rgb_timeseries = redgreen_contrast_timeseries(simulation)
print(np.shape(rgb_timeseries))

show_rgbvideo(rgb_timeseries,parID)
