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
circuit_n=2
variant=0
# parametersets_n = 1000
boundary_coef=1
shape='ca'
mechanism = 'general'
parID=int(sys.argv[1])
# parID=26
save_figure=False
L,T,J,N = [8,120,80,23880]

#load simulation file
filename = 'timeseries_circuit%r_variant%r_boundary%r_%s_%sID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundary_coef, shape,mechanism,parID,L,J,T,N)
simulation = pickle.load( open( modelling_ephemeral + '/3954/numerical_confocal/results/simulation/1M_colony_%s/2D/2D%s.pkl'%(shape,filename), "rb" ) )
# rgb_timeseries = redgreen_contrast_timeseries(simulation)
# show_rgbvideo(rgb_timeseries)
print(np.shape(simulation[0][int(80/2),:,0]))
lenght_list = []
for i in range(T):
    lenght = np.count_nonzero(simulation[0][int(80/2),:,i])
    lenght_list.append(lenght)
lenght_list = np.array(lenght_list)/10
plt.scatter(np.linspace(0,T,T),lenght_list, c='k',s=1)
plt.xlabel('Time (hours)')
plt.ylabel('Colony diameter (mm)')
plt.show()
