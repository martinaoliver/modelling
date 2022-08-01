#############################
#IMPORTS#
#############################
import sys
import os

from numpy import searchsorted
pwd = os.getcwd()
root = pwd.rpartition("mo2016")[0] + pwd.rpartition("mo2016")[1] #/Volumes/mo2016/ or '/Users/mo2016/' or '/rds/general/mo2016/'
if root == '/Users/mo2016':
    print(root)
    modelling_ephemeral = '/Volumes/mo2016/ephemeral/Documents/modelling'
    modelling_home = '/Volumes/mo2016/home/Documents/modelling'
    modelling_local = root + '/Documents/modelling'


if root == '/Volumes/mo2016' or root=='/rds/general/user/mo2016': #'/rds/general' or root=='/Volumes':
        modelling_ephemeral = root + '/ephemeral/Documents/modelling'
        modelling_home = root  + '/home/Documents/modelling'
        modelling_local = modelling_home

if root == '/Users/mo2016' or  root == '/Volumes/mo2016':
    import matplotlib as mpl
    mpl.use('tkagg')

modulepath = modelling_local + '\lib'

sys.path.append(modulepath)



from numerical.plotting_numerical import *

import pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable

import pandas as pd
import numpy as np
import matplotlib.animation as animation
from tqdm import tqdm
# %matplotlib inline

#############################
#Opening list with parID's
# file = open(modelling_ephemeral + '/3954/numerical_confocal/results/simulation/1M_colony_ca/2D/parID_list_8x10T120.txt')
folder = 'fullcircuit/1M_turingI'#'fullcircuit/1M'#'fullcircuit/1M_turingI'
circuit_n=2
variant=0#9
shape = 'caNodilution'
mechanism = 'fullcircuit'
boundarycoeff=1.5
L=10; x_gridpoints =15; J = L*x_gridpoints
T =120; t_gridpoints = 10; N = T*t_gridpoints

# L=2; x_gridpoints =50; J = L*x_gridpoints
# T =24; t_gridpoints = 100; N = T*t_gridpoints

# L=2; x_gridpoints =50; J = L*x_gridpoints
# T =120; t_gridpoints = 100; N = T*t_gridpoints
details = '1MturingI'

seed=1;p_division=0.5#0.147
## var=0.23
# data_path = modelling_ephemeral + '/3954/numerical_confocal/results/simulation/ca/2D/full_circuit/%s'%folder
data_path = modelling_home + '/3954/numerical_confocal/results/simulation/ca/%s'%folder
parID=198492
parID=383779

filename = 'circuit%r_variant%s_bc%s_%s_%sID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundarycoeff, shape,mechanism,parID,L,J,T,N)
final_concentration = pickle.load( open(data_path + '/2Dfinal_%s.pkl'%filename, 'rb' ) )
# plot_redgreen_contrast(final_concentration,L,mechanism,shape,filename,'',parID=parID,scale_factor=x_gridpoints,save_figure=False)

rgb = plot_redgreen_contrast_nonorm1(final_concentration,L,mechanism,shape,filename,modelling_ephemeral,parID=parID,dimension='2D',scale_factor=x_gridpoints,save_figure=False)
# plt.imshow(rgb.astype('uint8', casting='safe'))
plt.imshow(rgb.astype('uint8'))
plt.show()
# rgb_timeseries = redgreen_contrast_timeseries(U_record)
# show_rgbvideo(rgb_timeseries,parID)
# plt.plot(rgb[int(150/2)])
# plt.show()
# # plt.imshow(final_concentration[-2][150/2])
# # plt.show()
# plt.plot(rgb[int(150/2)])
# plt.show()