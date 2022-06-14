#############################
#IMPORTS#
#############################
import sys
import os
import matplotlib

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

modulepath = modelling_local + '/3954/modules/new_CN'

sys.path.append(modulepath)


import numpy as np
from plotting_numerical import plot_redgreen_contrast
from tqdm import tqdm
import matplotlib.pyplot as plt
import pickle

# %matplotlib inline

#execution parameters
circuit_n=2
variant='5716gaussian'
parametersets_n = 1000
save_figure = False
tqdm_disable = False #disable tqdm
n_species=6
# open parameter dictionaries
#chose parameter sets to analyse
parID = 54#int(sys.argv[1])

data_path = modelling_home + '/3954/numerical_confocal/results/simulation/square/fullcircuit_5716gaussian/var0.23'
# parID_list = pickle.load( open(data_path + '/parID_list_L5_J50_T150_N15000.pkl', "rb" ) )
parID_list = pickle.load( open(data_path + '/parID_list_L5_J50_T150_N15000.pkl', "rb" ) )
start = int(0)
# stop = int(len(parID_list)-1)
stop = int(10)


parID_list = [int(i) for i in parID_list] #turn string list into integer list

for parID in tqdm(parID_list):
    # print('parID = ' + str(parID))
    mechanism = 'fullcircuit'
    boundary_coef = 1 #1 is open boundary and 0 is closed boundary
    shape = 'square'
    var=0.23
    growth = True

    # L,J,T,N = [int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])]
    L,J,T,N = [5,50,150,15000]
    t_gridpoints = int(N/T)
    x_gridpoints = int(J/L)

    initial_condition = [0.001]*n_species
    filename = '2Dfinal_circuit%r_variant%svar%s_%s_%sID%s_L%r_J%r_T%r_N%r'%(circuit_n,variant,var,shape,mechanism, parID,L,J,T,N)
    final_concentration = pickle.load( open(modelling_home + '/3954/numerical_confocal/results/simulation/square/fullcircuit_5716gaussian/var0.23/%s.pkl'%filename, 'rb' ) )
    final_concentration = np.round(final_concentration,4)
    save_path =modelling_home + '/3954/numerical_confocal/results/simulation/1M_colony_ca/2D/full_circuit_newCN'
    plt.imshow(final_concentration[-1])
    plt.xticks([], labels=" ")    # plt.colorbar()
    plt.yticks([], labels=" ")    # plt.colorbar()
    # plt.show()
    plt.savefig(modelling_home + '/3954/numerical_confocal/results/figures/square/fullcircuit_5716gaussian/var0.23''/%s.jpeg' %filename)
    plt.close()