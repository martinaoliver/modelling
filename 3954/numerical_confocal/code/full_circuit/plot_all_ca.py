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
    print('yyyyy')
    import matplotlib as mpl
    mpl.use('tkagg')

modulepath = modelling_local + '/3954/modules/new_CN'

sys.path.append(modulepath)



from plotting_numerical import *

import pickle
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib.animation as animation
from tqdm import tqdm

#############################
#Opening list with parID's
# file = open(modelling_ephemeral + '/3954/numerical_confocal/results/simulation/1M_colony_ca/2D/parID_list_8x10T120.txt')
folder = 'fullcircuit/1M'
circuit_n=2
variant=9
shape='ca'
mechanism = 'fullcircuit'
boundarycoeff=1.5
# L=10; x_gridpoints =15; J = L*x_gridpoints
# T =120; t_gridpoints = 10; N = T*t_gridpoints

L=2; x_gridpoints =50; J = L*x_gridpoints
T =24; t_gridpoints = 100; N = T*t_gridpoints
details = '1M'

seed=1;p_division=0.147
## var=0.23
# print(var)
# data_path = modelling_ephemeral + '/3954/numerical_confocal/results/simulation/ca/2D/full_circuit/%s'%folder
data_path = modelling_home + '/3954/numerical_confocal/results/simulation/ca/%s'%folder
# parID_list = pickle.load( open(data_path + '/parID_list_L5_J50_T150_N15000.pkl', "rb" ) )
# parID_list = pickle.load( open(data_path + '/parID_list_variant5716gaussian_ca_fullcircuit_L10J150T120N1200.pkl', "rb" ) )
parID_list = pickle.load( open(data_path + '/parID_list_circuit%r_variant%s_bc%s_%s_%s_L%r_J%r_T%r_N%r.pkl'%(circuit_n,variant,boundarycoeff,shape,mechanism,L,J,T,N), "rb" ) )
# start = int(sys.argv[1])
start = 0
stop = int(len(parID_list)-1)
# stop = int(sys.argv[2])

parID_list = [int(i) for i in parID_list[start:stop]] #turn string list into integer list
print(len(parID_list))
parID_list.sort() #sort from lower to higher values


num=len(parID_list)
n_col = int(np.sqrt(num))
# for i in range(0,k):
# n_col = 10
#     row = np.where(fit==i)[0]
#     print(row) # row in Z for elements of cluster i
#     num = row.shape[0]       #  number of elements for each cluster
n_row = int(np.floor(num/n_col)+1)    # number of rows in the figure of the cluster
#     print("cluster "+str(i))
#     print(str(num)+" elements")
fig = plt.figure(figsize=(n_col/10+2, n_row/10+2))
dx = float(L)/float(J-1)
grid = np.array([j*dx for j in range(J)])
for count,parID in tqdm(enumerate(parID_list),disable=False):
    ax=plt.subplot(n_row,n_col, count+1)
    #rgb_timeseries=timeseries_unstacked_list[row[n]] # Read the numpy matrix with images in the rows
    # par_ID = parID_list[row[n]]
    # parID = parID_list[count]
    filename = 'circuit%r_variant%s_bc%s_%s_%sID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundarycoeff, shape,mechanism,parID,L,J,T,N)
    # final_concentration = pickle.load( open(modelling_ephemeral + '/3954/numerical_confocal/results/simulation/1M_colony_ca/2D/full_circuit_newCN/2Dfinal_%s.pkl'%filename, 'rb' ) )
    final_concentration = pickle.load( open(data_path + '/2Dfinal_%s.pkl'%filename, 'rb' ) )

    # ax.pcolormesh(grid, grid, final_concentration[2], shading='auto')
    # rgb = plot_redgreen_contrast(final_concentration,L,mechanism,shape,filename,modelling_ephemeral,parID=parID,dimension=dimension,scale_factor=x_gridpoints,save_figure='LargeImage')
    mask=pickle.load( open( modelling_home + "/3954/numerical_confocal/code/cellular_automata_templates/masks/caMask_seed%s_pdivision%s_L%s_J%s_T%s_N%s_fast.pkl"%(seed,p_division,L,J,T,N), "rb" ) )
    rgb = plot_redgreenblue_contrast(final_concentration,L,mechanism,shape,filename,parID=parID,mask=mask,scale_factor=x_gridpoints,save_figure='LargeImage')

    # # rgb_timeseries=timeseries_unstacked # Read the numpy matrix with images in the rows
    # ax.set_title(parID,size=0.1)
    ax.set(yticks=[],xticks=[],yticklabels=[],xticklabels=[])
    ax.imshow(rgb.astype('uint8'), origin= 'lower')
    ax.set_ylabel(parID,size= 1,c='y', labelpad=0.35)



# plt.title('1M numerical search 0-%r'%num)
filename = 'circuit%r_variant%s_bc%s_%s_%s_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundarycoeff, shape,mechanism,L,J,T,N)
# plt.savefig(modelling_home + '/3954/numerical_confocal/results/figures/%s/large_images/%s_%s-%s.png'%(shape,filename,start,stop), dpi=2000)
plt.savefig(modelling_home + '/3954/numerical_confocal/results/figures/%s/large_images/%s_%s.png'%(shape,filename,details), dpi=2000)

# plt.savefig('h.png',dpi=2000)
print('gh')
# plt.clf()
# plt.close(fig)
