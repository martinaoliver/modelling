
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
from numerical_solvers_variableboundary import plot_redgreen_contrast
from tqdm import tqdm
#Opening list with parID's
# file = open(modelling_ephemeral + '/3954/numerical_confocal/results/simulation/1M_colony_ca/2D/parID_list_8x10T120.txt')
parID_list = pickle.load( open(modelling_ephemeral + '/3954/numerical_confocal/results/simulation/1M_colony_ca/2D/nodeB_dele/parID_list_8x10T120.pkl', "rb" ) )
print(parID_list)

start = int(sys.argv[1])
stop = int(sys.argv[2])
parID_list = [int(i) for i in parID_list[start:stop]] #turn string list into integer list
circuit_n=9
variant=0
boundary_coef=1
shape='ca'
mechanism = 'nodeBdele'
L,T,J,N = [8,120,80,23880]
dimension='2D'
x_gridpoints=10
# k=20
num=len(parID_list)
print(num)
n_col = int(np.sqrt(num))
# for i in range(0,k):
# n_col = 10
#     row = np.where(fit==i)[0]
#     print(row) # row in Z for elements of cluster i
#     num = row.shape[0]       #  number of elements for each cluster
n_row = np.floor(num/n_col)+1    # number of rows in the figure of the cluster
#     print("cluster "+str(i))
#     print(str(num)+" elements")
fig = plt.figure(figsize=(n_col/10+2, n_row/10+2))
for count,n in enumerate(parID_list):
    print(n)
    # print(count)
    ax=plt.subplot(n_row,n_col, count+1)
    #rgb_timeseries=timeseries_unstacked_list[row[n]] # Read the numpy matrix with images in the rows
    # par_ID = parID_list[row[n]]
    parID = parID_list[count]
    filename = 'circuit%r_variant%r_boundary%r_%s_%sID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundary_coef, shape,mechanism,parID,L,J,T,N)
    final_concentration = pickle.load( open(modelling_ephemeral + '/3954/numerical_confocal/results/simulation/1M_colony_ca/2D/nodeB_dele/2Dfinal_%s.pkl'%filename, 'rb' ) )
    rgb = plot_redgreen_contrast(final_concentration,L,mechanism,shape,filename,modelling_ephemeral,parID=parID,dimension=dimension,scale_factor=x_gridpoints,save_figure='results')
    # rgb_timeseries=timeseries_unstacked # Read the numpy matrix with images in the rows
    # ax.set_title(parID,size=0.1)
    ax.set(yticks=[],xticks=[],yticklabels=[],xticklabels=[])
    ax.imshow(rgb.astype('uint8'), origin= 'lower')
    ax.set_ylabel(parID,size= 1,c='y', labelpad=0.35)


# plt.title('1M numerical search 0-%r'%num)
plt.savefig(modelling_home + '/3954/numerical_confocal/results/figures/redgreen/large_images/1M_numerical_search_%s-%s_%s.png'%(start,stop,mechanism), dpi=2000)
# plt.show()
# plt.savefig('h.png',dpi=2000)
print('gh')
# plt.clf()
# plt.close(fig)
