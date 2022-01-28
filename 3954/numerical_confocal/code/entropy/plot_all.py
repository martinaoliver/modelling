#############################
#IMPORTS#
#############################
import sys
import os
pwd = os.getcwd()
root = pwd.rpartition("mo2016")[0] + pwd.rpartition("mo2016")[1] #/Volumes/mo2016/ or '/Users/mo2016/' or '/rds/general/mo2016/'


if root == '/Users/mo2016':
    modelling_ephemeral = '/Volumes/mo2016/ephemeral/Documents/modelling'
    modelling_home = '/Volumes/mo2016/home/Documents/modelling'
    modelling_local = root + '/Documents/modelling'
    modelling_local = root + '/Documents/modelling'
    import matplotlib as mpl
    mpl.use('tkagg')

if root == '/Volumes/mo2016' or '/rds/general': #'/rds/general' or root=='/Volumes':
        modelling_ephemeral = root + '/ephemeral/Documents/modelling'
        modelling_home = root  + '/home/Documents/modelling'
        modelling_local = modelling_home

modulepath = modelling_local + '/3954/modules/new_CN'
sys.path.append(modulepath)


from plotting_numerical import plot_redgreen_contrast
from tqdm import tqdm
import pickle
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib.animation as animation
import matplotlib as mpl

#############################
# df= pickle.load( open(modelling_home + '/3954/parameter_space_search/parameterfiles/df_circuit%r_variant%s_%rparametersets.pkl'%(2,0,1000000), "rb" ) )

#Opening list with parID's
# file = open(modelling_ephemeral + '/3954/numerical_confocal/results/simulation/1M_colony_ca/2D/parID_list_8x10T120.txt')
parID_dict = pickle.load( open(modelling_home + '/3954/numerical_confocal/results/entropy/EntropyDicts/EntropyDictdictGreen_v1_variant0_ca_fullcircuit_L10J150T120N1200.pkl', "rb" ) )
len(parID_dict)
parID_list = []
entropy_list = []
for key in parID_dict:
    parID_list.append(key)
    entropy_list.append(parID_dict[key])
parID_list
start = int(sys.argv[1])
stop = int(sys.argv[2])
# parID_list = [int(i) for i in parID_list[start:stop]] #turn string list into integer list
parID_list = [int(i) for i in parID_list] #turn string list into integer list
circuit_n=2
variant='0'
boundary_coef=1
shape='ca'
mechanism = 'fullcircuit'
L,J,T,N = [10,150,120,1200]
dimension='2D'
x_gridpoints=L
details = 'EntropyDictdictGreen_v1_variant0_ca_fullcircuit_L10J150T120N1200'
num=len(parID_list)
print(num)
n_col = int(np.sqrt(num))
n_row = np.floor(num/n_col)+1    # number of rows in the figure of the cluster

print(parID_list)
fig = plt.figure(figsize=(n_col/10+2, n_row/10+2))
for count,n in tqdm(enumerate(parID_list)):
    # print(count)
    ax=plt.subplot(n_row,n_col, count+1)
    #rgb_timeseries=timeseries_unstacked_list[row[n]] # Read the numpy matrix with images in the rows
    # par_ID = parID_list[row[n]]
    parID = parID_list[count]
    filename = '2Dfinal_circuit%r_variant%s_%s_%sID%s_L%r_J%r_T%r_N%r'%(circuit_n,variant,shape,mechanism, parID,L,J,T,N)
    final_concentration = pickle.load( open(modelling_home + '/3954/numerical_confocal/results/simulation/1M_colony_ca/2D/full_circuit_newCN/%s.pkl'%filename, 'rb' ) )
    rgb = plot_redgreen_contrast(final_concentration,L,mechanism,shape,filename,modelling_ephemeral,parID=parID,scale_factor=x_gridpoints,save_figure='LargeImage')
    # rgb_timeseries=timeseries_unstacked # Read the numpy matrix with images in the rows
    # ax.set_title(parID,size=0.1)
    ax.set(yticks=[],xticks=[],yticklabels=[],xticklabels=[])
    ax.imshow(rgb.astype('uint8'), origin= 'lower')
    ax.set_ylabel(parID,size= 1,c='y', labelpad=0.35)




# plt.title('1M numerical search 0-%r'%num)
plt.savefig(modelling_home + '/3954/numerical_confocal/results/entropy/LargeImages/1M_numerical_search_%s-%s_%s.png'%(start,stop,details), dpi=2000)
# plt.show()
# plt.savefig('h.png',dpi=2000)
print('gh')
# plt.clf()
# plt.close(fig)
