#%%
#############
###paths#####
#############
import os
import sys


pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
modellingephemeral = '/rds/general/ephemeral/user/mo2016/ephemeral/Documents/modelling'
modellingephemeral = '/rds/general/ephemeral/user/mo2016/ephemeral/Documents/modelling'
sys.path.append(modellingpath + '/lib')

import pickle
import time

import matplotlib.pyplot as plt
import numpy as np

#############
###Imports#####
#############
from numerical.adi_ca_function_openclosed_nodilution import \
    adi_ca_openclosed_nodilution
from numerical.adi_ca_function_openclosed_nodilution_preMask import \
    adi_ca_openclosed_nodilution_preMask
from numerical.adi_ca_function_openclosed_nodilution_preMask_numba import \
    adi_ca_openclosed_nodilution_preMask as \
    adi_ca_openclosed_nodilution_preMask_numba
from numerical.adi_square_function import adi
from numerical.plotting_numerical import *
from numerical.cn_plot import *
from colonyMaskCreation import *
#%%
#############
###execution parameters#####
#############
# %matplotlib inline
shape = 'ca'
# circuit_n=14;variant='2nd';n_species=6
circuit_n=14;variant='fitted7';n_species=6
# Specifiy number of parameter sets in parameterset file to be loaded
n_samples = 13700000
# balance = 'balanced'
# folder = 'circuit14variant2ndBalancedTuring'
folder = f'circuit{circuit_n}variant{variant}'

# folder = 'circuit14variant2nd_turing'
# folder = 'circuit14fitted3balancedSemibalanced'
save_figure = False
tqdm_disable = False #disable tqdm
# boundarycoeff = float(sys.argv[6])


# open parameter dictionaries
df= pickle.load( open(modellingpath + "/3954/paper/input/fitted_parameterfiles/df_circuit%s_variant%s_%rparametersets.pkl"%(circuit_n,variant,n_samples), "rb"))

# L=20; dx =0.1; J = int(L/dx)
# T =50; dt = 0.1; N = int(T/dt)
# # T =50; dt = 0.1; N = int(T/dt)

# L=30; dx =0.2; J = int(L/dx)
# T =200; dt = 0.1; N = int(T/dt)

# boundarycoeff = 1
# divisionTimeHours=1
# p_division=0.25;seed=1





# L=20; dx =0.1; J = int(L/dx)
# T =50; dt = 0.3; N = int(T/dt)
# boundarycoeff = 1
# divisionTimeHours=0.3
# p_division=0.5;seed=1


L=20; dx =0.1; J = int(L/dx)
T =50; dt = 0.02; N = int(T/dt)
boundarycoeff = 1
divisionTimeHours=0.5
p_division=1;seed=1



shape = 'ca'
x_gridpoints=int(1/dx)

try:
    cell_matrix_record = pickle.load( open(modellingpath + "/3954/paper/out/numerical/masks/caMask_seed%s_pdivision%s_L%s_J%s_T%s_N%s.pkl"%(seed,p_division,L,J,T,N), "rb" ) )
    daughterToMotherDictList = pickle.load( open(modellingpath + "/3954/paper/out/numerical/masks/caMemory_seed%s_pdivision%s_L%s_J%s_T%s_N%s.pkl"%(seed,p_division,L,J,T,N), "rb" ) )
    print('fileNotCreated')

except:
    #file does not exist
    FileNotFoundError
    print('fileCreation')

    maskFunction(L=L,dx=dx, T=T, dt=dt, divisionTimeHours=divisionTimeHours, p_division=p_division, plot1D=True, plotScatter=True)
    cell_matrix_record = pickle.load( open(modellingpath + "/3954/paper/out/numerical/masks/caMask_seed%s_pdivision%s_L%s_J%s_T%s_N%s.pkl"%(seed,p_division,L,J,T,N), "rb" ) )
    daughterToMotherDictList = pickle.load( open(modellingpath + "/3954/paper/out/numerical/masks/caMemory_seed%s_pdivision%s_L%s_J%s_T%s_N%s.pkl"%(seed,p_division,L,J,T,N), "rb" ) )
    

# maskFunction(L=L,dx=dx, T=T, dt=dt, divisionTimeHours=divisionTimeHours, p_division=p_division, plot1D=True, plotScatter=True)
# cell_matrix_record = pickle.load( open(modellingpath + "/3954/paper/out/numerical/masks/caMask_seed%s_pdivision%s_L%s_J%s_T%s_N%s.pkl"%(seed,p_division,L,J,T,N), "rb" ) )
# daughterToMotherDictList = pickle.load( open(modellingpath + "/3954/paper/out/numerical/masks/caMemory_seed%s_pdivision%s_L%s_J%s_T%s_N%s.pkl"%(seed,p_division,L,J,T,N), "rb" ) )

filename= lambda parID: 'circuit%r_variant%s_bc%s_%s_ID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundarycoeff, shape,parID,L,J,T,N)
# parID=12837401
parID=12837401


#%%
# test = bool(sys.argv[7])
test = False
tqdm_disable=False
if test==True:
    print('test')
    T=1;N =1
    tqdm_disable=False


# degDiv = 1
# par_dict['muASV'] =par_dict['muASV']/degDiv
# par_dict['muLVA'] = par_dict['muLVA'] /degDiv
# degDiv = 1
# par_dict['muASV'] =par_dict['muASV']/degDiv

# 195238
print('parID = ' + str(parID))
par_dict = df.loc[parID].to_dict()
# paramList = [1,1.5,2]

par_dict = df.loc[parID].to_dict()
# par_dict = df.loc[(parID,0)].to_dict()
print(par_dict)

D = np.zeros(n_species)
Dr = float(par_dict['Dr'])
D[:2] = [1,Dr ]

print(par_dict)


# U_record,U_final =  adi_ca_openclosed_nodilution_preMask(par_dict,L,dx,J,T,dt,N, circuit_n, n_species,D,cell_matrix_record, daughterToMotherDictList,tqdm_disable=False, p_division=0.5,stochasticity=0, seed=1,growth='Slow', boundarycoeff=boundarycoeff)
# U_record,U_final =  adi(par_dict,L,L,J,J,T,N, circuit_n, n_species,D,tqdm_disable=False,stochasticity=0, steadystates=0)
# get the start time
st = time.time()
U_record,U_final =  adi_ca_openclosed_nodilution_preMask_numba(par_dict,L,dx,J,T,dt,N, circuit_n, n_species,D,cell_matrix_record, daughterToMotherDictList,tqdm_disable=tqdm_disable,divisionTimeHours=divisionTimeHours, stochasticity=0, seed=1, boundarycoeff=boundarycoeff)
elapsed_time = time.time() - st
print('Execution time numba:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
# plt.imshow(U_final[-1])
# plt.show()



pickle.dump(U_final, open(modellingpath + '/3954/paper/out/numerical/colonies/simulation/%s/2Dfinal_%s.pkl'%(folder,filename(parID)), "wb" ) )
pickle.dump(U_record, open(modellingpath + '/3954/paper/out/numerical/colonies/simulation/%s/2Drecord_%s.pkl'%(folder,filename(parID)), 'wb'))

print('saved')
saveFigPath = modellingpath + '/3954/paper/out/numerical/colonies/figures/%s'%folder

rgb = plot_redgreen_contrast(U_final,L,parID=filename(parID),filename = filename(parID), path = saveFigPath,scale_factor=x_gridpoints,save_figure=True)

# %%
# plt.imshow(U_final[-1])
# plt.colorbar()
# plt.imshow(U_final[-1])
# plt.colorbar()

# %%



U_final = pickle.load(open(modellingpath + '/3954/paper/out/numerical/colonies/simulation/%s/2Dfinal_%s.pkl'%(folder,filename(parID)), "rb" ) )
U_record = pickle.load(open(modellingpath + '/3954/paper/out/numerical/colonies/simulation/%s/2Drecord_%s.pkl'%(folder,filename(parID)), 'rb'))
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
# plt.rcParams['animation.ffmpeg_path'] = '~/Documents/virtualEnvironments/env1/lib/python3.8/site-packages/ffmpeg'
# plt.rcParams['animation.ffmpeg_path'] = '/usr/local/bin/'
rgb_timeseries = redgreen_contrast_timeseries(U_record)
# show_rgbvideo(rgb_timeseries,parID)
saveVideoPath = modellingpath + '/3954/paper/out/numerical/colonies/videos/%s/'%folder

def save_rgbvideo(timeseries_unstacked, saveVideoPath, filename, interval=10000):
    fig = plt.figure()
    ims = []
    rgb_timeseries=timeseries_unstacked # Read the numpy matrix with images in the rows
    im=plt.imshow(rgb_timeseries[0].astype('uint8'), origin= 'lower')

    for i in range(len(rgb_timeseries)):
        im=plt.imshow(rgb_timeseries[i].astype('uint8'), origin= 'lower')
        plt.title(str(filename) + str(i))
        plt.xlabel(f'Time: {i}h')
        
        ims.append([im])
    ani = animation.ArtistAnimation(fig, ims,interval=50000000)
    
    # FOR GIF
    # writergif = animation.PillowWriter(fps=10)
    # ani.save(saveVideoPath + filename + '.gif',writer=writergif)
    
    # FOR MP4
    mywriter = animation.FFMpegWriter()
    ani.save(saveVideoPath + '/%s.mp4' %filename,writer=mywriter)
    print('Video saved', filename)

save_rgbvideo(rgb_timeseries, saveVideoPath, filename(parID))
print('Video saved', filename(parID))
# %%
import seaborn as sns 
# def black_to_green_palette(num_colors):
#     cmap = sns.color_palette("b", "g", num_colors)  # Generate color palette
#     return cmap
cm = sns.color_palette('blend:green,black', as_cmap=True)
plt.imshow(-U_final[-1], cmap=cm)
plt.show()
cm = sns.color_palette('blend:red,black', as_cmap=True)
plt.imshow(-U_final[-2], cmap=cm)
plt.show()

    # %%
