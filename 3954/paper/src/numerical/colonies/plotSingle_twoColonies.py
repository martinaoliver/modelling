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
from colonyMaskCreation_twoColonies import *
#%%
#############
###execution parameters#####
#############
# %matplotlib inline
shape = 'caTwoColonies'

# shape = 'caTwoColonies'
circuit_n=14;variant='2nd';n_species=6
# circuit_n=14;variant='fitted1';n_species=6
# Specifiy number of parameter sets in parameterset file to be loaded
n_param_sets = 10
balance = 'balanced'
# folder = 'circuit14variantfitted1'
folder = 'circuit14variant2ndBalancedTuring'
# nsamples =  2000000
nsamples =  1000000
save_figure = False
tqdm_disable = False #disable tqdm
# boundarycoeff = float(sys.argv[6])


# open parameter dictionaries


df= pickle.load( open(modellingpath + '/3954/paper/input/balanced_parameterfiles/df_circuit%r_variant%s_%rparametersets_balanced.pkl'%(circuit_n,variant,nsamples), "rb" ) )
# df= pickle.load( open(modellingpath + '/3954/paper/input/fitted_parameterfiles/df_circuit%r_variant%s_%rparametersets.pkl'%(circuit_n,variant,nsamples), "rb" ) )
# instabilities_df= pickle.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/instabilities_dataframes/instability_df_circuit%s_variant%s_%rparametersets.pkl'%(circuit_n,variant,nsamples), "rb" ) )
# with open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/turing_dataframes/turing_df_circuit%s_variant%s_%rparametersets_balanced.pkl'%(circuit_n,variant,nsamples), "rb" ) as f:
#     df = pickle.load(f)
#solver parameters
# specify dimensions of system

L=25; dx =0.1; J = int(L/dx)
T =50; dt = 0.02; N = int(T/dt)
T =50; dt = 0.5; N = int(T/dt)
# T =3; dt = 0.5; N = int(T/dt)
boundarycoeff = 1
divisionTimeHours=0.5
p_division=1;seed=1

x_gridpoints=int(1/dx)
parID=195238

try:
    cell_matrix_record = pickle.load( open(modellingpath + "/3954/paper/out/numerical/masks/%sMask_seed%s_pdivision%s_L%s_J%s_T%s_N%s.pkl"%(shape,seed,p_division,L,J,T,N), "rb" ) )
    daughterToMotherDictList = pickle.load( open(modellingpath + "/3954/paper/out/numerical/masks/%sMemory_seed%s_pdivision%s_L%s_J%s_T%s_N%s.pkl"%(shape,seed,p_division,L,J,T,N), "rb" ) )
    print('fileNotCreated')

except:
    #file does not exist
    FileNotFoundError
    print('fileCreation')

    cell_matrix_record, memory_matrix_record, daughterToMotherDictList = maskFunction_twoColonies(L=L,dx=dx, T=T, dt=dt, divisionTimeHours=divisionTimeHours, p_division=p_division, plot1D=True, plotScatter=True)
    pickle.dump( cell_matrix_record,open(modellingpath + "/3954/paper/out/numerical/masks/%sMask_seed%s_pdivision%s_L%s_J%s_T%s_N%s.pkl"%(shape,seed,p_division,L,J,T,N), "wb" ) )
    pickle.dump( daughterToMotherDictList, open(modellingpath + "/3954/paper/out/numerical/masks/%sMemory_seed%s_pdivision%s_L%s_J%s_T%s_N%s.pkl"%(shape,seed,p_division,L,J,T,N), "wb" ) )

    # cell_matrix_record = pickle.load( open(modellingpath + "/3954/paper/out/numerical/masks/caTwoColoniesMask_seed%s_pdivision%s_L%s_J%s_T%s_N%s.pkl"%(seed,p_division,L,J,T,N), "rb" ) )
    # daughterToMotherDictList = pickle.load( open(modellingpath + "/3954/paper/out/numerical/masks/caTwoColoniesMemory_seed%s_pdivision%s_L%s_J%s_T%s_N%s.pkl"%(seed,p_division,L,J,T,N), "rb" ) )

filename= lambda parID: 'circuit%r_variant%s_bc%s_%s_ID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundarycoeff, shape,parID,L,J,T,N)
#%%
# test = bool(sys.argv[7])
test = False
tqdm_disable=False
if test==True:
    print('test')
    T=1;N = int(T/dt); N=3
    tqdm_disable=False
# parID=int(sys.argv[8])
print('parID = ' + str(parID))
par_dict = df.loc[parID].to_dict()
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
plt.imshow(U_final[-1])
plt.show()



pickle.dump(U_final, open(modellingpath + '/3954/paper/out/numerical/colonies/simulation/%s/2Dfinal_%s.pkl'%(folder,filename(parID)), "wb" ) )
pickle.dump(U_record, open(modellingpath + '/3954/paper/out/numerical/colonies/simulation/%s/2Drecord_%s.pkl'%(folder,filename(parID)), 'wb'))


print('saved')

rgb = plot_redgreen_contrast(U_final,L,parID=parID,scale_factor=x_gridpoints,save_figure=False)
# rgb = plot_redgreen_contrast(U_final,L,parID=parID,scale_factor=x_gridpoints,save_figure=False)



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
#snapshots

plt.imshow(rgb_timeseries[1].astype('uint8'), origin= 'lower')
# %%
