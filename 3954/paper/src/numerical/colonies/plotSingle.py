#%%
#############
###paths#####
#############
import os
import sys


pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
modellingephemeral = '/rds/general/ephemeral/user/mo2016/ephemeral/Documents/modelling'
sys.path.append(modellingpath + '/lib')

import pickle
import time

import matplotlib.pyplot as plt
import numpy as np
#############
###Imports#####
#############
from numerical.plotting_numerical import *
from numerical.cn_plot import *
#############
###execution parameters#####
#############

# Specify name of circuit and variant investigated
circuit_n=14;variant='2nd';n_species=6
# Specifiy number of parameter sets in parameterset file to be loaded
n_param_sets = 1000000
# balance = 'balanced'
folder = 'circuit14variant2ndBalancedTuring'
modelArgs = [circuit_n,variant,n_species,folder]

# Specifiy number of parameter sets in parameterset file to be loaded
nsamples = 1000000

# specify dimensions of system
L=9; dx =0.05; J = int(L/dx)
T =50; dt = 0.05; N = int(T/dt)
boundarycoeff = 1
divisionTimeHours=0.5
p_division=1;seed=1



shape = 'ca'
x_gridpoints=int(1/dx)



save_figure = False


parID=20766
# parID = 651894
# parID_list = [89407, 706280]
# parID = parID_list[1]
# parID = 62509

#%%
#load image
# U_final = pickle.load( open(modellingpath + '/3954/paper/out/numerical/colonies/simulation/%s/2Dfinal_%s.pkl'%(folder,filename(parID)), 'rb'))
degDiv=4
filename= lambda parID: 'circuit%r_variant%s_bc%s_%s_ID%r_L%r_J%r_T%r_N%r._degDiv%s'%(circuit_n,variant,boundarycoeff, shape,parID,L,J,T,N,degDiv)

# U_record = pickle.load( open(modellingpath + '/3954/paper/out/numerical/colonies/simulation/%s/2Drecord_%s.pkl'%(folder,filename(parID)), 'rb'))
U_record = pickle.load( open(modellingephemeral + '/3954/paper/out/numerical/colonies/simulation/%s/2Drecord_%spkl'%(folder,filename(parID)), 'rb'))
# U_record = pickle.load( open(modellingpath + '/3954/paper/out/numerical/colonies/simulation/%s/2Drecord_%s.pkl'%(folder,filename(parID)), 'rb'))



# pickle.dump(U_final, open(modellingpath + '/3954/paper/out/numerical/colonies/simulation/%s/2Dfinal_%s.pkl'%(folder,filename(parID)), "wb" ) )
# pickle.dump(U_record, open(modellingpath + '/3954/paper/out/numerical/colonies/simulation/%s/2Drecord_%s.pkl'%(folder,filename(parID)), 'wb'))


savefig_path  = ''
# rgb = plot_redgreen_contrast(U_final,L,parID=parID,scale_factor=x_gridpoints,save_figure='LargeImage')


# plt.imshow(rgb.astype('uint8'), origin= 'lower')


#%%
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
# plt.rcParams['animation.ffmpeg_path'] = '~/Documents/virtualEnvironments/env1/lib/python3.8/site-packages/ffmpeg'
plt.rcParams['animation.ffmpeg_path'] = '/usr/local/bin/'
rgb_timeseries = redgreen_contrast_timeseries(U_record)
# show_rgbvideo(rgb_timeseries,parID)
saveVideoPath = modellingpath + '/3954/paper/out/numerical/colonies/videos/%s/'%folder


def save_rgbvideo(timeseries_unstacked, saveVideoPath, filename, interval=100):
    fig = plt.figure()
    ims = []
    rgb_timeseries=timeseries_unstacked # Read the numpy matrix with images in the rows
    im=plt.imshow(rgb_timeseries[0].astype('uint8'), origin= 'lower')

    for i in range(len(rgb_timeseries)):
        im=plt.imshow(rgb_timeseries[i].astype('uint8'), origin= 'lower')
        plt.title(str(filename) + str(i))
        plt.xlabel(f'Time: {i}h')
        ims.append([im])
    ani = animation.ArtistAnimation(fig, ims, interval=interval)
    
    # ani.save(saveVideoPath + '/%s.mp4' %filename)
    print('Video saved')
    
    #FOR GIF
    writergif = animation.PillowWriter(fps=100)
    ani.save(saveVideoPath + filename + '.gif',writer=writergif)

    # FOR MP4
    # mywriter = animation.FFMpegWriter()
    # ani.save('mymovie.mp4',writer=mywriter)
    
save_rgbvideo(rgb_timeseries, saveVideoPath, filename(parID))




# U_record[-1][:,:,-1].shape
# plt.plot(U_record[-2][:,40,-1],C='r')
# plt.plot(U_record[-1][:,40,-1],C='g')



# plot1D([U_record[-2][:,40,-1], U_record[-1][:,40,-1]])




# %%
