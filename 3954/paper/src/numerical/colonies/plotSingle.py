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
circuit_n=14;variant='1nd';n_species=6
# Specifiy number of parameter sets in parameterset file to be loaded
n_param_sets = 1000000
balance = 'balanced'
folder = 'circuit14variant1nd_turing'
modelArgs = [circuit_n,variant,n_species,folder]

# Specifiy number of parameter sets in parameterset file to be loaded
nsamples = 1000000

# specify dimensions of system
L=4; dx =0.05; J = int(L/dx)
T =125; dt = 0.05; N = int(T/dt)
boundarycoeff = 1
shape = 'ca'

divisionTimeHours=1
p_division=0.22;seed=1
x_gridpoints=int(1/dx)



save_figure = False
nsamples=1000000


parID = 990343

filename= lambda parID: 'circuit%r_variant%s_bc%s_%s_ID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundarycoeff, shape,parID,L,J,T,N)


#%%

#load image
U_final = pickle.load( open(modellingephemeral + '/3954/paper/out/numerical/colonies/simulation/%s/2Dfinal_%s.pkl'%(folder,filename(parID)), 'rb'))
U_record = pickle.load( open(modellingephemeral + '/3954/paper/out/numerical/colonies/simulation/%s/2Drecord_%s.pkl'%(folder,filename(parID)), 'rb'))
# plt.imshow(U_final[-1])
# plt.colorbar()
# plt.show()
# plt.imshow(U_final[-2])
# plt.colorbar()
# plt.show()
savefig_path  = ''
# savefig_path = modellingpath + '/3954/paper/out/numerical/colonies/figures/%s'%(folder)
# rgb = plot_redgreen_contrast(U_final,L,path = savefig_path,filename=filename(parID),parID=parID,scale_factor=int(1/dx),save_figure='Large ')
rgb = plot_redgreen_contrast(U_final,L,parID=parID,scale_factor=x_gridpoints,save_figure='LargeImage')
# def plot_redgreen_contrast(final_concentration, mm,filename=None, path=None, parID=0, scale_factor=10, save_figure=False, dimension='2D'):


# plt.set(yticks=[],xticks=[],yticklabels=[],xticklabels=[])
plt.imshow(rgb.astype('uint8'), origin= 'lower')
# ax.set_ylabel(parID,size= 1,c='y', labelpad=0.35)


#%%

rgb_timeseries = redgreen_contrast_timeseries(U_record)
# show_rgbvideo(rgb_timeseries,parID)
saveVideoPath = modellingephemeral + '/3954/paper/out/numerical/colonies/videos/%s'%folder
save_rgbvideo(rgb_timeseries, saveVideoPath, filename(parID))




# U_record[-1][:,:,-1].shape
# plt.plot(U_record[-2][:,40,-1],C='r')
# plt.plot(U_record[-1][:,40,-1],C='g')



plot1D([U_record[-2][:,40,-1], U_record[-1][:,40,-1]])



