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
# circuit_n=14;variant='fitted1';n_species=6
# circuit_n=14;variant='2nd';n_species=6
circuit_n=14;variant='2nd';n_species=6
# circuit_n=14;variant=195238;n_species=6;nsr=0.05

# Specifiy number of parameter sets in parameterset file to be loaded
n_param_sets = 1000000
# balance = 'balanced'\
# circuit_n=14;variant='fitted1';n_species=6

# folder = 'circuit14variantfitted1'
# folder = 'circuit14variant2ndBalancedTuring'
folder = 'circuit14variant2ndBalancedKce100'
# folder = f'circuit14variant{variant}'

modelArgs = [circuit_n,variant,n_species,folder]

# Specifiy number of parameter sets in parameterset file to be loaded
nsamples =  1000000

# specify dimensions of system

L=20; dx =0.1; J = int(L/dx)
T =25; dt = 0.02; N = int(T/dt)
boundarycoeff = 1
divisionTimeHours=0.2
p_division=0.7;seed=1

shape = 'ca'
x_gridpoints=int(1/dx)



save_figure = False

parID=361095
# parID=94



#%%
filename= lambda parID: 'circuit%r_variant%s_%sparametersets_balanced_Kce%s_bc%s_%s_ID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant,nsamples,Kce,boundarycoeff, shape,parID,L,J,T,N)

# filename= lambda parID: 'circuit%r_variant%snsr%s_bc%s_%s_ID%s_L%r_J%r_T%r_N%r'%(circuit_n,variant,nsr,boundarycoeff, shape,parID,L,J,T,N)



U_record = pickle.load( open(modellingpath + '/3954/paper/out/numerical/colonies/simulation/%s/2Drecord_%s.pkl'%(folder,filename(parID)), 'rb'))
U_final = pickle.load( open(modellingpath + '/3954/paper/out/numerical/colonies/simulation/%s/2Dfinal_%s.pkl'%(folder,filename(parID)), 'rb'))
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
    writergif = animation.PillowWriter(fps=10)
    ani.save(saveVideoPath + filename + '.gif',writer=writergif)
    
    # FOR MP4
    # mywriter = animation.FFMpegWriter()
    # ani.save(saveVideoPath + '/%s.mp4' %filename,writer=mywriter)
    # print('Video saved', filename)

save_rgbvideo(rgb_timeseries, saveVideoPath, filename(parID))
print('Video saved', filename(parID))

#%%
#load image
# U_final = pickle.load( open(modellingpath + '/3954/paper/out/numerical/colonies/simulation/%s/2Dfinal_%s.pkl'%(folder,filename(parID)), 'rb'))
# filename= lambda parID: 'circuit%r_variant%s_bc%s_%s_ID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundarycoeff, shape,parID,L,J,T,N)
# perturbation = 'Kda2'
filename= lambda parID: 'circuit%r_variant%s_%sparametersets_balanced_Kce%s_bc%s_%s_ID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant,nsamples,Kce,boundarycoeff, shape,parID,L,J,T,N)

# filename= lambda parID: 'circuit%r_variant%snsr%s_bc%s_%s_ID%s_L%r_J%r_T%r_N%r'%(circuit_n,variant,nsr,boundarycoeff, shape,parID,L,J,T,N)



U_record = pickle.load( open(modellingpath + '/3954/paper/out/numerical/colonies/simulation/%s/2Drecord_%s.pkl'%(folder,filename(parID)), 'rb'))
U_final = pickle.load( open(modellingpath + '/3954/paper/out/numerical/colonies/simulation/%s/2Dfinal_%s.pkl'%(folder,filename(parID)), 'rb'))


savefig_path  = ''
# rgb = plot_redgreen_contrast(U_final,L,parID=parID,scale_factor=x_gridpoints,save_figure='LargeImage')


# plt.imshow(rgb.astype('uint8'), origin= 'lower')



##%



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

#%%



#%%
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
# plt.rcParams['animation.ffmpeg_path'] = '~/Documents/virtualEnvironments/env1/lib/python3.8/site-packages/ffmpeg'
# plt.rcParams['animation.ffmpeg_path'] = '/usr/local/bin/'
plt.rcParams['animation.ffmpeg_path'] = '/rds/general/user/mo2016/home/anaconda3/envs/env1/lib/python3.8/site-packages'

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
    ani = animation.ArtistAnimation(fig, ims,interval=50000000, repeat=False)
    
    # ani.save(saveVideoPath + '/%s.mp4' %filename)
    print('Video saved')
    
    # FOR GIF
    writergif = animation.PillowWriter(fps=10)
    ani.save(saveVideoPath + filename + '.gif',writer=writergif)
    ani.save('mygif' + '.gif',writer=writergif)

    # # FOR MP4
    # mywriter = animation.FFMpegWriter()
    # ani.save('mymovie.mp4',writer=mywriter)
    
save_rgbvideo(rgb_timeseries, saveVideoPath, filename(parID))
import moviepy.editor as mp

clip = mp.VideoFileClip(saveVideoPath + filename + '.gif')
clip.write_videofile(saveVideoPath + filename + '.mp4')

# %%
def plot1D(U_final, savefig=False,filename='',savefigpath='',pad=0.001):
    cutGreenUfinal = U_final[-1][int(J/2)]
    cutRedUfinal = U_final[-2][int(J/2)]

    fig,ax = plt.subplots()
    ax.plot(cutGreenUfinal, color='green')
    ax.set_ylim(np.amin(cutGreenUfinal)-pad, np.amax(cutGreenUfinal)+pad)
    ax.ticklabel_format(useOffset=False)
    ax.set_ylabel('GFP', fontsize=15)

    ax2=ax.twinx()
    ax2.plot(cutRedUfinal,color='red')
    ax2.set_ylim(np.amin(cutRedUfinal)-pad, np.amax(cutRedUfinal)+pad)
    ax2.set_ylabel('mCherry', fontsize=15)

    ax.ticklabel_format(useOffset=False)
    ax.tick_params(axis='x', labelsize=15)
    ax.tick_params(axis='y', labelsize=15)
    ax2.tick_params(axis='y', labelsize=15)
    ax.set_xlabel('Space', fontsize=15)
    if savefig==True:
        plt.savefig('%s%s.jpeg'%(savefigpath,filename))
    else:
        plt.show()

plot1D(U_final)

#%%

x_grid = np.array([j*dx for j in range(J)])
t_grid = np.array([n*dt for n in range(T)])

grids = [x_grid,t_grid]

def matrix_rgb_normalisation(matrix):
    row_n = 0
    NewMatrix = np.zeros(matrix.shape)

    OldMin = np.min(matrix[np.nonzero(matrix)])+1e-10
    OldMax = np.amax(matrix[np.nonzero(matrix)]) #Add 0.0001 so that patterns with no var dont give errors

    OldMin = np.min(matrix)+1e-10
    OldMax = np.amax(matrix) #Add 0.0001 so that patterns with no var dont give errors


    if OldMin < 0:
        print('WARNING: Negative numbers!!!!!')
    # NewMin = 1 #make newmin 1 instead of zero so 1 can represent cells
    NewMin = 0
    NewMax = 255
    OldRange = (OldMax - OldMin)
    NewRange = (NewMax - NewMin)
    for value in matrix:
        NewMatrix[row_n] = int((((value- OldMin) * NewRange) / OldRange) + NewMin)
        row_n += 1
    return NewMatrix, OldMin, OldMax






def redgreen_contrast_timeseries(records):
    rgb_timeseries = []
    simulation_time = len(records[0][0][0])
    for time in range (simulation_time):
        red_timeseries,green_timeseries = records[-2],records[-1]
        red = red_timeseries[int(J/2),:,time]
        green = green_timeseries[int(J/2),:,time]
        normalised_red = matrix_rgb_normalisation(red)[0]
        normalised_green = matrix_rgb_normalisation(green)[0]
        # normalised_red =preprocessing.normalize([red])*255
        # normalised_green = preprocessing.normalize([green])*255
        zeros = np.zeros(normalised_red.shape)

        # rgb = np.dstack((red,green,zeros))[0]
        rgb = np.dstack((normalised_red,normalised_green,zeros))[0]
    
        # print(normalised_red)

        # rgb = np.rot90(rgb)
        rgb_timeseries.append(rgb)
    return rgb_timeseries


U_record_rgb = np.array(redgreen_contrast_timeseries(U_record))
print(np.shape(U_record_rgb))
# fig= plt.figure(figsize=(10,10))
plt.imshow(U_record_rgb.astype('uint8'), origin='lower', aspect='4')
plt.ylabel('Time', fontsize=15)
plt.xlabel('Space', fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.show()


# %%

# %%
