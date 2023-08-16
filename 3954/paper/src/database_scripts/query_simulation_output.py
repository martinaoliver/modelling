
#%%
import sys
import os
pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')
from database.databaseFunctions import *
import pickle
#############
#############
###paths#####
#############

import sys
import os
import pickle
import psycopg2
import matplotlib.pyplot as plt
from matplotlib import animation
from numerical.plotting_numerical import *
import seaborn as sns 


# %%

# slow
L=20; dx =0.1; J = int(L/dx)
T =100; dt = 0.02; N = int(T/dt)
boundaryCoeff = 1
division_time_hours=0.5
p_division=0.38;seed=1


# # medium
# # medium
# L=20; dx =0.1; J = int(L/dx)
# T =50; dt = 0.02; N = int(T/dt)
# T =50; dt = 0.02; N = int(T/dt)
# boundaryCoeff = 1
# division_time_hours=0.5
# p_division=1;seed=1





# fast
L=20; dx =0.1; J = int(L/dx)
T =25; dt = 0.02; N = int(T/dt)
boundaryCoeff = 1
division_time_hours=0.2
p_division=0.7;seed=1
x_gridpoints=int(1/dx)


shape='ca'
simulation_param_dict = {'L':L, 'dx':dx, 'J':J, 'T':T, 'dt':dt, 'N':N, 
            'boundaryCoeff':boundaryCoeff, 
            'shape':'ca', 'p_division': p_division, 'seed':seed, 'division_time_hours':division_time_hours}




ssID=0
circuit_n=14
#  circuit_n='circuit14'
variant='2nd' #variant='fitted7'#f'fitted7_gaussian4187715_nsr{0.01}'#
balance='Balanced'
Kce=100
n_samples = 1000000 #n_samples = 13700000
folder = 'circuit14variant2ndBalancedKce100'
# folder = f'circuit14variantfitted7_gaussian4187715'

# model_param_dict = {'parID':1, 'circuit_n':circuit_n,'variant':variant, 'n_samples':n_samples, 'balance':balance}
parID=141318

model_param_dict = {'parID':parID, 'circuit_n':circuit_n,'variant':variant, 'n_samples':n_samples}

filename= lambda parID: 'circuit%r_variant%s_bc%s_%s_ID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundaryCoeff, shape,parID,L,J,T,N)




#%%
simulationOutput = query_simulationOutput_single_from_sql(simulation_param_dict, model_param_dict,query_column = 'U_final_1D', ssID=ssID)
saveFigurePath = modellingpath + '/3954/paper/out/numerical/colonies/figures/%s/'%folder

plot_redgreen_contrast(simulationOutput,L, filename=filename(parID), save_figure=True, path=saveFigurePath)





#%%
U_record = query_simulationOutput_single_from_sql(simulation_param_dict, model_param_dict,query_column = 'U_record_1D', ssID=ssID)

#%%
def kymograph_U_record(U_record, J, folder, filename, parID):
    U_record_rgb = np.array(redgreen_contrast_timeseries(U_record))

    # print(np.shape(U_record_rgb))
    crop_U_record_rgb = U_record_rgb[:,int(J/2)]
    # plt.figure(figsize=(10, 2))
    plt.imshow(crop_U_record_rgb.astype('uint8'), origin='lower',aspect=2)
    plt.ylabel('Time', fontsize=15)
    plt.xlabel('Space', fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    filename= lambda parID: 'circuit%s_variant%s_%sparametersets_balanced_Kce%s_bc%s_%s_ID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant,n_samples,Kce,boundaryCoeff, shape,parID,L,J,T,N)

    plt.savefig(modellingpath + f'/3954/paper/out/numerical/colonies/figures/{folder}/kymograph_{filename(parID)}.pdf',dpi=2000)
    plt.show()

kymograph_U_record(U_record, J, folder, filename, parID)
# %%
# %%

rgb_timeseries = redgreen_contrast_timeseries(U_record)
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
    # ani.save('/%s.mp4' %filename,writer=mywriter)
    print('Video saved', filename)

save_rgbvideo(rgb_timeseries, saveVideoPath, filename(parID))
print('Video saved', filename(parID))

# %%
saveFigPath = modellingpath + '/3954/paper/out/numerical/colonies/figures/%s/'%folder
show_red_green_separately(U_final, savePath=saveFigPath, filename=filename(parID), saveFig=True)
# %%
plot_redgreen_contrast(U_final,L,parID=6 ,save_figure=True, path=saveFigPath, filename=filename(parID))

# def plot_redgreen_contrast(final_concentration, mm,filename=None, path=None, parID=0, scale_factor=10, save_figure=False):

# %%
time=20
plot_redgreen_contrast(U_record[:][:,:,:,time],L,parID=parID ,save_figure=True, path=saveFigurePath, filename=filename(parID)+('_time%r'%time))
# time=10
# plot_redgreen_contrast(U_record[:][:,:,:,time],L,parID=parID ,save_figure=True, path=saveFigurePath, filename=filename(parID)+('_time%r'%time))

# %%
