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

from numerical.cn_plot import plot1D, surfpattern, surfpattern_growth
from database.databaseFunctions import *

#%%


# 1763202_circuit:turinghill_variant:0_samples:2000000


# L=25; dx =0.05; J = int(L/dx)
# T =2000; dt = 0.005; N = int(T/dt)

L=100; dx =0.2; J = int(L/dx)
T =18000; dt = 0.05; N = int(T/dt)


rate=L/T
suggesteddt = float(dx*dx*2)

parID = 2366481#1073797

circuit_n='turinghill'
variant= '11'
n_samples=1000000
ssID = 0
folder = f'{circuit_n}_variant{variant}'
model_param_dict = {'parID':parID, 'circuit_n':circuit_n,'variant':variant, 'n_samples':n_samples}

#%%
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm

from matplotlib.colors import ListedColormap
import seaborn as sns
my_cmap = ListedColormap(sns.color_palette("Spectral",256))   
# from sklearn import preprocessing
import seaborn as sns
def surfpattern(results,L,dx,J,T, record_every_x_hours=10,growth='linear', rate=0, morphogen = 0,savefig=False,filename='',savefigpath='',logResults=False, normalize=False, cmap=my_cmap, space_crop=None):
    

    dx = float(L)/float(J-1)
    x_grid = np.array([j*dx for j in range(J)])
    t_grid = np.arange(0,T,10) 
    
    if normalize == True:
        print('NEEDS NORMALIZATION')
    results = results[morphogen]

    values = results.reshape(len(x_grid),len(t_grid))
    x, t = np.meshgrid(x_grid, t_grid)

    # t,x = np.meshgrid(t_grid, x_grid)
    # plt.contourf(t,x,results, cmap=cmap)
    plt.contourf(x,t,results, levels=100, cmap=cmap)
    if logResults==True:
        plt.colorbar(label='Concentration (logscale)')
    else:
        plt.colorbar()


    plt.ylabel('Time')
    plt.xlabel('Space')
    if savefig==True:
        plt.savefig('%s%s.pdf'%(savefigpath,filename))
        plt.savefig('%s%s.png'%(savefigpath,filename), dpi=1000)
        plt.show()
        plt.close()

    else:
        plt.show()





def surfpattern_growth(results,L,dx,J,T, masking=False, record_every_x_hours=10,growth='linear', rate=0, morphogen = 0,savefig=False,filename='',savefigpath='',logResults=False, normalize=False, cmap=my_cmap, space_crop=None):
    
    def create_growth_mask(shape):
        height, width = shape
        middle_point = width // 2
        mask = np.zeros((height, width))

        for t in range(height):
            growth_extent = int(1 + t*width/2/height)
            start = max(middle_point - growth_extent, 0)
            end = min(middle_point + growth_extent, width)
            mask[t, start:end] = 1
            
        return mask

    mask = create_growth_mask(np.shape(results[0]))




    dx = float(L)/float(J-1)
    x_grid = np.array([j*dx for j in range(J)])
    t_grid = np.arange(0,T,10) 
    
    if normalize == True:
        print('NEEDS NORMALIZATION')
    results = results[morphogen]
    if masking == True:
        results = results * mask
        # Create a masked array where zeros are masked
        results = np.ma.masked_where(results == 0, results)


    values = results.reshape(len(x_grid),len(t_grid))
    x, t = np.meshgrid(x_grid, t_grid)

    # t,x = np.meshgrid(t_grid, x_grid)
    # plt.contourf(t,x,results, cmap=cmap)]

    
    print(np.shape(x), np.shape(t), np.shape(results))
    plt.contourf(x,t,results, levels=100, cmap=my_cmap)
    if logResults==True:
        plt.colorbar(label='Concentration (logscale)')
    else:
        plt.colorbar()


    plt.ylabel('Time')
    plt.xlabel('Space')
    if savefig==True:
        plt.savefig('%s%s.pdf'%(savefigpath,filename))
        plt.savefig('%s%s.png'%(savefigpath,filename), dpi=1000)
        plt.show()
        plt.close()

    else:
        plt.show()
#%%
sns.set_context('poster')

for mechanism, boundaryCoeff in zip(['nogrowth', 'openboundary', 'edgegrowth2'], [1,2,2]):
    filename= lambda mechanism, parID: 'circuit%s_variant%s_bc%s_%s_rate%s_ID%s_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundaryCoeff, mechanism,rate,parID,L,J,T,N)
    simulation_param_dict = {'L':L, 'dx':dx, 'J':J, 'T':T, 'dt':dt, 'N':N, 'boundaryCoeff':boundaryCoeff, 'mechanism':mechanism, 'growth_rate': rate}

    U_final = query_simulationOutput_single_from_sql(simulation_param_dict,model_param_dict,'U_final_1D', ssID=0)
    plot1D(U_final, savefig=True,savefigpath = '/Users/mo2016/Documents/modelling/growth/out/paper_figures/',filename='1D_%s'%filename(mechanism,parID))

    U_record = query_simulationOutput_single_from_sql(simulation_param_dict,model_param_dict,'U_record_1D', ssID=0)
    if mechanism == 'edgegrowth2':
        surfpattern_growth(U_record,L,dx,J,T,masking=True, morphogen=0, savefig=True,savefigpath = '/Users/mo2016/Documents/modelling/growth/out/paper_figures/',filename='2D_morphogen0_%s'%filename(mechanism,parID), cmap='coolwarm')
        surfpattern_growth(U_record,L,dx,J,T,masking=True, morphogen=1, savefig=True,savefigpath = '/Users/mo2016/Documents/modelling/growth/out/paper_figures/',filename='2D_morphogen1_%s'%filename(mechanism,parID), cmap='coolwarm')

    else:
        surfpattern(U_record,L,dx,J,T, morphogen=0,savefig=True,savefigpath = '/Users/mo2016/Documents/modelling/growth/out/paper_figures/',filename='2D_morphogen0_%s'%filename(mechanism,parID), cmap='coolwarm')
        surfpattern(U_record,L,dx,J,T, morphogen=1,savefig=True,savefigpath = '/Users/mo2016/Documents/modelling/growth/out/paper_figures/',filename='2D_morphogen1_%s'%filename(mechanism,parID), cmap='coolwarm')
    plt.show()
    plt.close('all')


#%%
lsa_df = pickle.load(open( modellingpath + f'/growth/out/analytical/lsa_dataframes/lsa_df_turinghill_variant{variant}_{n_samples}parametersets.pkl', "rb" ) )
# lsa_df = pickle.load(open( modellingpath + f'/growth/out/analytical/lsa_dataframes/lsa_df_turinghill_variant0_{n_samples}parametersets.pkl', "rb" ) )
# lsa_df = pickle.load(open( modellingpath + f'/growth/out/analytical/instability/instability_df_circuitturinghill_variant8-9_combinedparametersets.pkl', "rb" ) )
lsa_df.loc[parID]

#%%
#linear stability analysis of parID
from analytical.linear_stability_analysis import big_turing_analysis_df, detailed_turing_analysis_dict
from randomfunctions import *

# par_dict = lsa_df.loc[parID].iloc[1].to_dict() #converts a dataframe row into a dictionary outputing a dictionary for a specific parameter set
par_dict = lsa_df.loc[parID,0].to_dict() #converts a dataframe row into a dictionary outputing a dictionary for a specific parameter set
#Run analysis on 1M parameter sets
out = detailed_turing_analysis_dict(par_dict, circuit_n, 2)
# plot_all_dispersion(out[4][0],2, crop=20)
plot_highest_dispersion_noticks(out[4][0], crop=65)
plt.savefig('/Users/mo2016/Documents/modelling/growth/out/paper_figures/lsa_%s.pdf'%filename(mechanism,parID))

# %%
model_param_id = 'model_param_id:2366481_circuit:turinghill_variant:11_samples:1000000'
par_dict = model_param_dict_from_model_param_id(model_param_id); ssID=0
out = detailed_turing_analysis_dict(par_dict, 'turinghill', 2)
plot_highest_dispersion_noticks(out[4][ssID],crop = 400, top = 2000)
plt.show()
# %%
