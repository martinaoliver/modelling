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
    import matplotlib as mpl
    mpl.use('tkagg')

modulepath = modelling_local + '/3954/modules/new_CN'

sys.path.append(modulepath)


from datetime import timedelta
import time
import scipy.io
import numpy as np
import pickle
from numpy import linalg as LA
from parametercombination_analysis import *
from randomfunctions import *
import datetime
from numerical_solvers_variableboundary import *
from dispersionrelation_functions import *
import sys
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import cm
import decimal
from scipy import signal
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.ticker as ticker
from math import isclose
path =  sixeqpath + '/numerical_confocal'
PCApath = sixeqpath + '\PCA'
parametersearchpath = sixeqpath + '/parameter_space_search'


def plot_1D_final_concentration(final_concentration,grids,par_ID, L,J,T,N,boundary_coef,mechanism = 'turing', circuit_n =2):
    A = final_concentration[0]
    B = final_concentration[1]
    C = final_concentration[2]
    D = final_concentration[3]
    E = final_concentration[4]
    F = final_concentration[5]

    x_grid = grids[0]


    fig, (ax1,ax2,ax3,ax4,ax5,ax6) = plt.subplots(6, figsize=(10,10))
    fig.suptitle('Crank Nicolson 1D Numerical solution  par_ID%r'%par_ID)
    ax1.plot(x_grid, A, 'r', label = 'A')
    ax1.set(title = 'A', xlabel = 'x (space units)', ylabel = '[A]', ylim=(0, round(np.amax(A), 3) + 0.15*np.amax(A)))
    ax2.plot(x_grid, B,'g', label = 'B')
    ax2.set(title = 'B', xlabel = 'x (space units)', ylabel = '[B]', ylim=(0, round(np.amax(B), 3) + 0.15*np.amax(B)))
    ax3.plot(x_grid, C, 'b',label = 'C')
    ax3.set(title = 'C', xlabel = 'x (space units)', ylabel = '[C]', ylim=(0, round(np.amax(C), 3) + 0.15*np.amax(C)))
    ax4.plot(x_grid, D, 'c', label = 'D')
    ax4.set(title = 'D', xlabel = 'x (space units)', ylabel = '[D]', ylim=(0, round(np.amax(D), 3) + 0.15*np.amax(D)))
    ax5.plot(x_grid, E,'r', label = 'E')
    ax5.set(title = 'E', xlabel = 'x (space units)', ylabel = '[E]', ylim=(0, round(np.amax(E), 3) + 0.15*np.amax(E)))
    ax6.plot(x_grid, F, 'y',label = 'F')
    ax6.set(title = 'F', xlabel = 'x (space units)', ylabel = '[F]', ylim=(0, round(np.amax(F), 3) + 0.15*np.amax(F)))

    fig.tight_layout()
    fig.subplots_adjust(top = 0.85)

    filename = 'CN_circuit%rboundary%r_%sID%r_L%r_J%r_T%r_N%r' % (circuit_n,boundary_coef, mechanism, par_ID, L, J, T, N)
    plt.savefig(path + '/results/figures/1D/nonturing_numerical/growing_colony/%s.jpeg' % filename)
    plt.close()

def plot_2D_final_concentration(final_concentration,grids,filename):
    A = final_concentration[0]
    B = final_concentration[1]
    C = final_concentration[2]
    D = final_concentration[3]
    E = final_concentration[4]
    F = final_concentration[5]

    x_grid = grids[0]
    y_grid = grids[1]


    fig, axs = plt.subplots(2,3,figsize=(7.5,4))

    black_yellow= [(45/255,45/255,45/255),(255/255,255/255,0/255)]
    black_red = [(45/255,45/255,45/255),(255/255,0/255,0/255)]
    black_green = [(45/255,45/255,45/255),(50/255,205/255,50/255)]
    yellow_cmap = LinearSegmentedColormap.from_list('black_yellow',black_yellow)
    red_cmap = LinearSegmentedColormap.from_list('black_red',black_red)
    green_cmap = LinearSegmentedColormap.from_list('black_green',black_green)


    i=0
    im1 = axs[0,0].pcolormesh(x_grid, y_grid, A,vmin = 0, vmax = round(np.amax(A), 3), label = 'A', shading = 'flat',cmap=yellow_cmap)
    im2 = axs[0,1].pcolormesh(x_grid, y_grid, B,vmin = 0, vmax = round(np.amax(B), 3), label = 'B', shading = 'flat',cmap=yellow_cmap)
    im3 = axs[0,2].pcolormesh(x_grid, y_grid, C,vmin = 0, vmax = round(np.amax(C), 3), label = 'C', shading = 'flat',cmap=yellow_cmap)
    im4 = axs[1,0].pcolormesh(x_grid, y_grid, D,vmin = 0, vmax = round(np.amax(D), 3), label = 'D', shading = 'flat',cmap=green_cmap)
    im5 = axs[1,1].pcolormesh(x_grid, y_grid, E,vmin = 0, vmax = round(np.amax(E), 3), label = 'E', shading = 'flat',cmap=red_cmap)
    im6 = axs[1,2].pcolormesh(x_grid, y_grid, F,vmin = 0, vmax = round(np.amax(F), 3), label = 'F', shading = 'flat',cmap=green_cmap)





    for ax in axs.flat:
        ax.label_outer()

    count1=0
    morphogens = ('A','B','C','D','E','F')
    ims = (im1,im2,im3,im4,im5,im6)
    for ax in axs.flat:
        ax.set(title=morphogens[count1])
        fig.colorbar(ims[count1], ax=ax)

        count1+=1

    plt.tight_layout()
    plt.savefig(path + '/results/figures/2D/nonturing_numerical/robustness4/growing_colony/%s.jpeg' % filename)
    plt.close()

dimension = str(sys.argv[1])
df_n = int(sys.argv[5])
# dimension = '2D'

# open parameter dictionaries
# turing_df= pickle.load( open(parametersearchpath +  '/results/turing_dataframes/circuit2robustness4_turingI_monostable_df.pkl', "rb" ) )
nonturing_df = pickle.load(open(parametersearchpath + '/parameterfiles/df_parameterfile_circuit2_robustness4_2020-12-09/df_parameterfile_2020-12-09_ID%r.pkl'%df_n, "rb"))
#turing_df= pickle.load( open(parametersearchpath +  '/results/turing_dataframes/circuit2_turingI_monostable_df.pkl', "rb" ) )

print(len(nonturing_df))
#for par_ID in range(368, len(turing_df)):
#interesting_list = [33, 158, 296, 332, 376, 378, 441, 457, 607,693, 690 ,755, 782, 866]
interesting_list = [int(sys.argv[2])]
for par_ID in interesting_list:
    start = time.time()
    print('parID = ' + str(par_ID))
    circuit_n = 2
    mechanism = 'nonturing'
    boundary_coef = 1
    n_species=6

    par_dict = nonturing_df.iloc[par_ID].to_dict()
    lsa_results = detailed_turing_analysis_dict(par_dict, circuit_n,n_species)
    pickle.dump( lsa_results, open( "lsa_results%r.pkl"%par_ID, "wb" ) )
    if dimension == '1D':
        # Define 1D numerical parameters
        L = 8 # total leght of space simulated
        J = L *int(sys.argv[3])  # number of equally spaced gridpoints in space domain (larger J means more spatial precision(tends towards continuum solution) )
        T = int(sys.argv[4])  # total lenght of time simulated
        t_gridpoints = t_gridpoints_stability(L, J, T)  # number of equally spaced gridpoints in domain (larger N means more temporal precision (tends towards continuum solution) )
        N = T * t_gridpoints
        print('N = ' + str(T) + '*' + str(t_gridpoints) + '=' + str(N))
        initial_condition = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001]

        # Run 1D simulation
        records, final_concentration, grids = crank_nicolson(par_dict, initial_condition, L, J, T, N, circuit_n,boundary_coef=boundary_coef)
        plot_1D_final_concentration(final_concentration, grids,par_ID, L, J, T, N)

    if dimension == '2D':
     # Define 2D numerical parameters
        L_x = 8# total leght of space x simulated
        L_y = 8# total leght of space y simulated
        x_gridpoints = int(sys.argv[3])
        J = L_x * x_gridpoints  # number of equally spaced gridpoints in space x domain (larger J means more spatial precision(tends towards continuum solution) )
        I = L_y * x_gridpoints  # number of equally spaced gridpoints in space y domain (larger J means more spatial precision(tends towards continuum solution) )

        T =int(sys.argv[4])# 400#total lenght of time simulated
        t_gridpoints = t_gridpoints_stability(L_x, J,
                                              T)  # number of equally spaced gridpoints in domain (larger N means more temporal precision (tends towards continuum solution) )
        N = T * t_gridpoints

        print('N = ' + str(T) + '*' + str(x_gridpoints) + '=' + str(N))
        stability_test(L_x, J, T, N)  # if stability test fails, oscillations might occur
        shape_name = 'growing_colony'
        filename = 'ADI_circuit%rboundary%r_%s_%sID%r_df%r_L%r_J%r_T%r_N%r'%(circuit_n,boundary_coef, shape_name, mechanism,par_ID,df_n,L_x,J,T,N)
        initial_condition = [0.001,0.001,0.001,0.001,0.001,0.001]
        # Run 2D simulation
        records,final_concentration,grids= adi_shape(par_dict,initial_condition,L_x,L_y,J,I,T,N, circuit_n,shape_name, boundary_coef=boundary_coef)
        end = time.time()
        seconds = end - start
        execution_time = timedelta(seconds=seconds)
        print('Execution time =' + str(execution_time))
        # records,final_concentration,grids = adi_growth(par_dict,initial_condition,L_x,L_y,J,I,T,N, circuit_n, boundary_coef=boundary_coef)
        # records,final_concentration,grids = adi(par_dict,initial_condition,L_x,L_y,J,I,T,N, circuit_n, boundary_coef=0)
        #plot_2D_final_concentration(final_concentration,grids,filename)
        plot_redgreen_contrast(final_concentration,x_gridpoints)
        plt.ylabel('y axis (mm)', size=16)
        plt.xlabel('x axis (mm)', size=16)
        plt.yticks(size=15)
        plt.xticks(size=15)
        plt.tight_layout()
        plt.savefig(path + '/results/figures/redgreen/nonturing_numerical/robustness4/growing_colony/8x%r_T%r/redgreen_%s.png' % (x_gridpoints,T,filename))
        plt.close()
        ephemeral_path = '/rds/general/user/mo2016/ephemeral/Documents/modelling/6eq/numerical_confocal'
        #save_numerical_results(records, grids, filename,ephemeral_path)

        rgb_timeseries = redgreen_contrast_timeseries(records)
        pickle_out = open(ephemeral_path + '/results/simulation/redgreen/nonturing/8x%r_T%r/redgreen_%s.pkl' % (x_gridpoints,T,filename),'wb')
        pickle.dump(rgb_timeseries, pickle_out)
        pickle_out.close()
