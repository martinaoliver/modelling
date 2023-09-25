#%%
#############
###paths#####
#############
import sys
import os

pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')

import copy
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
import pickle
import numba
from numba import cuda, float32

@numba.jit(nopython=True)
def check_neighbours(cell_matrix,y_pos,x_pos): #returns grid with the neighbouring points
    top_array = [cell_matrix[y_pos-1, x_pos-1], cell_matrix[y_pos-1,x_pos], cell_matrix[y_pos-1,x_pos+1]]
    middle_array = [cell_matrix[y_pos, x_pos-1], np.nan, cell_matrix[y_pos,x_pos+1]]
    bottom_array = [cell_matrix[y_pos+1, x_pos-1], cell_matrix[y_pos+1,x_pos], cell_matrix[y_pos+1,x_pos+1]]
    neighbours_cellmatrix = np.array([top_array,middle_array,bottom_array])
    return neighbours_cellmatrix

def cell_automata_colony(cell_matrix, p_division):
    daughterToMotherDict = {}
    cell_matrix_new = copy.deepcopy(cell_matrix)
    for y_pos in np.linspace(1,len(cell_matrix)-2,len(cell_matrix)-2):
        for x_pos in np.linspace(1,len(cell_matrix)-2,len(cell_matrix)-2):
            y_pos = int(y_pos)
            x_pos = int(x_pos)
            if cell_matrix[y_pos, x_pos]!=0:
                neighbours_cellmatrix = check_neighbours(cell_matrix,y_pos,x_pos)

                if 0 in neighbours_cellmatrix:
                    cell_division=np.random.choice([1,0],p=[p_division,1-p_division])
                    if cell_division==1:
                        #choose cell to divide to
                        index_nocells=np.where(np.array(neighbours_cellmatrix )== 0) #find empty cells
                        divided_cell_index = np.random.choice(range(len(index_nocells[0]))) #select empty cell to divide to
                        index_newcell_y, index_newcell_x = (index_nocells[n][divided_cell_index] for n in range(2))
                        
                        #create mask
                        cell_matrix_new[index_newcell_y+y_pos-1,index_newcell_x+x_pos-1]=1

                        #create cellular tissue memory
                        daughterToMotherDict[(index_newcell_y+y_pos-1,index_newcell_x+x_pos-1)] = (y_pos,x_pos)
    return cell_matrix_new, daughterToMotherDict
    
# numba.jit(nopython=True)
def adi_ca(L,dx,J,T,dt,N,n_species,divisionTimeHours,tqdm_disable=False, p_division=0.3,stochasticity=0, seed=1, growth='Fast'):
    I=J
    print(f'dt{dt}')


    #Define initial conditions and cell matrix
    U0 = []
    perturbation=0.001
    steadystates=[0.1]*n_species
    np.random.seed(seed)

    cell_matrix = np.zeros(shape=(I,J))
    cell_matrix[int(I/2), int(J/2)] = 1
    # for index in range(n_species):
    #     U0.append(np.ones((I, J)))
    # U0 = U0*cell_matrix

    # U = copy.deepcopy(U0)

    divisionTimeUnits=divisionTimeHours/dt
    print(f'divisionTimeUnits{divisionTimeHours}')
    # cell_matrix_record = np.zeros([J, I, int(T/divisionTimeHours)])
    cell_matrix_record = np.zeros([J, I, N])
    memory_matrix_record = np.zeros([J, I, N], dtype='i,i')
    cell_matrix_record[:, :, 0] = cell_matrix #issue in this line
    memory_matrix_record[int(I/2), int(J/2), :] = int(I/2), int(J/2)#issue in this line
    memory_matrix = memory_matrix_record[:,:,0]
    # print('Initial Mask')
    # plt.imshow(cell_matrix)
    # plt.show()
    divide_counter=0
    daughterToMotherDictList = []

    for ti in tqdm(range(1,N), disable = tqdm_disable):
        # print(ti)
        daughterToMotherDict = {}


        # U_new = copy.deepcopy(U)

        hour = ti / (N / T)

        if (ti%divisionTimeUnits==0):
            cell_matrix_new ,daughterToMotherDict = cell_automata_colony( cell_matrix, p_division)
            cell_matrix = copy.deepcopy(cell_matrix_new)
            
            # cell_matrix_record[:, :, divide_counter] = cell_matrix #issue in this line
            divide_counter+=1
        else:
            memory_matrix = np.zeros([J, I], dtype='i,i')
        
        cell_matrix_record[:, :, ti] = cell_matrix #issue in this line
        memory_matrix_record[:, :, ti] = memory_matrix #issue in this line
        daughterToMotherDictList.append(daughterToMotherDict)
    # U = copy.deepcopy(U_new)
    # print(np.shape(cell_matrix_record))
    return cell_matrix_record,memory_matrix_record, daughterToMotherDictList

def maskFunction(L=9, dx=0.05, T=50, dt=0.05, divisionTimeHours=1, p_division=0.5,seed=1,plot1D=False, plotScatter = False, plotVideo=False):
        
    #execution parameters
    n_species=6

    J = int(L/dx)
    N = int(T/dt)

    suggesteddt = float(dx*dx*2)
    print(f'suggested dt = {suggesteddt}')
    # p_division=float(sys.argv[5]);seed=1

    cell_matrix_record,memory_matrix_record, daughterToMotherDictList = adi_ca(L,dx,J,T,dt,N,n_species,divisionTimeHours,tqdm_disable=False,p_division=p_division,seed=seed)
    print(np.shape(cell_matrix_record))
    print(np.shape(cell_matrix_record))
    pickle.dump( cell_matrix_record, open(modellingpath + "/3954/paper/out/numerical/masks/caMask_seed%s_pdivision%s_L%s_J%s_T%s_N%s.pkl"%(seed,p_division,L,J,T,N), "wb" ) )
    pickle.dump( memory_matrix_record, open(modellingpath + "/3954/paper/out/numerical/masks/caMemory_seed%s_pdivision%s_L%s_J%s_T%s_N%s.pkl"%(seed,p_division,L,J,T,N), "wb" ) )
    pickle.dump( daughterToMotherDictList, open(modellingpath + "/3954/paper/out/numerical/masks/caMemory_seed%s_pdivision%s_L%s_J%s_T%s_N%s.pkl"%(seed,p_division,L,J,T,N), "wb" ) )
    # print(np.shape(cell_matrix_record))
    # for ti in range(N):
        # print(np.shape(memory_matrix_record[:,:,ti]))
        # print(memory_matrix_record[:,:,ti])
    #%%
    if plot1D == True:
        print('Final Mask')
        plt.imshow(cell_matrix_record[:,:,-1], cmap='Greys')# plot_2D_final_concentration(final_concentration,grids,filename,n_species=n_species)
        tick_positions = np.linspace(0, L/dx, 4)
        tick_labels = np.linspace(0, L, 4).round(decimals=2)
        plt.xticks(tick_positions, tick_labels)
        plt.yticks(tick_positions, tick_labels)
        plt.savefig(modellingpath + "/3954/paper/out/numerical/masks/caMask_seed%s_pdivision%s_L%s_J%s_T%s_N%s.png"%(seed,p_division,L,J,T,N))
        plt.savefig(modellingpath + "/3954/paper/out/numerical/masks/caMask_seed%s_pdivision%s_L%s_J%s_T%s_N%s.pdf"%(seed,p_division,L,J,T,N))
        plt.show()
        plt.close()
    #%%
    if plotVideo==True:
        def show_rgbvideo(timeseries_unstacked):
            time=0

            fig = plt.plot()
            rgb_timeseries=timeseries_unstacked # Read the numpy matrix with images in the rows
            N = len(rgb_timeseries)
            im=plt.imshow(rgb_timeseries[:,:,0].astype('uint8'), origin= 'lower', cmap='Greys')
            for time in range(len(rgb_timeseries[0,0,:])):
                im.set_data(rgb_timeseries[:,:,time].astype('uint8'))
                
                plt.xlabel(time)
                plt.pause(0.00000000000000001)
            plt.show()

        show_rgbvideo(cell_matrix_record)

    #%%
    if plotScatter==True:
        print('Scatter')
        lenght_list = []
        for i in range(len(cell_matrix_record[0,0,:])):
            lenght = np.count_nonzero(cell_matrix_record[:,int(J/2),i])
            lenght_list.append(lenght)
        lenght_list = np.array(lenght_list)*dx

        plt.scatter(range(0,N),lenght_list, c='k',s=1)
        plt.xlabel('Time (hours)',fontsize=15)
        plt.ylabel('Colony diameter (mm)', fontsize=15)
        tick_positions = np.linspace(0, T/dt, 4)
        tick_labels = np.linspace(0, T , 4).round(decimals=2)
        plt.xticks(tick_positions, tick_labels,fontsize=15)
        plt.savefig(modellingpath + "/3954/paper/out/numerical/masks/growthScatter_seed%s_pdivision%s_L%s_J%s_T%s_N%s.png"%(seed,p_division,L,J,T,N))
        plt.savefig(modellingpath + "/3954/paper/out/numerical/masks/growthScatter_seed%s_pdivision%s_L%s_J%s_T%s_N%s.pdf"%(seed,p_division,L,J,T,N))
        
        print('y')
        plt.show()

    return cell_matrix_record, memory_matrix_record, daughterToMotherDictList

    # %%
    if plotVolume == True:
        plt.rcParams["figure.figsize"] = [7.00, 3.50]
        plt.rcParams["figure.autolayout"] = True
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        z, x, y = cell_matrix_record.nonzero()
        my_cmap = plt.get_cmap('magma')

        ax.plot_trisurf(z,x,y, linewidth = 100,
                        antialiased = True,alpha=1,cmap=my_cmap)
        ax.set_xlabel('y', labelpad=-10,size=20)
        ax.set_ylabel('x', labelpad=-10,size=20)
        ax.set_zlabel('Time', labelpad=-10,size=20, rotation=90)
        # ax.set_ylabel('Y', rotation=0, labelpad=10)
        color_tuple = (1.0, 1.0, 1.0, 0.0)

        ax.tick_params(color=color_tuple, labelcolor=color_tuple)
        ax.grid(False)
        ax.view_init(20,40, 0)

        plt.show()


    return cell_matrix_record, memory_matrix_record, daughterToMotherDictList



#%%

#slowgrowth
L=20; dx =0.1; J = int(L/dx)
T =100; dt = 0.02; N = int(T/dt)
boundaryCoeff = 1
division_time_hours=0.5
p_division=0.38;seed=1
# 
# #mediumgrowth
# L=20; dx =0.1; J = int(L/dx)
# T =50; dt = 0.02; N = int(T/dt)
# boundaryCoeff = 1
# division_time_hours=0.5
# p_division=1;seed=1

# # fastgrowth
L=20; dx =0.1; J = int(L/dx)
T =25; dt = 0.02; N = int(T/dt)
boundaryCoeff = 1
division_time_hours=0.2
p_division=0.7;seed=1

dt = 0.2; N = int(T/dt)


# cell_matrix_record, memory_matrix_record, daughterToMotherDictList = maskFunction(L=L,dx=dx, T=T ,dt=dt, divisionTimeHours=division_time_hours, p_division=p_division, plot1D=True, plotScatter=True)
# # %%


# %%
