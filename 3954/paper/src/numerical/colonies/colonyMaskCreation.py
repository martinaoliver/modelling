#############
###paths#####
#############
import sys
import os

from importlib_metadata import distribution
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
    memory_matrix_new =  np.zeros((len(cell_matrix),len(cell_matrix)),dtype='i,i')
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
def adi_ca(initial_condition,L,dx,J,T,dt,N,n_species,tqdm_disable=False, p_division=0.3,stochasticity=0, seed=1, growth='Fast'):
    L_x=L;L_y=L;dy=dx;I=J

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

    divisionTimeHours=1
    divisionTimeUnits=int(divisionTimeHours/dt)
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

#execution parameters

save_figure = False
tqdm_disable = False #disable tqdm
n_species=6

# L=5; x_gridpoints =5; J = L*x_gridpoints
# T =10; t_gridpoints = 10; N = T*t_gridpoints
# L=8; dx =0.05; J = int(L/dx)
# T =125; dt = 0.05; N = int(T/dt)
# T =125; dt = 0.04; N = int(T/dt)

L=8; dx =0.5; J = int(L/dx)
T =125; dt = 0.5; N = int(T/dt)

# L=int(sys.argv[1]); x_gridpoints =int(sys.argv[2]); J = L*x_gridpoints
# T =int(sys.argv[3]); t_gridpoints = int(sys.argv[4]); N = T*t_gridpoints
L_x = L
L_y = L
I = J

p_division=0.05;seed=1
# p_division=float(sys.argv[5]);seed=1
initial_condition = [1000]*n_species
cell_matrix_record,memory_matrix_record, daughterToMotherDictList = adi_ca(initial_condition,L,dx,J,T,dt,N,n_species,tqdm_disable=False,p_division=p_division,seed=seed)
pickle.dump( cell_matrix_record, open(modellingpath + "/3954/paper/out/numerical/masks/caMask_seed%s_pdivision%s_L%s_J%s_T%s_N%s.pkl"%(seed,p_division,L,J,T,N), "wb" ) )
pickle.dump( memory_matrix_record, open(modellingpath + "/3954/paper/out/numerical/masks/caMemory_seed%s_pdivision%s_L%s_J%s_T%s_N%s.pkl"%(seed,p_division,L,J,T,N), "wb" ) )
pickle.dump( daughterToMotherDictList, open(modellingpath + "/3954/paper/out/numerical/masks/caMemory_seed%s_pdivision%s_L%s_J%s_T%s_N%s.pkl"%(seed,p_division,L,J,T,N), "wb" ) )

# print(np.shape(cell_matrix_record))
# for ti in range(N):
    # print(np.shape(memory_matrix_record[:,:,ti]))
    # print(memory_matrix_record[:,:,ti])

plot1D=True
if plot1D == True:
    print('Final Mask')
    plt.imshow(cell_matrix_record[:,:,-1], cmap='Greys')# plot_2D_final_concentration(final_concentration,grids,filename,n_species=n_species)
    tick_positions = np.linspace(0, L/dx, 4)
    tick_labels = np.linspace(0, L, 4).round(decimals=2)
    plt.xticks(tick_positions, tick_labels)
    plt.yticks(tick_positions, tick_labels)
    plt.savefig(modellingpath + "/3954/paper/out/numerical/masks/caMask_seed%s_pdivision%s_L%s_J%s_T%s_N%s.png"%(seed,p_division,L,J,T,N))
    plt.show()
    plt.close()

plotVideo=False
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
            plt.pause(0.0001)
        plt.show()

    show_rgbvideo(cell_matrix_record)


plotScatter = True
if plotScatter==True:
    print('Scatter')
    lenght_list = []
    for i in range(len(cell_matrix_record[0,0,:])):
        lenght = np.count_nonzero(cell_matrix_record[:,int(J/2),i])
        lenght_list.append(lenght)
    lenght_list = np.array(lenght_list)*dx

    plt.scatter(range(0,N),lenght_list, c='k',s=1)
    plt.xlabel('Time (hours)')
    plt.ylabel('Colony diameter (mm)')
    tick_positions = np.linspace(0, T/dt, 4)
    tick_labels = np.linspace(0, T , 4).round(decimals=2)
    plt.xticks(tick_positions, tick_labels)
    plt.savefig(modellingpath + "/3954/paper/out/numerical/masks/growthScatter_seed%s_pdivision%s_L%s_J%s_T%s_N%s.png"%(seed,p_division,L,J,T,N))
    plt.show()

