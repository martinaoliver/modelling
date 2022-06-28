import sys
import os
import numpy as np
import copy
import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
from tqdm import tqdm
import pickle


pwd = os.getcwd()
root = pwd.rpartition("mo2016")[0] + pwd.rpartition("mo2016")[1] #/Volumes/mo2016/ or '/Users/mo2016/' or '/rds/general/mo2016/'
if root == '/Users/mo2016':
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

# def adi_ca(par_dict,L_x,L_y,J,I,T,N, circuit_n, n_species,D,tqdm_disable=False, p_division=0.5,stochasticity=0):
def adi_ca(initial_condition,L_x,L_y,J,I,T,N, n_species,tqdm_disable=False, p_division=0.5,stochasticity=0, seed=1, growth='Fast'):


    #spatial variables
    dx = float(L_x)/float(J-1); dy = float(L_y)/float(I-1)
    x_grid = np.array([j*dx for j in range(J)]); y_grid = np.array([i*dy for i in range(I)])


    dt = float(T)/float(N-1)
    t_grid = np.array([n*dt for n in range(N)])


    #Define initial conditions and cell matrix
    U0 = []
    perturbation=0.001
    steadystates=[0.1]*n_species
    np.random.seed(seed)

    cell_matrix = np.zeros(shape=(I,J))
    cell_matrix[int(I/2), int(J/2)] = 1
    for index in range(n_species):
        U0.append(np.random.uniform(low=steadystates[index] - perturbation, high=steadystates[index] + perturbation, size=(I, J)))
    U0 = U0*cell_matrix


    def check_neighbours(cell_matrix,y_pos,x_pos): #returns grid with the neighbouring points
        top_array = [cell_matrix[y_pos-1, x_pos-1], cell_matrix[y_pos-1,x_pos], cell_matrix[y_pos-1,x_pos+1]]
        middle_array = [cell_matrix[y_pos, x_pos-1], np.nan, cell_matrix[y_pos,x_pos+1]]
        bottom_array = [cell_matrix[y_pos+1, x_pos-1], cell_matrix[y_pos+1,x_pos], cell_matrix[y_pos+1,x_pos+1]]
        neighbours_cellmatrix = np.array([top_array,middle_array,bottom_array])
        return neighbours_cellmatrix
    def cell_automata_colony(species_list,cell_matrix, p_division):
        new_species_list = copy.deepcopy(species_list)
        original_species_list = copy.deepcopy(species_list)
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
                            index_nocells=np.where(np.array(neighbours_cellmatrix )== 0)
                            divided_cell_index = np.random.choice(range(len(index_nocells[0])))
                            index_newcell_y, index_newcell_x = (index_nocells[n][divided_cell_index] for n in range(2))
                            for count,species in enumerate(original_species_list):
                                new_species_list[count][index_newcell_y+y_pos-1,index_newcell_x+x_pos-1] += species[y_pos,x_pos]/2
                                new_species_list[count][y_pos,x_pos] += species[y_pos,x_pos]/2
                            cell_matrix_new[index_newcell_y+y_pos-1,index_newcell_x+x_pos-1]=1


        return new_species_list, cell_matrix_new


    U = copy.deepcopy(U0)
    cell_matrix_record = np.zeros([J, I, T])


    unittime=0
    for ti in tqdm(range(N), disable = tqdm_disable):

        U_new = copy.deepcopy(U)

        hour = ti / (N / T)
        if hour % 1 == 0:  #only consider division at unit time (hour)
            if growth=='Slow':
                #predict if division occurs based on the p_division, the current cell matrix
                #return new cell matrix and updated concentrations with dilution
                U_new, cell_matrix_new = cell_automata_colony(U_new, cell_matrix, p_division)
                cell_matrix = copy.deepcopy(cell_matrix_new)
                cell_matrix_record[:, :, int(hour)] = cell_matrix #issue in this line

        if growth=='Fast':
            #predict if division occurs based on the p_division, the current cell matrix
            #return new cell matrix and updated concentrations with dilution
            U_new, cell_matrix_new = cell_automata_colony(U_new, cell_matrix, p_division)
            cell_matrix = copy.deepcopy(cell_matrix_new)
            cell_matrix_record[:, :, int(hour)] = cell_matrix #issue in this line
        U = copy.deepcopy(U_new)
    print(np.shape(cell_matrix_record))
    return cell_matrix_record

import pickle
#execution parameters

save_figure = False
tqdm_disable = False #disable tqdm
n_species=6

L=5; x_gridpoints =10; J = L*x_gridpoints
T =10; t_gridpoints = 100; N = T*t_gridpoints

L=int(sys.argv[1]); x_gridpoints =int(sys.argv[2]); J = L*x_gridpoints
T =int(sys.argv[3]); t_gridpoints = int(sys.argv[4]); N = T*t_gridpoints
L_x = L
L_y = L
I = J
p_division=float(sys.argv[5]);seed=1
initial_condition = [1000]*n_species

cell_matrix_record= adi_ca(initial_condition,L_x,L_y,J,I,T,N, n_species,tqdm_disable=tqdm_disable,p_division=p_division,seed=seed)
plt.imshow(cell_matrix_record[:,:,-1], cmap='Greys')# plot_2D_final_concentration(final_concentration,grids,filename,n_species=n_species)
plt.savefig('masks/caMask_seed%s_pdivision%s_L%s_J%s_T%s_N%s.png'%(seed,p_division,L,J,T,N))# plot_redgreen_contrast(final_concentration,L,mechanism,shape,filename,modelling_ephemeral,parID=parID,dimension=dimension,scale_factor=x_gridpoints,save_figure=save_figure)
pickle.dump( cell_matrix_record, open( "masks/caMask_seed%s_pdivision%s_L%s_J%s_T%s_N%s.pkl"%(seed,p_division,L,J,T,N), "wb" ) )
plt.show()

def show_rgbvideo(timeseries_unstacked):
    time=0

    fig = plt.plot()
    rgb_timeseries=timeseries_unstacked # Read the numpy matrix with images in the rows
    T = len(rgb_timeseries)
    im=plt.imshow(rgb_timeseries[:,:,0].astype('uint8'), origin= 'lower', cmap='Greys')
    for time in range(len(rgb_timeseries[0,0,:])):
        im.set_data(rgb_timeseries[:,:,time].astype('uint8'))
        plt.xlabel(time)
        plt.pause(0.01)
    plt.show()

show_rgbvideo(cell_matrix_record)


    # print(np.shape(cell_matrix_record[0][0]))
#         # print(np.shape(cell_grid_list(int(80/2),:,0)))
# lenght_list = []
# for i in range(len(cell_matrix_record)):
#     lenght = np.count_nonzero(cell_matrix_record[i][int(80/2)][:])
#     lenght_list.append(lenght)
# lenght_list = np.array(lenght_list)/10
# plt.scatter(np.linspace(0,T,int(T/1)),lenght_list, c='k',s=1)
# plt.xlabel('Time (hours)')
# plt.ylabel('Colony diameter (mm)')
# plt.show()
#             # else:
#             #     plt.show()
#         # except ValueError:
#         #     print('!!!!!!!!!!!!!')
#         #     print('ValueError --> unstable solution')
#         #     print('!!!!!!!!!!!!!')
#         #     print()
#         #
#         #     pass
