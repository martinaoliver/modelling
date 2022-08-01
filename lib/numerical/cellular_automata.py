import numpy
import matplotlib as mpl
mpl.use('tkagg')

import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import spdiags, diags
import matplotlib.pyplot as plt
from tqdm import tqdm
import copy
from scipy.sparse import linalg
from scipy.linalg import solve_banded
from class_circuit_eq import *
import sys

def show_rgbvideo(timeseries_unstacked,parID):
    time=0

    fig = plt.plot()
    rgb_timeseries=timeseries_unstacked # Read the numpy matrix with images in the rows
    im=plt.imshow(rgb_timeseries[0].astype('uint8'), origin= 'lower')
    for time in range(len(rgb_timeseries)):
        im.set_data(rgb_timeseries[time].astype('uint8'))
        plt.title(parID)
        plt.pause(0.01)
    plt.show()

n_species=3
parID= int(sys.argv[1]);

t_gridpoints = int(sys.argv[2])
x_gridpoints = int(sys.argv[3])
T =int(sys.argv[4])
L_x =int(sys.argv[5])

J = L_x*x_gridpoints;  L_y=L_x; I=J
N = T*t_gridpoints
U0 = []
perturbation=0.001
steadystates=[0.1]*n_species
cell_matrix = np.zeros(shape=(I,J))
cell_matrix[int(I/2), int(J/2)] = 1
for index in range(n_species):
    U0.append(np.random.uniform(low=steadystates[index] - perturbation, high=steadystates[index] + perturbation, size=(I, J)))
U_new = U0*cell_matrix
cell_matrix_list = []
def check_neighbours(cell_matrix,y_pos,x_pos): #returns grid with the neighbouring points
    top_array = [cell_matrix[y_pos-1, x_pos-1], cell_matrix[y_pos-1,x_pos], cell_matrix[y_pos-1,x_pos+1]]
    middle_array = [cell_matrix[y_pos, x_pos-1], np.nan, cell_matrix[y_pos,x_pos+1]]
    # print('afe')
    # print(cell_matrix[2,1])
    bottom_array = [cell_matrix[y_pos+1, x_pos-1], cell_matrix[y_pos+1,x_pos], cell_matrix[y_pos+1,x_pos+1]]
    neighbours_cellmatrix = np.array([top_array,middle_array,bottom_array])
    return neighbours_cellmatrix
def cell_automata_colony(species_list,cell_matrix, p_division):
    new_species_list = copy.deepcopy(species_list)
    original_species_list = copy.deepcopy(species_list)
    cell_matrix_new = copy.deepcopy(cell_matrix)
    for y_pos in np.linspace(1,len(cell_matrix)-2,len(cell_matrix)-2):
    # for x_pos in np.linspace(0,len(cell_matrix)-2,len(cell_matrix)-2):
        for x_pos in np.linspace(1,len(cell_matrix)-2,len(cell_matrix)-2):
        # for y_pos in np.linspace(0,len(cell_matrix)-2,len(cell_matrix)-2):
            y_pos = int(y_pos)
            x_pos = int(x_pos)
            # print(y_pos,x_pos,cell_matrix)
            if cell_matrix[y_pos, x_pos]!=0:
                neighbours_cellmatrix = check_neighbours(cell_matrix,y_pos,x_pos)
                if 0 in neighbours_cellmatrix:
                    cell_division=np.random.choice([1,0],p=[p_division,1-p_division])
                    if cell_division==1:
                        index_nocells=np.where(np.array(neighbours_cellmatrix )== 0)
                        divided_cell_index = np.random.choice(range(len(index_nocells[0])))
                        index_newcell_y, index_newcell_x = (index_nocells[n][divided_cell_index] for n in range(2))
                        # index_newcell_y , index_newcell_x = 1,2

                        # for species,new_species in zip(original_species_list,new_species_list):
                        for count,species in enumerate(original_species_list):
                            new_species_list[count][index_newcell_y+y_pos-1,index_newcell_x+x_pos-1] += species[y_pos,x_pos]/2
                            new_species_list[count][y_pos,x_pos] += species[y_pos,x_pos]/2
                            # count+=1
                        cell_matrix_new[index_newcell_y+y_pos-1,index_newcell_x+x_pos-1]=1


    return new_species_list, cell_matrix_new
p_division=0.7
print (p_division)
np.random.seed(1)
for i in range(120):        #predict if division occurs based on the p_division, the current cell matrix
    #return new cell matrix and updated concentrations with dilution
    U_new, cell_matrix_new = cell_automata_colony(U_new, cell_matrix, p_division)
    cell_matrix = copy.deepcopy(cell_matrix_new)
    cell_matrix_list.append(cell_matrix)

plt.imshow(cell_matrix)
plt.colorbar()
plt.show()

# show_rgbvideo(cell_matrix_list,parID)
