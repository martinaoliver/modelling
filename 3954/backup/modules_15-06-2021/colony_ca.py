import numpy as np
import copy

def check_neighbours(grid,x_pos,y_pos): #returns grid with the neighbouring points
    top_array = [grid[y_pos-1, x_pos-1], grid[y_pos-1,x_pos], grid[y_pos-1,x_pos+1]]
    middle_array = [grid[y_pos, x_pos-1], np.nan, grid[y_pos,x_pos+1]]
    bottom_array = [grid[y_pos+1, x_pos-1], grid[y_pos+1,x_pos], grid[y_pos+1,x_pos+1]]
    neighbours_grid = np.array([top_array,middle_array,bottom_array])
    return neighbours_grid

def cell_automata_colony(species_list,p_division):
    new_species_list = copy.deepcopy(species_list)
    original_species_list = copy.deepcopy(species_list)
    grid=original_species_list[0]

    for x_pos in np.linspace(1,len(grid)-2,len(grid)-2):
        for y_pos in np.linspace(1,len(grid)-2,len(grid)-2):
            x_pos = int(x_pos)
            y_pos = int(y_pos)
            if grid[x_pos,y_pos]!=0:
                neighbours_grid = check_neighbours(grid,x_pos,y_pos)
                if 0 in neighbours_grid:
                    cell_division=np.random.choice([1,0],p=[p_division,1-p_division])
                    if cell_division==1:
                        index_nocells=np.where(np.array(neighbours_grid )== 0)

                        index_newcell_y, index_newcell_x = [np.random.choice(index_nocells[0]),np.random.choice(index_nocells[1])]
                        count=0
                        for species,new_species in zip(original_species_list,new_species_list):
                            new_species_list[count][index_newcell_y+y_pos-1,index_newcell_x+x_pos-1] = species[y_pos,x_pos]/2
                            new_species_list[count][y_pos,x_pos] = species[y_pos,x_pos]/2
                            count+=1


    return new_species_list

grid1 = np.zeros(shape=(3,3))
grid2 = np.zeros(shape=(3,3))
grid1[int(len(grid1)/2),int(len(grid1)/2)]=2
grid2[int(len(grid2)/2),int(len(grid1)/2)]=2
species_list = [grid1,grid1,grid2]


for t in range(1):
    species_list = cell_automata_colony(species_list,1)
    cell_grid = np.where(species_list)
    print(species_list)
# print (new_species_list[0], cell_grid_list[0])

np.array([[1,1,1],[1,0,1],[2,2,2]])
