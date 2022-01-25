import numpy as np
import copy
np.random.seed(2)

def ca(par_dict, initial_condition, L_x, L_y, J, I, T, N, circuit_n, boundary_coef=1, perturbation=0.001,timeresolution='fast',n_species = 6,tqdm_disable=True,division_index=2.3,seed=1):

    parent_list = [circuit1_eq, circuit2_eq,circuit3_eq,circuit4_eq,circuit5_eq,circuit6_eq,circuit7_eq]
    eq = parent_list[circuit_n - 1](par_dict)

    d_A = eq.d_A
    d_B = eq.d_B
    diffusion_rate_list = [d_A,d_B,0,0,0,0]

    date = datetime.date.today()

    # dx needs to be large and dt small
    dx = float(L_x) / float(J - 1) #size of cell
    dy = float(L_y) / float(I - 1)
    dt = float(T) / float(N - 1) #size of timestep

    x_grid = np.array([j * dx for j in range(J)])
    y_grid = np.array([i * dy for i in range(I)])
    t_grid = np.array([n * dt for n in range(N)])

    x_gridpoints = int(J/L_x) #x_gridpoints per space unit
    t_gridpoints = int(N/T) #t_gridpoints per timeunit

    species_list = []
    # if growth ==False:
    #     for index in range(n_species):
    #         species_list.append(np.random.uniform(low=steadystates[index] - perturbation, high=steadystates[index] + perturbation, size=(I, J)))
    for index in range(n_species):
        grid = np.zeros(shape=(I,J))
        grid[int(len(grid)/2),int(len(grid)/2)]=np.random.uniform(low=initial_condition[index] - perturbation, high=initial_condition[index] + perturbation, size=(1))
        species_list.append(grid)
    alpha_x = []
    alpha_y = []
    for index in range(n_species):
        alpha_x.append(diffusion_rate_list[index] * dt / (2. * dx * dx))
        alpha_y.append(diffusion_rate_list[index] * dt / (2. * dy * dy))


    X_new_list = []
    X_record_list = []
    for species_index in range(n_species):
        X_new_list.append(np.zeros([J, I]))
        if timeresolution != 'slow':
            X_record_list.append(np.zeros([J, I, T]))
        elif timeresolution == 'slow':
            X_record_list.append(np.zeros([J, I, N]))

    unittime = 0

    if n_species==6:
        A_new,B_new,C_new,D_new,E_new,F_new = X_new_list

    if n_species==2:
        A_new,B_new = X_new_list
    np.random.seed(seed)
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



    def p_division_evolution(final_time=158, t_gridpoints=1):
        ydata = [0.5,1,1.5,2,3.58,4.12,4.37,4.93,5.28,5.56,5.72,5.77]
        xdata = [0, 10, 20, 40, 65, 77, 102, 110, 124, 135, 149, 158]
        def sigmoid(x, L, x0, k, b):
            y = L / (1 +scipy.special.expit(-k * (x - x0))) + b
            # y = L / (1 + np.exp(-k * (x - x0))) + b
            return (y)



        p0 = [max(ydata), np.median(xdata), 1, min(ydata)]  # this is an mandatory initial guess
        popt, pcov = curve_fit(sigmoid, xdata, ydata, p0, method='lm')

        x = np.linspace(0,final_time,final_time*t_gridpoints)
        derivative_growth = derivative(sigmoid,x,args = (popt))

        return derivative_growth

    derivative_growth = p_division_evolution(final_time = T, t_gridpoints=t_gridpoints)
    cell_grid = np.where(species_list[0]!=0, 1, species_list[0])
    X_new_list=species_list
    cell_grid_list = []
    for n in tqdm(range(N),disable = tqdm_disable):

        # save results over time. these results can be saved a high time resolutions (save every time point) or low time resolution (save at every time unit).
        if timeresolution != 'slow':
            if n % (N / T) == 0:
                for species_index in range(n_species):
                    X_record_list[species_index][:, :, unittime] = species_list[species_index]
                unittime += 1

        elif timeresolution == 'slow':
            for species_index in range(n_species):
                X_record_list[species_index][:, :, unittime] = species_list[species_index]
            unittime += 1


        for species_index in range(n_species):
            X_new_list[species_index] = np.multiply(X_new_list[species_index],cell_grid)
        if n_species==6:
            A_new,B_new,C_new,D_new,E_new,F_new = X_new_list

        if n_species==2:
            A_new,B_new = X_new_list


        if n % (N / T) == 0: #only consider division at unit time (hour)

            p_division = derivative_growth[n]*2*x_gridpoints/division_index
            # print(p_division)
            X_new_list = cell_automata_colony(X_new_list,p_division)
            # print('dwfge')
            # print(cell_grid)
            # print(X_new_list[0])
            cell_grid = np.where(X_new_list[0]!=0, 1, X_new_list[0])
            cell_grid_list.append(cell_grid)
        species_list = X_new_list


    grids = (x_grid, y_grid, t_grid)
    records = X_record_list
    final_concentration = X_new_list


    return cell_grid,cell_grid_list
    # return records, final_concentration, grids, cell_grid
import sys
import os
pwd = os.getcwd()
root = pwd.split("home", 1)[0]
modelling_home = root + 'home/Documents/modelling'
modelling_ephemeral = root + 'ephemeral/Documents/modelling'
modulepath = modelling_home + '/3954/modules'
sys.path.append(modulepath)

from numerical_solvers_variableboundary import *
import pickle
#execution parameters
circuit_n=2
variant=0
parametersets_n = 1000
save_figure = False
tqdm_disable = False #disable tqdm
dimension = str(sys.argv[1])#str(sys.argv[1])
n_species=6
# open parameter dictionaries
general_df = pickle.load(open(modelling_home + '/3954/parameter_space_search/results/output_dataframes/lsa_df_circuit%r_variant%r_%rparametersets.pkl'%(circuit_n,variant,parametersets_n), "rb"))
#chose parameter sets to analyse
parID = str(sys.argv[2])
if parID == 'all':
    interesting_list = np.unique(general_df.index.get_level_values(0))
else:
    interesting_list = [int(parID)]

for parID in interesting_list:
    print('parID = ' + str(parID))
    mechanism = 'general'
    boundary_coef = 1 #1 is open boundary and 0 is closed boundary
    shape = 'ca'
    shape = 'growing_colony'
    growth = True

    # boundary_coef = 1 #1 is open boundary and 0 is closed boundary
    # shape = 'growing_colony'

    x_gridpoints = int(sys.argv[3])
    T =int(sys.argv[4])
    par_dict = general_df.loc[(parID,0)].to_dict()
    # print(par_dict)
    L=int(sys.argv[5])
    J = L *x_gridpoints  # number of equally spaced gridpoints in space domain (larger J means more spatial precision(tends towards continuum solution) )
    t_gridpoints = t_gridpoints_stability(L, J, T)  # numbesr of equally spaced gridpoints in domain (larger N means more temporal precision (tends towards continuum solution) )
    t_gridpoints = int(t_gridpoints)
    N = T * t_gridpoints
    initial_condition = [1000]*n_species
    filename = 'circuit%r_variant%r_boundary%r_%s_%sID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundary_coef, shape,mechanism,parID,L,J,T,N)

    if dimension == '2D':
     # Define 2D numerical parameters
        L_x = L
        L_y = L
        I = J
        try:

            # Run 2D simulation

            # records,final_concentration,grids = adi_ca(par_dict,initial_condition,L_x,L_y,J,I,T,N, circuit_n, boundary_coef=boundary_coef,tqdm_disable=tqdm_disable,n_species=n_species)
            cell_grid,cell_grid_list = ca(par_dict,initial_condition,L_x,L_y,J,I,T,N, circuit_n, shape,tqdm_disable=tqdm_disable,n_species=n_species,division_index=1.9,seed=1)
            plt.imshow(cell_grid)# plot_2D_final_concentration(final_concentration,grids,filename,n_species=n_species)
            plt.show()# plot_redgreen_contrast(final_concentration,L,mechanism,shape,filename,modelling_ephemeral,parID=parID,dimension=dimension,scale_factor=x_gridpoints,save_figure=save_figure)
            # rgb_timeseries = redgreen_contrast_timeseries(records)
            show_rgbvideo(cell_grid_list,parID)
            # if save_figure ==True:
            #     pickle.dump(final_concentration, open(modelling_ephemeral + '/3954/numerical_confocal/results/simulation/1M_colony_ca/2D/2Dfinal_%s.pkl'%filename, 'wb'))
            #     pickle.dump(records,open(modelling_ephemeral + '/3954/numerical_confocal/results/simulation/1M_colony_ca/2D/2Dtimeseries_%s.pkl'%filename, 'wb'))

            # else:
            #     plt.show()
        except ValueError:
            print('!!!!!!!!!!!!!')
            print('ValueError --> unstable solution')
            print('!!!!!!!!!!!!!')
            print()

            pass
