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
circuit_n=8
variant=6
parametersets_n = 1000
save_figure = True
tqdm_disable = False #disable tqdm
dimension = str(sys.argv[1])#str(sys.argv[1])
n_species=6
# open parameter dictionaries
# general_df = pickle.load(open(modelling_home + '/3954/parameter_space_search/results/output_dataframes/lsa_df_circuit%r_variant%r_%rparametersets.pkl'%(circuit_n,variant,parametersets_n), "rb"))
#chose parameter sets to analyse
parID = int(sys.argv[2])
# if parID == 'all':
#     interesting_list = np.unique(general_df.index.get_level_values(0))
# else:
interesting_list = [int(parID)]

for parID in interesting_list:
    print('parID = ' + str(parID))
    mechanism = 'nodeAdele'
    boundary_coef = 1 #1 is open boundary and 0 is closed boundary
    shape = 'ca'
    growth = True

    # boundary_coef = 1 #1 is open boundary and 0 is closed boundary
    # shape = 'growing_colony'

    x_gridpoints = int(sys.argv[3])
    T =int(sys.argv[4])
    L=int(sys.argv[5])
    p_division=float(sys.argv[6])
    J = L *x_gridpoints  # number of equally spaced gridpoints in space domain (larger J means more spatial precision(tends towards continuum solution) )
    t_gridpoints = t_gridpoints_stability(L, J, T)  # number of equally spaced gridpoints in domain (larger N means more temporal precision (tends towards continuum solution) )
    t_gridpoints = int(t_gridpoints)
    N = T * t_gridpoints
    initial_condition = [0.001]*n_species
    filename = 'circuit%r_variant%r_boundary%r_%s_%sID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundary_coef, shape,mechanism,parID,L,J,T,N)
    save_path = modelling_home + '/3954/numerical_confocal/results/simulation/1M_colony_ca/2D/nodeA_dele'
    final_concentration = pickle.load(open(modelling_ephemeral + '/3954/numerical_confocal/results/simulation/1M_colony_ca/2D/nodeA_dele/2Dfinal_%s.pkl'%filename, 'rb'))
    plot_redgreen_contrast(final_concentration,L,mechanism,shape,filename,save_path,parID=parID,dimension=dimension,scale_factor=x_gridpoints,save_figure=save_figure)
