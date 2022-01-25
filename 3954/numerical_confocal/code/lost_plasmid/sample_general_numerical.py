# %config Completer.use_jedi = False
import sys
import os
pwd = os.getcwd()
root = pwd.split("home", 1)[0]
modelling_home = root + 'home/Documents/modelling'
modelling_ephemeral = root + 'ephemeral/Documents/modelling'
modulepath = modelling_home + '/6eq/modules'
sys.path.append(modulepath)

from numerical_solvers_variableboundary import *
#execution parameters
circuit_n=4
variant=0
parametersets_n = 10000
save_figure = False
tqdm_disable = False #disable tqdm
dimension = str(sys.argv[1])
n_species=4

# open parameter dictionaries
general_df = pickle.load(open(modelling_home + '/6eq/parameter_space_search/parameterfiles/df_circuit2_variant0_10000parametersets.pkl', "rb"))
print(general_df)# general_df = pickle.load(open(modelling_home + '/6eq/parameter_space_search/results/output_dataframes/lsa_df_circuit%r_variant%r_%rparametersets.pkl'%(circuit_n,variant,parametersets_n), "rb"))
#chose parameter sets to analyse


parID = str(sys.argv[2])
if parID == 'all':
    interesting_list = np.unique(general_df.index.get_level_values(0))
else:
    interesting_list = [int(parID)]

for parID in interesting_list:
    print('parID = ' + str(parID))
    mechanism = 'general'
    boundary_coef = 0 #1 is open boundary and 0 is closed boundary
    shape = 'no_growth'
    growth = False


    # boundary_coef = 1 #1 is open boundary and 0 is closed boundary
    # shape = 'growing_colony'

    x_gridpoints = int(sys.argv[3])
    T =int(sys.argv[4])
    # par_dict = general_df.loc[(parID,0)].to_dict()
    par_dict = general_df.loc[parID].to_dict()
    print(par_dict)
    L=8
    J = L *x_gridpoints  # number of equally spaced gridpoints in space domain (larger J means more spatial precision(tends towards continuum solution) )
    t_gridpoints = t_gridpoints_stability(L, J, T)  # number of equally spaced gridpoints in domain (larger N means more temporal precision (tends towards continuum solution) )

    N = T * t_gridpoints
    initial_condition = [0.001]*n_species

    filename = 'circuit%r_variant%r_boundary%r_%s_%sID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundary_coef, shape,mechanism,parID,L,J,T,N)

    if dimension == '1D':
        try:
            records, final_concentration, grids = crank_nicolson(par_dict, initial_condition, L, J, T, N, circuit_n, n_species=n_species,boundary_coef=boundary_coef,growth = growth,tqdm_disable=tqdm_disable)
            # plot_1D_final_concentration(final_concentration, grids,mechanism,shape,filename,parID,save_figure=save_figure,path=modelling_hpc)
            plot_redgreen_contrast(final_concentration,L,mechanism,shape,filename,modelling_ephemeral,circuit_n=circuit_n,dimension=dimension,save_figure=False)
            if save_figure ==True:
                pickle.dump(final_concentration, open(modelling_ephemeral + '/6eq/numerical_confocal/results/simulation/1M_growingcolony/1Dfinal_%s.pkl'%filename, 'wb'))
                pickle.dump(records,open(modelling_ephemeral + '/6eq/numerical_confocal/results/simulation/1M_growingcolony/1Dtimeseries_%s.pkl'%filename, 'wb'))
                # timeseries_1D_rgb(records, L, J,T,filename,modelling_hpc,save_figure=save_figure)

        except ValueError:
            print('!!!!!!!!!!!!!')
            print('ValueError --> unstable solution')
            print('!!!!!!!!!!!!!')
            print()

            pass


    if dimension == '2D':
     # Define 2D numerical parameters
        L_x = L
        L_y = L
        I = J

        # Run 2D simulation
        try:
            records,final_concentration,grids = adi_shape(par_dict,initial_condition,L_x,L_y,J,I,T,N, circuit_n,shape, boundary_coef=boundary_coef,tqdm_disable=tqdm_disable)#,tqdm_disable=tqdm_disable)
            plot_redgreen_contrast(final_concentration,L,mechanism,shape,filename,modelling_ephemeral,dimension=dimension,scale_factor=x_gridpoints,save_figure=save_figure)
            if save_figure ==True:
                pickle.dump(final_concentration, open(modelling_ephemeral + '/6eq/numerical_confocal/results/simulation/1M_growingcolony/2Dfinal_%s.pkl'%filename, 'wb'))
                pickle.dump(records,open(modelling_ephemeral + '/6eq/numerical_confocal/results/simulation/1M_growingcolony/2Dtimeseries_%s.pkl'%filename, 'wb'))

            # else:
            #     plt.show()

        except ValueError:
            print('!!!!!!!!!!!!!')
            print('ValueError --> unstable solution')
            print('!!!!!!!!!!!!!')
            print()

            pass
