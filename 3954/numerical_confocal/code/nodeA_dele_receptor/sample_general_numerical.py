import sys
import os
pwd = os.getcwd()
root = pwd.split("home", 1)[0]
modelling_home = root + 'home/Documents/modelling'
modelling_ephemeral = root + 'ephemeral/Documents/modelling'
modulepath = modelling_home + '/3954/modules'
sys.path.append(modulepath)

from numerical_solver import *
import pickle
#execution parameters
circuit_n=10
variant=7
# parametersets_n = 1000000
save_figure = False
tqdm_disable = False #disable tqdm
dimension = str(sys.argv[1])#str(sys.argv[1])
n_species=5
# open parameter dictionaries
# general_df = pickle.load(open(modelling_home + '/3954/parameter_space_search/results/output_dataframes/lsa_df_circuit%r_variant%r_%rparametersets.pkl'%(circuit_n,variant,parametersets_n), "rb"))
# general_df = pickle.load(open(modelling_home + '/3954/parameter_space_search/results/output_dataframes/interesting_variations/ATC_24240.pkl', "rb"))
# general_df = pickle.load(open(modelling_home + '/3954/parameter_space_search/results/output_dataframes/interesting_variations/ATC_2852.pkl', "rb"))
# general_df = pickle.load(open(modelling_home + '/3954/parameter_space_search/code/interesting_parIDs/interesting_parIDs.pkl', "rb"))
# print(general_df['kce'])#chose parameter sets to analyse
parID = str(sys.argv[2])
# if parID == 'all':
#     interesting_list = np.unique(general_df.index.get_level_values(0))
# else:
#     interesting_list = [int(parID)]

# for parID in interesting_list:
for parID in [1]:
    print('parID = ' + str(parID))
    mechanism = 'nodeA_receptor'
    boundary_coef = 1 #1 is open boundary and 0 is closed boundary
    shape = 'ca'
    growth = True

    # boundary_coef = 1 #1 is open boundary and 0 is closed boundary
    # shape = 'growing_colony'

    x_gridpoints = int(sys.argv[3])
    T =int(sys.argv[4])
    # par_dict = general_df.loc[(parID,0)].to_dict()
    # par_dict = general_df.loc[parID].to_dict()
    par_dict = {'bb':0.95,'bd':0.95,'be':0.95,'bf':0.95, 'Vb': 0.91,'Vd': 0.91,'Ve': 0.91,'Vf': 0.91,'keb': 0.91, 'khd': 0.91,'kfe': 0.91,'kee': 0.91,'khf': 0.91,'muB':0.45,'muH':0.45,'muD':0.45, 'muE':0.45,'muF':0.45,'alphaH':0.45,'d_H':0.2,'n':0.45,'kon':2,'koff':2,'R':100} #circuit10

    # par_dict['bb']=10000
    print(par_dict)
    # print(par_dict)
    L=int(sys.argv[5])
    J = L *x_gridpoints  # number of equally spaced gridpoints in space domain (larger J means more spatial precision(tends towards continuum solution) )
    t_gridpoints = t_gridpoints_stability(L, J, T)  # number of equally spaced gridpoints in domain (larger N means more temporal precision (tends towards continuum solution) )
    t_gridpoints = int(t_gridpoints)
    N = T * t_gridpoints
    initial_condition = [100]*n_species#0.001
    p_division=1#0.41
    seed=1
    # filename = 'circuit%r_variant%r_boundary%r_%s_%sID%r_L%r_J%r_T%r_N%r_divisionindex%r'%(circuit_n,variant,boundary_coef, shape,mechanism,parID,L,J,T,N,division_index)
    filename = 'circuit%r_variant%r_boundary%r_%s_%sID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundary_coef, shape,mechanism,parID,L,J,T,N)
    print(filename)
    if dimension == '1D':
        try:
            records, final_concentration, grids = crank_nicolson_ca(par_dict, initial_condition, L, J, T, N, circuit_n,n_species=n_species,boundary_coef=boundary_coef,growth = growth,tqdm_disable=tqdm_disable)
            # plot_1D_final_concentration(final_concentration,grids,mechanism,shape,filename,n_species, parID,save_figure=save_figure,path=modelling_home)
            plot_redgreen_contrast(final_concentration,L,mechanism,shape,filename,modelling_ephemeral,dimension=dimension,save_figure=False)
            if save_figure ==True:
                pickle.dump(final_concentration, open(modelling_ephemeral + '/3954/numerical_confocal/results/simulation/1M_colony_ca/1D/1Dfinal_%s.pkl'%filename, 'wb'))
                pickle.dump(records,open(modelling_ephemeral + '/3954/numerical_confocal/results/simulation/1M_colony_ca/1D/1Dtimeseries_%s.pkl'%filename, 'wb'))
                timeseries_1D_rgb(records, L, J,T,filename,modelling_home,save_figure=save_figure)
        except ValueError:
            print('!!!!!!!!!!!!!')
            print('ValueError --> unstable solution')
            print('!!!!!!!!!!!!!')
            print()

            pass

        # plot_1D_final_concentration(final_concentration, grids,mechanism,shape,filename,parID,save_figure=save_figure,path=modelling_home)
        # plot_redgreen_contrast(final_concentration,grids,mechanism,shape,filename,parID,modelling_home,dimension=dimension)
        # plt.show()
    if dimension == '2D':
     # Define 2D numerical parameters
        L_x = L
        L_y = L
        I = J
        # try:

            # Run 2D simulation
        records,final_concentration,grids = adi_ca(par_dict,initial_condition,L_x,L_y,J,I,T,N, circuit_n, boundary_coef=boundary_coef,tqdm_disable=tqdm_disable,n_species=n_species,p_division=p_division,seed=seed)
        print(final_concentration[0])
        # plot_2D_final_concentration(final_concentration,grids,filename,n_species=n_species)
        plot_redgreen_contrast(final_concentration,L,mechanism,shape,filename,modelling_ephemeral,parID=0,dimension=dimension,scale_factor=x_gridpoints,save_figure=save_figure)
        # rgb_timeseries = redgreen_contrast_timeseries(records)
        # show_rgbvideo(rgb_timeseries)
        if save_figure ==True:
            pickle.dump(final_concentration, open(modelling_ephemeral + '/3954/numerical_confocal/results/simulation/1M_colony_ca/2D/2Dfinal_%s.pkl'%filename, 'wb'))
            pickle.dump(records,open(modelling_ephemeral + '/3954/numerical_confocal/results/simulation/1M_colony_ca/2D/2Dtimeseries_%s.pkl'%filename, 'wb'))

            # else:
            # else:
            #     plt.show()
        # except ValueError:
        #     print('!!!!!!!!!!!!!')
        #     print('ValueError --> unstable solution')
        #     print('!!!!!!!!!!!!!')
        #     print()
        #
        #     pass
