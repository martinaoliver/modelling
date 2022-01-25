
modulepath = '/Users/mo2016/Documents/modelling/6eq/modules'
import sys
sys.path.append(modulepath)
modelling_hpc = '/Volumes/mo2016/home/Documents/modelling'
modelling_ephemeral = '/Volumes/mo2016/ephemeral/Documents/modelling'
while True:
    try:
        from numerical_solvers_variableboundary import *
        break
    except ImportError:
        modelling_hpc = '/rds/general/user/mo2016/home/Documents/modelling'
        modulepath = modelling_hpc + '/6eq/modules'
        modelling_ephemeral = '/rds/general/user/mo2016/ephemeral/Documents/modelling'
        sys.path.append(modulepath)

from numerical_solvers_variableboundary import *
circuit_n=2
variant= 0

# Specifiy number of parameter sets in parameterset file to be loaded
n_param_sets = 1000

df= pickle.load( open(modelling_hpc + '/6eq/parameter_space_search/parameterfiles/df_circuit%r_variant%r_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), "rb" ) )

df

def numerical_check(start_batch_index,n_param_sets,df, dimension,x_gridpoints,T,circuit_n=2, variant=0, n_species=6):
    save_figure = True
    tqdm_disable = True #disable tqdm

    df_index = np.unique(df.index.get_level_values(0))
    for parID in df_index:
        print('parID = ' + str(parID))
        mechanism = 'general'
        boundary_coef = 1 #1 is open boundary and 0 is closed boundary
        shape = 'growing_colony'

        # x_gridpoints = int(sys.argv[3])
        # T =int(sys.argv[4])
        par_dict = df.loc[(parID,0)].to_dict()
        L=8
        J = L *x_gridpoints  # number of equally spaced gridpoints in space domain (larger J means more spatial precision(tends towards continuum solution) )
        t_gridpoints = t_gridpoints_stability(L, J, T)  # number of equally spaced gridpoints in domain (larger N means more temporal precision (tends towards continuum solution) )

        N = T * t_gridpoints
        initial_condition = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001]

        filename = 'circuit%r_variant%r_boundary%r_%s_%sID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundary_coef, shape,mechanism,parID,L,J,T,N)

        if dimension == '1D':
            try:
                records, final_concentration, grids = crank_nicolson(par_dict, initial_condition, L, J, T, N, circuit_n,boundary_coef=boundary_coef,tqdm_disable=tqdm_disable)
                # plot_1D_final_concentration(final_concentration, grids,mechanism,shape,filename,parID,save_figure=save_figure,path=modelling_hpc)
                plot_redgreen_contrast(final_concentration,L,mechanism,shape,filename,modelling_ephemeral,dimension=dimension,save_figure=save_figure)
                if save_figure ==True:
                    pickle.dump(final_concentration, open(modelling_hpc + '/6eq/numerical_confocal/results/simulation/1M_growingcolony/1Dfinal_%s.pkl'%filename, 'wb'))

            except ValueError:
                print('!!!!!!!!!!!!!')
                print('ValueError --> unstable solution')
                print('!!!!!!!!!!!!!')
                print()

                pass

            # plot_1D_final_concentration(final_concentration, grids,mechanism,shape,filename,parID,save_figure=save_figure,path=modelling_hpc)
            # plot_redgreen_contrast(final_concentration,grids,mechanism,shape,filename,parID,modelling_hpc,dimension=dimension)
            # plt.show()
        if dimension == '2D':
         # Define 2D numerical parameters
            L_x = L
            L_y = L
            I = J

            # Run 2D simulation
            try:
                records,final_concentration,grids = adi_shape(par_dict,initial_condition,L_x,L_y,J,I,T,N, circuit_n,shape, boundary_coef=boundary_coef,tqdm_disable=tqdm_disable)#,tqdm_disable=tqdm_disable)
                plot_redgreen_contrast(final_concentration,L,mechanism,shape,filename, modelling_ephemeral,dimension=dimension,scale_factor=x_gridpoints,save_figure=save_figure)
                if save_figure ==True:
                    # plt.savefig(modelling_hpc + '/6eq/numerical_confocal/results/figures/redgreen/%s/%s/redgreen_%s.png' % (mechanism,shape,filename))
                    pickle.dump(final_concentration, open(modelling_hpc + '/6eq/numerical_confocal/results/simulation/1M_growingcolony/2Dfinal_%s.pkl'%filename, 'wb'))
                    # plt.show()
                # else:
                #     plt.show()

            except ValueError:
                print('!!!!!!!!!!!!!')
                print('ValueError --> unstable solution')
                print('!!!!!!!!!!!!!')
                print()

                pass
