import sys
import os
pwd = os.getcwd()
root = pwd.split("home", 1)[0]
modelling_home = root + 'home/Documents/modelling'
modelling_ephemeral = root + 'ephemeral/Documents/modelling'
modulepath = modelling_home + '/3954/modules'
sys.path.append(modulepath)

# from numerical_solvers_variableboundary import *
from adi_v1_openclosed_ca_function import *
from plotting_numerical import *
import pickle
#execution parameters

mechanism = 'general'
shape = 'ca'
parID = int(sys.argv[1])
circuit_n=8
variant='5716gaussian'
parametersets_n = 30000 #1000000
save_figure = False
tqdm_disable = False #disable tqdm
n_species=3

# open parameter dictionaries
general_df= pickle.load( open(modelling_home + '/3954/parameter_space_search/parameterfiles/df_circuit%r_variant%s_%rparametersets.pkl'%(2,variant,parametersets_n), "rb" ) )
par_dict = general_df.loc[parID].to_dict()


#solver parameters
L_x=int(sys.argv[2]); x_gridpoints = int(sys.argv[3]); J = L_x*10;  L_y=L_x; I=J
T =int(sys.argv[4]); t_gridpoints = 5 ; N = T*t_gridpoints

filename = 'circuit%r_variant%s_%s_%sID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant, shape,mechanism,parID,L_x,J,T,N)#,p_division,kce)
savefig_path = ''
# Run 2D simulation
try:
    U_record,U_final = adi_ca(par_dict,L_x,L_y,J,I,T,N, circuit_n,n_species, tqdm_disable=tqdm_disable)#,p_division=p_division,seed=seed)
        # plot_2D_final_concentration(final_concentration,grids,filename,n_species=n_species)
        # savefig_path = modelling_ephemeral + '/3954/numerical_confocal/results/figures/1M_colony_ca/2D/5716gaussian'
    plot_redgreen_contrast(U_final,L_x,mechanism,shape,filename,savefig_path,parID=parID,scale_factor=x_gridpoints,save_figure=save_figure)
        # plot_redgreen_contrast(final_concentration,L,mechanism,shape,filename,modelling_ephemeral,parID='%r-%r'%(parID,kce),dimension=dimension,scale_factor=x_gridpoints,save_figure=save_figure)
        # rgb_timeseries = redgreen_contrast_timeseries(records)
        # show_rgbvideo(rgb_timeseries)
        # if save_figure ==True:
        #     pickle.dump(final_concentration, open(modelling_ephemeral + '/3954/numerical_confocal/results/simulation/1M_colony_ca/2D/5716gaussian/2Dfinal_%s.pkl'%filename, 'wb'))
        #     pickle.dump(records,open(modelling_ephemeral + '/3954/numerical_confocal/results/simulation/1M_colony_ca/2D/5716gaussian/2Dtimeseries_%s.pkl'%filename, 'wb'))

        # else:
        # else:
        #     plt.show()
except ValueError:
    print('!!!!!!!!!!!!!')
    print('ValueError --> unstable solution')
    print('!!!!!!!!!!!!!')
    print()

    pass
