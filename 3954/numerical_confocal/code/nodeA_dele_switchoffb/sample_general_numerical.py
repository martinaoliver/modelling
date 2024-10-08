import sys
import os
import sys
from adi_switchoffb import *

pwd = os.getcwd()
root = pwd.split("home", 1)[0]
modelling_home = root + 'home/Documents/modelling'
modelling_ephemeral = root + 'ephemeral/Documents/modelling'
modulepath = modelling_home + '/3954/modules/new_CN'
sys.path.append(modulepath)

if root == '/Volumes/mo2016/':
    import matplotlib
    matplotlib.use('TKAgg')
# from numerical_solvers_variableboundary import *
from plotting_numerical import *
import pickle
#execution parameters

mechanism = 'nodeAdele_switchoffb'
shape = 'ca'
parID = int(sys.argv[1])
circuit_n=8
variant=0
parametersets_n = 10000 #1000000
save_figure = True
tqdm_disable = False #disable tqdm
n_species=3

# open parameter dictionaries
general_df= pickle.load( open(modelling_home + '/3954/parameter_space_search/parameterfiles/df_circuit%r_variant%s_%rparametersets.pkl'%(2,variant,parametersets_n), "rb" ) )
par_dict = general_df.loc[parID].to_dict()
par_dict['bc']=1
par_dict['be']=1
par_dict['bf']=1

d_B = par_dict['d_B']
D = np.zeros(n_species)
D[0]=d_B

#solver parameters
L_x=int(sys.argv[2]); x_gridpoints = int(sys.argv[3]); J = L_x*x_gridpoints;  L_y=L_x; I=J
T =int(sys.argv[4]); t_gridpoints = int(sys.argv[5]) ; N = T*t_gridpoints
p_division=float(sys.argv[6])
filename = 'circuit%r_variant%s_%s_%sID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant, shape,mechanism,parID,L_x,J,T,N)#,p_division,kce)
savefig_path = modelling_ephemeral + '/3954/numerical_confocal/results/figures/1M_colony_ca/2D/nodeA_dele_newCN'
try:
    U_record,U_final = adi_ca(par_dict,L_x,L_y,J,I,T,N, circuit_n,n_species,D, tqdm_disable=tqdm_disable, p_division=p_division)#,p_division=p_division,seed=seed)
        # plot_2D_final_concentration(final_concentration,grids,filename,n_species=n_species)
        # savefig_path = modelling_ephemeral + '/3954/numerical_confocal/results/figures/1M_colony_ca/2D/5716gaussian'
    plot_redgreen_contrast(U_final,L_x,mechanism,shape,filename,savefig_path,parID=parID,scale_factor=x_gridpoints,save_figure=save_figure)
    plot_redgreen_contrast(U_final,L_x,mechanism,shape,filename,savefig_path,parID=parID,scale_factor=x_gridpoints,save_figure=False)

    rgb_timeseries = redgreen_contrast_timeseries(U_record)
    show_rgbvideo(rgb_timeseries,parID)
    if save_figure ==True:
        pickle.dump(U_final, open(modelling_ephemeral + '/3954/numerical_confocal/results/simulation/1M_colony_ca/2D/nodeA_dele_newCN/2Dfinal_%s.pkl'%filename, 'wb'))
        pickle.dump(U_record,open(modelling_ephemeral + '/3954/numerical_confocal/results/simulation/1M_colony_ca/2D/nodeA_dele_newCN/2Dtimeseries_%s.pkl'%filename, 'wb'))
        print(filename)
        # else:
        # else:
        #     plt.show()
    # print(np.shape(U_record))
    # print(np.shape(U_record[15 120 2][:,:,-1]))
    print(U_record[2][:,:,-1])
    print(U_final[2])
except ValueError:
    print('!!!!!!!!!!!!!')
    print('ValueError --> unstable solution')
    print('!!!!!!!!!!!!!')
    print()

    pass
