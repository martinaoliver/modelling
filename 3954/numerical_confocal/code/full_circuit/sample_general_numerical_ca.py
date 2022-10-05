#############
###paths#####
#############
import sys
import os

from importlib_metadata import distribution
pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')
#############



from numerical.adi_ca_function_openclosed_nodilution import *
from numerical.plotting_numerical import *
import pickle
#execution parameters

mechanism = 'fullcircuit'
shape = 'ca'
parID = int(sys.argv[1])

circuit_n=2
variant= 0#9
folder = 'fullcircuit/1M_turingI'#'fullcircuit/1M'#'fullcircuit/1M_turingI'
n_species = 6

parametersets_n = 30000 #1000000
save_figure = True
tqdm_disable = False #disable tqdm
boundarycoeff = float(sys.argv[6])
seed=1;p_division=0.41#0.147#0.5
# open parameter dictionaries
# general_df= pickle.load( open(modelling_home + '/3954/parameter_space_search/parameterfiles/5716gaussian/df_circuit%r_variant%s_%rparametersets.pkl'%(circuit_n,variant,parametersets_n), "rb" ) )
general_df= pickle.load( open(modelling_home + "/3954/parameter_space_search//results/turing_dataframes/turing_lsa_df_circuit%r_variant%r_%rparametersets.pkl"%(circuit_n,variant,1000000), "rb" ) )
par_dict = general_df.loc[parID].to_dict()
# par_dict = general_df.iloc[parID].to_dict()

d_A = par_dict['d_A']
d_B = par_dict['d_B']

D = np.zeros(n_species)
D[0]=d_A
d_B=1
D[1]=d_B
print(par_dict)
par_dict['mulva'] = par_dict['mulva'] + np.log(2)*p_division


#solver parameters
L_x=int(sys.argv[2]); x_gridpoints = int(sys.argv[3]); L=L_x; J = L*x_gridpoints;  L_y=L_x; I=J
T =int(sys.argv[4]); t_gridpoints = int(sys.argv[5]) ; N = T*t_gridpoints
suggested_tgridpoints = x_gridpoints**2

# filename = 'circuit%r_variant%s_%s_%sID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant, shape,mechanism,parID,L,J,T,N)
filename = 'circuit%r_variant%s_bc%s_%s_%sID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundarycoeff, shape,mechanism,parID,L,J,T,N)
savefig_path = 'test_results'
try:
    U_record,U_final = adi_ca(par_dict,L_x,L_y,J,I,T,N, circuit_n,n_species,D, seed=1, p_division=p_division, tqdm_disable=tqdm_disable, growth='Slow')#,p_division=p_division,seed=seed)
    # U_record,U_final = adi_ca_openclosed_nodilution(par_dict,L_x,L_y,J,I,T,N, circuit_n,n_species,D, seed=seed, p_division=p_division, tqdm_disable=tqdm_disable, growth='Fast', boundarycoeff=boundarycoeff)#,p_division=p_division,seed=seed)
    # plot_2D_final_concentration(U_final,L_x,J,filename,savefig_path,n_species=n_species,save_figure=save_figure)
    # mask=pickle.load( open( modelling_home + "/3954/numerical_confocal/code/cellular_automata_templates/masks/caMask_seed%s_pdivision%s_L%s_J%s_T%s_N%s.pkl"%(seed,p_division,L,J,T,N), "rb" ) )
    # plot_redgreenblue_contrast(U_final,L_x,mechanism,shape,filename,parID=parID,scale_factor=x_gridpoints,save_figure=save_figure)

  
    plot_redgreen_contrast(U_final,L_x,mechanism,shape,filename,savefig_path,parID=parID,scale_factor=x_gridpoints,save_figure=save_figure)
    # rgb_timeseries = redgreen_contrast_timeseries(records)
    # show_rgbvideo(rgb_timeseries)
    if save_figure ==True:
        print('saved')
        # pickle.dump(U_final, open('test_results/2Dfinal_%s.pkl'%(filename), 'wb'))
        # pickle.dump(U_record,open(modelling_home + '/3954/numerical_confocal/results/simulation/ca/2D/full_circuit/%s/2Dtimeseries_%s.pkl'%(folder,filename), 'wb'))

        # savefig_path = modelling_ephemeral + '/3954/numerical_confocal/results/figures/1M_colony_ca/2D/5716gaussian'
    # plot_redgreen_contrast(U_final,L_x,mechanism,shape,filename,savefig_path,parID=parID,scale_factor=x_gridpoints,save_figure=save_figure)
    # plot_redgreen_contrast(U_final,L_x,mechanism,shape,filename,savefig_path,parID=parID,scale_factor=x_gridpoints,save_figure=False)
    plt.imshow(U_final[-1])
    plt.show()
    plt.imshow(U_final[0])
    plt.show()
    rgb_timeseries = redgreen_contrast_timeseries(U_record)
    show_rgbvideo(rgb_timeseries,parID)
    # if save_figure ==True:
    #     pickle.dump(U_final, open(modelling_ephemeral + '/3954/numerical_confocal/results/simulation/1M_colony_ca/2D/nodeA_dele_newCN/2Dfinal_%s.pkl'%filename, 'wb'))
    #     pickle.dump(U_record,open(modelling_ephemeral + '/3954/numerical_confocal/results/simulation/1M_colony_ca/2D/nodeA_dele_newCN/2Dtimeseries_%s.pkl'%filename, 'wb'))
    #     print(filename)

except ValueError:
    print('!!!!!!!!!!!!!')
    print('ValueError --> unstable solution')
    print('!!!!!!!!!!!!!')
    print()

    pass




































