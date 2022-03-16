#############################
#IMPORTS#
#############################
import sys
import os
pwd = os.getcwd()
root = pwd.rpartition("mo2016")[0] + pwd.rpartition("mo2016")[1] #/Volumes/mo2016/ or '/Users/mo2016/' or '/rds/general/mo2016/'
print(root)
if root == '/Users/mo2016':
    modelling_ephemeral = '/Volumes/mo2016/ephemeral/Documents/modelling'
    modelling_home = '/Volumes/mo2016/home/Documents/modelling'
    modelling_local = root + '/Documents/modelling'
    import matplotlib as mpl
    mpl.use('tkagg')

if root == '/Volumes/mo2016' or root=='/rds/general/user/mo2016': #'/rds/general' or root=='/Volumes':
        modelling_ephemeral = root + '/ephemeral/Documents/modelling'
        modelling_home = root  + '/home/Documents/modelling'
        modelling_local = modelling_home

modulepath = modelling_local + '/3954/modules/new_CN'

sys.path.append(modulepath)


from adi_function import *
from plotting_numerical import *
import pickle
#execution parameters

mechanism = 'fullcircuit'
shape = 'square' #'ca'
parID = int(sys.argv[1])
circuit_n=2
variant='5716gaussian'
folder='fullcircuit_5716gaussian'
parametersets_n = 1000 #1000000
save_figure = True
tqdm_disable = False #disable tqdm
n_species=6
var=0.001
# open parameter dictionaries
general_df= pickle.load( open(modelling_home + '/3954/parameter_space_search/parameterfiles/5716gaussian/df_circuit%r_variant%s_%rparametersets_%rvar.pkl'%(circuit_n,variant,parametersets_n,var), "rb" ) )
# general_df= pickle.load( open(modelling_home + '/3954/parameter_space_search/results/output_dataframes/lsa_df_circuit2_variant0_20000parametersets.pkl', "rb" ) )
par_dict = general_df.loc[parID].to_dict()

d_A = par_dict['d_A']
d_B = par_dict['d_B']
D = np.zeros(n_species)
D[0]=d_A
d_B=1
D[1]=d_B


#solver parameters
L_x=int(sys.argv[2]); x_gridpoints = int(sys.argv[3]); J = L_x*x_gridpoints;  L_y=L_x; I=J
T =int(sys.argv[4]); t_gridpoints = int(sys.argv[5]) ; N = T*t_gridpoints
suggested_tgridpoints = x_gridpoints**2

filename = 'circuit%r_variant%s_%s_%sID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant, shape,mechanism,parID,L_x,J,T,N)#,p_division,kce)
savefig_path = modelling_ephemeral + '/3954/numerical_confocal/results/figures/square/%s'%folder
try:
    U_record,U_final = adi(par_dict,L_x,L_y,J,I,T,N, circuit_n,n_species,D, tqdm_disable=tqdm_disable)#,p_division=p_division,seed=seed)
    plot_2D_final_concentration(U_final,L_x,J,filename,savefig_path,n_species=n_species,save_figure=True)
        # savefig_path = modelling_ephemeral + '/3954/numerical_confocal/results/figures/1M_colony_ca/2D/5716gaussian'
    # plot_redgreen_contrast(U_final,L_x,mechanism,shape,filename,savefig_path,parID=parID,scale_factor=x_gridpoints,save_figure=save_figure)
    # plot_redgreen_contrast(U_final,L_x,mechanism,shape,filename,savefig_path,parID=parID,scale_factor=x_gridpoints,save_figure=False)
    #
    # rgb_timeseries = redgreen_contrast_timeseries(U_record)
    # show_rgbvideo(rgb_timeseries,parID)
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
