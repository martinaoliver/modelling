#############################
#IMPORTS#
#############################
import sys
import os

from numpy import searchsorted
pwd = os.getcwd()
root = pwd.rpartition("mo2016")[0] + pwd.rpartition("mo2016")[1] #/Volumes/mo2016/ or '/Users/mo2016/' or '/rds/general/mo2016/'
if root == '/Users/mo2016':
    print(root)
    modelling_ephemeral = '/Volumes/mo2016/ephemeral/Documents/modelling'
    modelling_home = '/Volumes/mo2016/home/Documents/modelling'
    modelling_local = root + '/Documents/modelling'


if root == '/Volumes/mo2016' or root=='/rds/general/user/mo2016': #'/rds/general' or root=='/Volumes':
        modelling_ephemeral = root + '/ephemeral/Documents/modelling'
        modelling_home = root  + '/home/Documents/modelling'
        modelling_local = modelling_home

if root == '/Users/mo2016' or  root == '/Volumes/mo2016':
    import matplotlib as mpl
    mpl.use('tkagg')

modulepath = modelling_local + '/3954/modules/new_CN'

sys.path.append(modulepath)



from adi_square_function import *
from plotting_numerical import *
import pickle
#execution parameters

mechanism = 'fullcircuit'
shape = 'square' #'ca'
parID = int(sys.argv[1])
circuit_n=2
variant='0'
folder='fullcircuit2var0'
parametersets_n = 1000000
save_figure = False
tqdm_disable = False #disable tqdm
n_species=6
# var=0.001
# open parameter dictionaries
# general_df= pickle.load( open(modelling_home + '/3954/parameter_space_search/parameterfiles/5716gaussian/df_circuit%r_variant%s_%rparametersets_%rvar.pkl'%(circuit_n,variant,parametersets_n), "rb" ) )
general_df= pickle.load( open(modelling_home + '/3954/parameter_space_search/results/turing_dataframes/turing_lsa_df_circuit2_variant0_1000000parametersets.pkl', "rb" ) )
# par_dict = general_df.loc[parID].to_dict()
par_dict = general_df.iloc[parID].to_dict()
print(par_dict)
d_A = par_dict['d_A']
d_B = par_dict['d_B']
D = np.zeros(n_species)
D[0]=d_A
# d_B=1
D[1]=d_B
steadystates=par_dict['ss_list']
# var_list=[0.01,0.02,0.04, 0.06,0.08, 0.1,0.23]
#solver parameters
L_x=int(sys.argv[2]); x_gridpoints = int(sys.argv[3]); L=L_x; J = L*x_gridpoints;  L_y=L_x; I=J
T =int(sys.argv[4]); t_gridpoints = int(sys.argv[5]) ; N = T*t_gridpoints
# suggested_tgridpoints = x_gridpoints**2

# filename = 'circuit%r_variant%svar%r_%s_%sID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant,var, shape,mechanism,parID,L,J,T,N)
filename = 'circuit%r_variant%s_%s_%sID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant, shape,mechanism,parID,L,J,T,N)
savefig_path = modelling_ephemeral + '/3954/numerical_confocal/results/figures/square/%s'%(folder)
# savefig_path = modelling_ephemeral + '/3954/numerical_confocal/results/figures/square/%s/var%r'%(folder,var)
# try:
U_record,U_final = adi(par_dict,L_x,L_y,J,I,T,N, circuit_n,n_species,D, tqdm_disable=tqdm_disable, steadystates=steadystates)#,p_division=p_division,seed=seed)
plot_2D_final_concentration(U_final,L_x,J,filename,savefig_path,n_species=n_species,save_figure=False)
print('asfdg')
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

# except ValueError:
#     print('!!!!!!!!!!!!!')
#     print('ValueError --> unstable solution')
#     print('!!!!!!!!!!!!!')
#     print()

#     pass
