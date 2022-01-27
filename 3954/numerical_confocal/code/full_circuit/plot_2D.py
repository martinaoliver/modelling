#############################
#IMPORTS#
#############################
import sys
import os
pwd = os.getcwd()
root = pwd.rpartition("mo2016")[0] + pwd.rpartition("mo2016")[1] #/Volumes/mo2016/ or '/Users/mo2016/' or '/rds/general/mo2016/'
if root == '/Users/mo2016':
    print('fasdghgfsgth')
    modelling_ephemeral = '/Volumes/mo2016/ephemeral/Documents/modelling'
    modelling_home = '/Volumes/mo2016/home/Documents/modelling'
else:
    modelling_ephemeral = root + 'ephemeral/Documents/modelling'
    modelling_home = root  + '/Documents/modelling'

modelling_path_local = root + '/Documents/modelling'

modulepath = modelling_path_local + '/3954/modules/new_CN'
sys.path.append(modulepath)
from plotting_numerical import plot_redgreen_contrast
from tqdm import tqdm
import matplotlib as mpl
mpl.use('tkagg')
import pickle


#execution parameters
circuit_n=2
variant='5716gaussian'
parametersets_n = 1000
save_figure = False
tqdm_disable = False #disable tqdm
n_species=6
# open parameter dictionaries
general_df = pickle.load(open(modelling_home + '/3954/parameter_space_search/results/output_dataframes/lsa_df_circuit%r_variant%r_%rparametersets.pkl'%(2,variant,parametersets_n), "rb"))
#chose parameter sets to analyse
parID = int(sys.argv[1])
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

    L,J,T,N = [int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])]
    t_gridpoints = int(N/T)
    x_gridpoints = int(J/L)

    initial_condition = [0.001]*n_species
    filename = '2Dfinal_circuit%r_variant%s_%s_%sID%s_L%r_J%r_T%r_N%r'%(circuit_n,variant,shape,mechanism, parID,L,J,T,N)
    final_concentration = pickle.load( open(modelling_home + '/3954/numerical_confocal/results/simulation/1M_colony_ca/2D/full_circuit_newCN/%s.pkl'%filename, 'rb' ) )
    save_path =modelling_home + '/3954/numerical_confocal/results/simulation/1M_colony_ca/2D/full_circuit_newCN'
    rgb = plot_redgreen_contrast(final_concentration,L,mechanism,shape,filename,save_path,parID=parID,scale_factor=x_gridpoints,save_figure=save_figure)
    print(rgb)
