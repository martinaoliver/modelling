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



from numerical.adi_ca_function_openclosed_nodilution import adi_ca_openclosed_nodilution
from numerical.adi_ca_function_openclosed_nodilution_numba import adi_ca_openclosed_nodilution as adi_ca_openclosed_nodilution_numba
from numerical.plotting_numerical import *
import pickle
#execution parameters

mechanism = 'fullcircuit'
shape = 'ca'
# parID = int(sys.argv[1])
parID = int(1)

circuit_n=2
variant= 9
folder = 'fullcircuit/1M_turingI'#'fullcircuit/1M'#'fullcircuit/1M_turingI'
n_species = 6

parametersets_n = 10 #1000000
save_figure = False
tqdm_disable = False #disable tqdm
# boundarycoeff = float(sys.argv[6])
boundarycoeff = 1.5
seed=1;p_division=0.41#0.147#0.5
# open parameter dictionaries
general_df= pickle.load( open(modellingpath + '/3954/parameter_space_search/parameterfiles/df_circuit%r_variant%s_%rparametersets.pkl'%(circuit_n,variant,parametersets_n), "rb" ) )
# general_df= pickle.load( open(modellingpath + "/3954/parameter_space_search//results/turing_dataframes/turing_lsa_df_circuit%r_variant%r_%rparametersets.pkl"%(circuit_n,variant,10), "rb" ) )
par_dict = general_df.loc[parID].to_dict()

d_A = par_dict['d_A']
d_B = par_dict['d_B']

D = np.zeros(n_species)
D[0]=d_A
d_B=1
D[1]=d_B
print(par_dict)
par_dict['mulva'] = par_dict['mulva'] + np.log(2)*p_division


#solver parameters
# L_x=int(sys.argv[2]); x_gridpoints = int(sys.argv[3]); L=L_x; J = L*x_gridpoints;  L_y=L_x; I=J
# T =int(sys.argv[4]); t_gridpoints = int(sys.argv[5]) ; N = T*t_gridpoints

L_x=int(8); x_gridpoints = int(20); L=L_x; J = L*x_gridpoints;  L_y=L_x; I=J
T =int(1/0.05); t_gridpoints = int(150) ; N = T*t_gridpoints

suggested_tgridpoints = x_gridpoints**2

# filename = 'circuit%r_variant%s_%s_%sID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant, shape,mechanism,parID,L,J,T,N)
filename = 'circuit%r_variant%s_bc%s_%s_%sID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundarycoeff, shape,mechanism,parID,L,J,T,N)
savefig_path = 'test_results'
# U_record,U_final = adi_ca_openclosed_nodilution_numba(par_dict,L_x,L_y,J,I,T,N, circuit_n,n_species,D, seed=1, p_division=p_division, tqdm_disable=tqdm_disable, growth='Slow')#,p_division=p_division,seed=seed)
U_record,U_final = adi_ca_openclosed_nodilution(par_dict,L_x,L_y,J,I,T,N, circuit_n,n_species,D, seed=1, p_division=p_division, tqdm_disable=tqdm_disable, growth='Slow')#,p_division=p_division,seed=seed)


plot_redgreen_contrast(U_final,L_x,mechanism,shape,filename,savefig_path,parID=parID,scale_factor=x_gridpoints,save_figure=save_figure)

# plt.imshow(U_final[-1])
# plt.show()
# plt.imshow(U_final[0])
# plt.show()
rgb_timeseries = redgreen_contrast_timeseries(U_record)
show_rgbvideo(rgb_timeseries,parID)



































