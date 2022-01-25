import sys
import os
import sys
print('fg')
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
from adi_v1_openclosed_ca_function import *
from plotting_numerical import *
import pickle


#system parameters
circuit_n=2
variant='5716gaussian'
shape='ca'
mechanism = 'nodeAdele'
parID=int(sys.argv[1])
save_figure=False
L=int(sys.argv[2]); x_gridpoints = int(sys.argv[3]); J = L*x_gridpoints;  I=J
T =int(sys.argv[4]); t_gridpoints =int(sys.argv[5])  ; N = T*t_gridpoints# L,T,J,N = [8,300,80,59700]
#load simulation file
filename = 'circuit%r_variant%s_%s_%sID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant, shape,mechanism,parID,L,J,T,N)#,p_division,kce)
U_record = pickle.load( open( modelling_ephemeral + '/3954/numerical_confocal/results/simulation/1M_colony_%s/2D/full_circuit_newCN/2Dtimeseries_%s.pkl'%(shape,filename), "rb" ) )
U_final = pickle.load( open( modelling_ephemeral + '/3954/numerical_confocal/results/simulation/1M_colony_%s/2D/full_circuit_newCN/2Dfinal_%s.pkl'%(shape,filename), "rb" ) )
print('g')
#
plot_redgreen_contrast(U_final,L,mechanism,shape,filename,'',parID=parID,scale_factor=x_gridpoints,save_figure=save_figure)


rgb_timeseries = redgreen_contrast_timeseries(U_record)
show_rgbvideo(rgb_timeseries,parID)
