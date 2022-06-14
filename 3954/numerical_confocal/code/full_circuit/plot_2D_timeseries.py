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


from adi_function import *
from plotting_numerical import *
import pickle

#system parameters

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
save_figure=False
#solver parameters
L_x=int(sys.argv[2]); x_gridpoints = int(sys.argv[3]); J = L_x*x_gridpoints;  L_y=L_x; I=J
T =int(sys.argv[4]); t_gridpoints = int(sys.argv[5]) ; N = T*t_gridpoints

#load simulation file
filename = 'circuit%r_variant%s_%s_%sID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant, shape,mechanism,parID,L_x,J,T,N)#,p_division,kce)
savefig_path = modelling_ephemeral + '/3954/numerical_confocal/results/figures/square/%s'%folder
U_record = pickle.load( open( modelling_ephemeral + '/3954/numerical_confocal/results/simulation/square/%s/2Dtimeseries_%s.pkl'%(folder,filename), "rb" ) )
print('g')
#
# plot_redgreen_contrast(U_final,L,mechanism,shape,filename,'',parID=parID,scale_factor=x_gridpoints,save_figure=save_figure)


rgb_timeseries = redgreen_contrast_timeseries(U_record)
show_rgbvideo(rgb_timeseries,parID)
