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



from plotting_numerical import plot_redgreen_contrast
from tqdm import tqdm

import matplotlib.pyplot as plt
from scipy.fft import fft, ifft
import matplotlib as mpl
from PIL import Image, ImageDraw
import numpy as np
from numpy import asarray
import pickle
from tqdm import tqdm
# from send_email import *
from scipy.ndimage import laplace
#############################
# %matplotlib inline 



def compute_LaplaceSum(M1, plot=True):
    l = laplace(M1)
    l_abs = np.absolute(l)
    sum_l = np.sum(l_abs)
    if plot==True:
        plt.imshow(l_abs)
        plt.colorbar()
        plt.show()

    
    return sum_l


folder = 'fullcircuit_5716gaussian/1M_turingI'
# var=0.23
circuit_n=2
variant=0#'5716gaussian'
shape='square'
mechanism = 'fullcircuit'
# L=5; x_gridpoints =10; J = L*x_gridpoints
# T =150; t_gridpoints = 100; N = T*t_gridpoints
L=5; x_gridpoints =10; J = L*x_gridpoints
T =2000; t_gridpoints = 10; N = T*t_gridpoints

data_path = modelling_home + '/3954/numerical_confocal/results/simulation/square/%s'%(folder)
parID_list = pickle.load( open(data_path + '/parID_list_L5_J50_T2000_N20000.pkl', "rb" ) )


plot=False
parID_lpl = {}
for parID in tqdm(parID_list, disable=False):
# for parID in tqdm([805,686,472,252,688], disable=True):
    # parID=435376
    # print(parID)

    # filename = 'circuit%r_variant%svar%r_%s_%sID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant,var, shape,mechanism,int(parID),L,J,T,N)
    filename = 'circuit%r_variant%s_%s_%sID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant, shape,mechanism,int(parID),L,J,T,N)
    final_concentration = pickle.load( open(data_path + '/2Dfinal_%s.pkl'%filename, 'rb' ) )
    final_concentration0 = pickle.load( open(data_path + '/2Dfinal_%s.pkl'%filename, 'rb' ) )
    if not np.isnan(final_concentration0).any():
        final_concentration = np.round(final_concentration0,4)    
        lpl_sum = compute_LaplaceSum(final_concentration[5], plot=plot)
    # print('kSI',kSI, 'HKS',HKS, 'IKS_real',IKS_real, 'IKS_im',IKS_im)
    # print(lpl_sum)
        parID_lpl[parID] = lpl_sum

    # if plot==True:

    #     plt.imshow(final_concentration[5])
    #     plt.colorbar()
    #     plt.show()
# print(parID_entropy)
# filename = 'circuit%r_variant%svar%r_%s_%s_L%r_J%r_T%r_N%r'%(circuit_n,variant,var, shape,mechanism,L,J,T,N)
filename = 'circuit%r_variant%s_%s_%s_L%r_J%r_T%r_N%r'%(circuit_n,variant, shape,mechanism,L,J,T,N)

pickle.dump( parID_lpl, open( modelling_home + "/3954/numerical_confocal/results/entropy/EntropyDicts/lpl_dict_%s.pkl"%filename, "wb" ) )
