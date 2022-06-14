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
from send_email import *
#############################

#input is final concentration, output is entropy of image
def entropy(final_concentration):

    pixels=final_concentration[5]
    pixels = np.round(pixels,4)
    pixels[pixels == 0] = 0.0000000001

    #normalise so total sum=1
    pixels/=np.sum(pixels)

    pixels_vector = np.reshape(pixels,-1)
    HKS= np.sum([-p*np.log2(p) for p in pixels_vector if p!=0])
    return HKS


#definition of variables and path
circuit_n=2
variant=0#'5716gaussian'
shape='ca'
mechanism = 'fullcircuit'#nodeAdele' #not really nodeAdele, mistake in filename
L=10; x_gridpoints =15; J = L*x_gridpoints
T =120; t_gridpoints = 10; N = T*t_gridpoints
data_path = modelling_home + '/3954/numerical_confocal/results/simulation/1M_colony_ca/2D/full_circuit_newCN'
# parID_list = pickle.load( open(data_path + '/parID_list_5716gaussian_L10J150T120N1200.pkl', "rb" ) )

parID_list = pickle.load( open(data_path + '/parID_list_variant0_ca_fullcircuit_L10J150T120N1200.pkl', "rb" ) )

print(parID_list)

parID_entropy ={}
plot=False
for parID in tqdm(parID_list, disable=False):

    #open concentration file

    filename = 'circuit%r_variant%s_%s_%sID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant, shape,mechanism,int(parID),L,J,T,N)
    final_concentration = pickle.load( open(data_path + '/2Dfinal_%s.pkl'%filename, 'rb' ) )
    final_concentration = np.round(final_concentration,4)
    
    #calculate entropy and save to dictionary
    HKS = entropy(final_concentration)
    parID_entropy[parID]=HKS
    
    #if informative, plot image
    if plot==True:
        plt.imshow(final_concentration[5])
        plt.colorbar()
        plt.show()

#save dictionary
filename = 'circuit%r_variant%s_%s_%s_L%r_J%r_T%r_N%r_test'%(circuit_n,variant, shape,mechanism,L,J,T,N)
pickle.dump( parID_entropy, open( modelling_home + "/3954/numerical_confocal/results/entropy/EntropyDicts/HKSdict_%s.pkl"%filename, "wb" ) )
if root == '/rds/general':
    sendemail()


