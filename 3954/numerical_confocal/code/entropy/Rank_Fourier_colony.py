#############################
#IMPORTS#
#############################
import sys
import os
pwd = os.getcwd()
root = pwd.rpartition("mo2016")[0] + pwd.rpartition("mo2016")[1] #/Volumes/mo2016/ or '/Users/mo2016/' or '/rds/general/mo2016/'


if root == '/Users/mo2016':
    modelling_ephemeral = '/Volumes/mo2016/ephemeral/Documents/modelling'
    modelling_home = '/Volumes/mo2016/home/Documents/modelling'
    modelling_local = root + '/Documents/modelling'
    import matplotlib as mpl
    mpl.use('tkagg')

if root == '/Volumes/mo2016' or root=='/rds/general': #'/rds/general' or root=='/Volumes':
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



# %matplotlib inline
def entropy_fourier(final_concentration):

    pixels=final_concentration[5]
    # plt.imshow(final_concentration[5])
    # plt.colorbar()
    # plt.show()
    pixels[pixels == 0] = 0.0000000001

    #normalise so total sum=1
    pixels/=np.sum(pixels)

    H=0
    binwidth = 0.00001

    f = fft(pixels)
    f_real = np.real(f)
    f_im= np.imag(f)
    # plt.imshow(f_real)
    # plt.show()
    # plt.imshow(f_im)
    # plt.show()
    pixels_vector = np.reshape(pixels,-1)
    hist_i = plt.hist(pixels_vector, bins=np.linspace(np.amin(pixels_vector), np.amax(pixels_vector) + binwidth, 15))
    # plt.show()

    f_real_vector = np.reshape(f_real,-1)
    hist_a = plt.hist(f_real_vector, bins=np.linspace(np.amin(f_real_vector), np.amax(f_real_vector) + binwidth, 15))[0]
    # plt.show()

    f_im_vector = np.reshape(f_im,-1)
    hist_b = plt.hist(f_im_vector, bins=np.linspace(np.amin(f_im_vector), np.amax(f_im_vector) + binwidth, 15))[0]

    # plt.show()
    def softmax(x):
    # """Compute softmax values for each sets of scores in x."""
        e_x = np.exp(x - np.max(x))
        return e_x / e_x.sum()
    f_real_vector = softmax(f_real_vector)
    f_im_vector = softmax(f_im_vector)

    hist_a/=np.sum(hist_a)
    hist_b/=np.sum(hist_b)
    return pixels_vector, hist_a,hist_b

def entropy(histogram):
    x= np.sum([-p*np.log2(p) for p in histogram if p!=0])
    # print(x)
    return x

def compute_metrics(i,a,b): #entropy - (fourier entropy re + fourier entropy imag)
    HKS = entropy(i)
    IKS_real = entropy(a)
    IKS_im = entropy(b)
    kSI = HKS - (IKS_real + IKS_im) 
    return kSI,HKS, IKS_real, IKS_im


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

parID_kSI ={}
parID_HKS ={}
parID_IKS_real ={}
parID_IKS_im ={}

plot=True

for parID in tqdm(parID_list, disable=False):
# for parID in tqdm([805,686,472,252,688], disable=True):


    filename = 'circuit%r_variant%s_%s_%sID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant, shape,mechanism,int(parID),L,J,T,N)
    final_concentration = pickle.load( open(data_path + '/2Dfinal_%s.pkl'%filename, 'rb' ) )
    final_concentration = np.round(final_concentration,4)
    i, a,b = entropy_fourier(final_concentration)
    # kSI = entropy(final_concentration)
    # print(np.shape(i))
    kSI, HKS, IKS_real, IKS_im = compute_metrics(i,a,b)
    # print('kSI',kSI, 'HKS',HKS, 'IKS_real',IKS_real, 'IKS_im',IKS_im)

    parID_kSI[parID] = kSI
    parID_HKS[parID] = HKS
    parID_IKS_real[parID] = IKS_real 
    parID_IKS_im[parID] = IKS_im

    # parID_entropy[parID]=entropy(i)
    # if plot==True:
    #     plt.imshow(final_concentration[5])
    #     plt.colorbar()
    #     plt.show()
# print(parID_entropy)
filename = 'circuit%r_variant%s_%s_%s_L%r_J%r_T%r_N%r_test'%(circuit_n,variant, shape,mechanism,L,J,T,N)

pickle.dump( parID_kSI, open( modelling_home + "/3954/numerical_confocal/results/entropy/EntropyDicts/kSIdict_%s.pkl"%filename, "wb" ) )
pickle.dump( parID_HKS, open( modelling_home + "/3954/numerical_confocal/results/entropy/EntropyDicts/HKSdict_%s.pkl"%filename, "wb" ) )
pickle.dump( parID_IKS_real, open( modelling_home + "/3954/numerical_confocal/results/entropy/EntropyDicts/IKS_real_dict_%s.pkl"%filename, "wb" ) )
pickle.dump( parID_IKS_im, open( modelling_home + "/3954/numerical_confocal/results/entropy/EntropyDicts/IKS_im_dict_%s.pkl"%filename, "wb" ) )
if root == '/rds/general':
    sendemail()

