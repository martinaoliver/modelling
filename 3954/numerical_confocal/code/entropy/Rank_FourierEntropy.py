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

# def plot(parID,filename,results_path,L=10,mechanism='general',shape='ca',savefig_path='',x_gridpoints=8,save_figure=False):
# #     filename = '2Dfinal_circuit2_variant0_boundary1_ca_generalID%r_L8_J80_T120_N23880.pkl'%parID
#     final_concentration = pickle.load( open( results_path + '/' + filename, "rb" ) )
#     plot_redgreen_contrast(final_concentration,L,mechanism,shape,filename,savefig_path,parID=parID,scale_factor=x_gridpoints,save_figure=save_figure)
%matplotlib inline


def entropy_fourier(final_concentration):

    pixels=final_concentration[0]
    pixels = np.round(pixels,4)
    # plt.imshow(pixels)
    # plt.colorbar()
    # plt.show()
    pixels[pixels == 0] = 0.0000000001

    #normalise so total sum=1
    pixels/=np.sum(pixels)

    H=0
    binwidth = 0.001

    f = fft(pixels)
    f_real = np.real(f)
    f_im= np.imag(f)

    pixels_vector = np.reshape(pixels,-1)
    # hist_i = plt.hist(pixels_vector, bins=np.linspace(np.amin(pixels_vector), np.amax(pixels_vector) + binwidth, 15))
    # plt.show()

    f_real_vector = np.reshape(f_real,-1)
    # hist_a = plt.hist(f_real_vector, bins=np.linspace(np.amin(f_real_vector), np.amax(f_real_vector) + binwidth, 15))[0]
    # plt.show()

    f_im_vector = np.reshape(f_im,-1)
    # hist_b = plt.hist(f_im_vector, bins=np.linspace(np.amin(f_im_vector), np.amax(f_im_vector) + binwidth, 15))[0]

    # plt.show()
    # plt.close()
    def softmax(x):
    # """Compute softmax values for each sets of scores in x."""
        e_x = np.exp(x - np.max(x))
        return e_x / e_x.sum()
    f_real_vector = softmax(f_real_vector)
    f_im_vector = softmax(f_im_vector)

    HKS= np.sum([-p*np.log2(p) for p in pixels_vector if p!=0])
    IKS_real= np.sum([-p*np.log2(p) for p in f_real_vector if p!=0])
    IKS_im= np.sum([-p*np.log2(p) for p in f_im_vector if p!=0])
    print(HKS, IKS_real, IKS_im)
    kSI = HKS - (IKS_real+IKS_im)
    return kSI



folder = 'fullcircuit_5716gaussian'
var=0.23
circuit_n=2
variant='5716gaussian'
shape='square'
mechanism = 'fullcircuit'
L=5; x_gridpoints =10; J = L*x_gridpoints
T =150; t_gridpoints = 100; N = T*t_gridpoints


data_path = modelling_ephemeral + '/3954/numerical_confocal/results/simulation/square/%s/var%r'%(folder,var)
parID_list = pickle.load( open(data_path + '/parID_list_L5_J50_T150_N15000.pkl', "rb" ) )

parID_entropy ={}

plot=True

for parID in tqdm(parID_list[:10], disable=True):


    filename = 'circuit%r_variant%svar%r_%s_%sID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant,var, shape,mechanism,int(parID),L,J,T,N)
    final_concentration = pickle.load( open(data_path + '/2Dfinal_%s.pkl'%filename, 'rb' ) )
    final_concentration = np.round(final_concentration,4)
    # # i, a,b = entropy_fourier(final_concentration)
    kSI = entropy_fourier(final_concentration)
    # print(np.shape(i))
    parID_entropy[parID]=kSI
    # parID_entropy[parID]=entropy(i)
    if plot==True:
        plt.imshow(final_concentration[5])
        plt.colorbar()
        plt.show()
print(parID_entropy)
filename = 'circuit%r_variant%svar%r_%s_%s_L%r_J%r_T%r_N%r'%(circuit_n,variant,var, shape,mechanism,L,J,T,N)

# pickle.dump( parID_entropy, open( modelling_home + "/3954/numerical_confocal/results/entropy/EntropyDicts/HKSdict_%s.pkl"%filename, "wb" ) )
if root == '/rds/general':
    sendemail()
