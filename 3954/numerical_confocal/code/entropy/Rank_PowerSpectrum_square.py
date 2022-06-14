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
# from sendmail import *
from tqdm import tqdm

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from scipy.fft import fft, ifft
import matplotlib as mpl
from PIL import Image, ImageDraw
import numpy as np
from numpy import asarray
import pickle
from tqdm import tqdm
#############################



# %matplotlib inline
def fourier2d(final_concentration):
    # print('Normalised PowerSpectrum of original')
    pixels=final_concentration[5]
    pixels[pixels == 0] = 0.0000000001
    pixels/=np.sum(pixels) #line added



    #fourier transform
    f = np.fft.ifftshift(pixels)
    f = np.fft.fft2(f)
    f = np.fft.fftshift(f)
    max_coordinates=list(zip(np.where(f == np.amax(f))[0], np.where(f == np.amax(f))[1]))[0]
    PowerSpectrum = np.real(f)**2 +  np.imag(f)**2
    PowerSpectrum[max_coordinates]=0
    if np.sum(PowerSpectrum)!=0:
        PowerSpectrum/=np.sum(PowerSpectrum)

    # print('sum',np.sum(PowerSpectrum))

    # print('max', np.amax(PowerSpectrum), 'min', np.amin(PowerSpectrum))

    # plt.subplot(132)
    # plt.imshow(PowerSpectrum)#,norm=LogNorm())
    # plt.colorbar()


    # binwidth = 0.00001
    PowerSpectrum_vector = np.reshape(PowerSpectrum,-1)
    # print(np.sum(PowerSpectrum_vector))
    # plt.subplot(133)
    # hist_PowerSpectrum= plt.hist(PowerSpectrum_vector, bins=np.linspace(np.amin(PowerSpectrum_vector), np.amax(PowerSpectrum_vector) + binwidth, 300))
    # plt.yscale('log')
    # plt.show()



    return PowerSpectrum_vector

def entropy(histogram):
    x= np.sum([-p*np.log2(p) for p in histogram if p!=0])
    return x

lhs=True

if lhs==True:
    folder = 'fullcircuit_5716gaussian/1M_turingI'
    variant=0
else:
    var=float(sys.argv[1])
    # var=0.23
    folder = 'fullcircuit_5716gaussian/var%s'%var
    variant='5716gaussian'

circuit_n=2
shape='square'
mechanism = 'fullcircuit'

L=5; x_gridpoints =10; J = L*x_gridpoints
T =2000; t_gridpoints = 10; N = T*t_gridpoints

data_path = modelling_home + '/3954/numerical_confocal/results/simulation/square/%s'%(folder)
parID_list = pickle.load( open(data_path + '/parID_list_L%r_J%r_T%r_N%r.pkl'%(L,J,T,N), "rb" ) )

parID_ps = {}

plot=False
# parID_list=[497]
for parID in tqdm(parID_list, disable=False):
# for parID in tqdm(parID_list[:20], disable=False):
# for parID in tqdm([805,686,472,252,688], disable=True):
    # print('parID',parID)
    if lhs==True:
        filename = 'circuit%r_variant%s_%s_%sID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant, shape,mechanism,int(parID),L,J,T,N)
    else:
        filename = 'circuit%r_variant%svar%r_%s_%sID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant,var, shape,mechanism,int(parID),L,J,T,N)

    final_concentration0 = pickle.load( open(data_path + '/2Dfinal_%s.pkl'%filename, 'rb' ) )
    if not np.isnan(final_concentration0).any():
        final_concentration = np.round(final_concentration0,4)
        if plot==True:
            plt.figure(figsize=[14,2])
            plt.subplot(131)
            plt.imshow(final_concentration[5])
            plt.colorbar()
        hist_PowerSpectrum  = fourier2d(final_concentration)
        ps = entropy(hist_PowerSpectrum)

        parID_ps[parID] = ps
    # print(final_concentration0)


if lhs==True:
    filename = 'circuit%r_variant%s_%s_%s_L%r_J%r_T%r_N%r'%(circuit_n,variant, shape,mechanism,L,J,T,N)
else:
    filename = 'circuit%r_variant%svar%s_%s_%s_L%r_J%r_T%r_N%r'%(circuit_n,variant,var, shape,mechanism,L,J,T,N)

pickle.dump( parID_ps, open( modelling_home + "/3954/numerical_confocal/results/entropy/EntropyDicts/psprenorm_dict_%s.pkl"%filename, "wb" ) )


# if root=='/rds/general/user/mo2016':
#     sendmail('general_filename')