from scipy.fft import fft, ifft
import matplotlib as mpl
mpl.use('tkagg')
import sys
import os
pwd = os.getcwd()
root = pwd.split("home", 1)[0]
modelling_home = root + 'home/Documents/modelling'
modelling_ephemeral = root + 'ephemeral/Documents/modelling'
modulepath = modelling_home + '/3954/modules'
sys.path.append(modulepath)
from numerical_solvers_variableboundary import *
from PIL import Image, ImageDraw
import numpy as np
from numpy import asarray
import pickle
#
# def plot(parID,filename,results_path,L=10,mechanism='general',shape='ca',savefig_path='',x_gridpoints=8,save_figure=False):
# #     filename = '2Dfinal_circuit2_variant0_boundary1_ca_generalID%r_L8_J80_T120_N23880.pkl'%parID
#     final_concentration = pickle.load( open( results_path + '/' + filename, "rb" ) )
#     plot_redgreen_contrast(final_concentration,L,mechanism,shape,filename,savefig_path,parID=parID,scale_factor=x_gridpoints,save_figure=save_figure)

def entropy_fourier(final_concentration):

    pixels_red,pixels_green=final_concentration[-2],final_concentration[-1]
    pixels_red[pixels_red == 0] = 0.0000000001
    pixels_green[pixels_green == 0] = 0.0000000001

    #normalise so total sum=1
    pixels_red/=np.sum(pixels_red)
    pixels_green/=np.sum(pixels_green)

    H_red=0
    H_green=0
    binwidth_real = 0.001
    binwidth_im = 0.001

    f_red = fft(pixels_red)
    f_red_real = np.real(f_red)
    f_red_im = np.imag(f_red)

    pixels_red_vector = np.reshape(pixels_red,-1)
    hist_i = plt.hist(pixels_red_vector, bins=np.linspace(np.amin(pixels_red_vector), np.amax(pixels_red_vector) + binwidth_real, 15))

    f_red_real_vector = np.reshape(f_red_real,-1)
    hist_a = plt.hist(f_red_real_vector, bins=np.linspace(np.amin(f_red_real_vector), np.amax(f_red_real_vector) + binwidth_real, 15))

    f_red_im_vector = np.reshape(f_red_im,-1)
    hist_b = plt.hist(f_red_im_vector, bins=np.linspace(np.amin(f_red_im_vector), np.amax(f_red_im_vector) + binwidth_im, 15))


    plt.close()
    return hist_i, hist_a,hist_b

def entropy(histogram):
    return np.sum([p*np.log2(p) for p in histogram[0] if p!=0])

def kSI(HSK,IKS_real,IKS_imag): #entropy - (fourier entropy re + fourier entropy imag)
    kSI = HKS - (IKS_real + IKS_im)
    return(kSI)

results_path = modelling_home + '/3954/numerical_confocal/results/simulation/1M_colony_ca/2D/full_circuit_newCN/'
parID_list = pickle.load( open( results_path + '/' + 'parID_list_5716gaussian_L10J150T120N1200.pkl', "rb" ) )
parID_entropy_dict = {}
for parID in parID_list[:10]:

    filename = '2Dfinal_circuit2_variant5716gaussian_ca_nodeAdeleID%s_L10_J150_T120_N1200.pkl'%parID
    final_concentration = pickle.load( open( results_path + '/' + filename, "rb" ) )
    i, a,b = entropy_fourier(final_concentration)


    HKS = entropy(i)
    IKS_real = entropy(a)
    IKS_im = entropy(b)





    parID_entropy_dict[parID]=kSI(HKS,IKS_real,IKS_im)
print(parID_entropy_dict)
parID_entropy_dict = dict(sorted(parID_entropy_dict.items(), key=lambda item: item[1]))
print(parID_entropy_dict)
