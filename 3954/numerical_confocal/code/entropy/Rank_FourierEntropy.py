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
from tqdm import tqdm
# def plot(parID,filename,results_path,L=10,mechanism='general',shape='ca',savefig_path='',x_gridpoints=8,save_figure=False):
# #     filename = '2Dfinal_circuit2_variant0_boundary1_ca_generalID%r_L8_J80_T120_N23880.pkl'%parID
#     final_concentration = pickle.load( open( results_path + '/' + filename, "rb" ) )
#     plot_redgreen_contrast(final_concentration,L,mechanism,shape,filename,savefig_path,parID=parID,scale_factor=x_gridpoints,save_figure=save_figure)

def entropy_fourier(final_concentration, channel):

    if channel == 'red':
        pixels=final_concentration[-2]
    if channel == 'green':
            pixels=final_concentration[-1]

    pixels[pixels == 0] = 0.0000000001

    #normalise so total sum=1
    pixels/=np.sum(pixels)

    H=0
    binwidth = 0.001

    f = fft(pixels)
    f_real = np.real(f)
    f_im= np.imag(f)

    pixels_vector = np.reshape(pixels,-1)
    hist_i = plt.hist(pixels_vector, bins=np.linspace(np.amin(pixels_vector), np.amax(pixels_vector) + binwidth, 15))

    f_real_vector = np.reshape(f_real,-1)
    hist_a = plt.hist(f_real_vector, bins=np.linspace(np.amin(f_real_vector), np.amax(f_real_vector) + binwidth, 15))

    f_im_vector = np.reshape(f_im,-1)
    hist_b = plt.hist(f_im_vector, bins=np.linspace(np.amin(f_im_vector), np.amax(f_im_vector) + binwidth, 15))


    plt.close()
    return hist_i, hist_a,hist_b

def entropy(histogram):
    return np.sum([p*np.log2(p) for p in histogram[0] if p!=0])

def kSI(i,a,b): #entropy - (fourier entropy re + fourier entropy imag)
    HKS = entropy(i)
    IKS_real = entropy(a)
    IKS_im = entropy(b)
    kSI = HKS - (IKS_real + IKS_im) #v1
    return(kSI)

results_path = modelling_home + '/3954/numerical_confocal/results/simulation/1M_colony_ca/2D/full_circuit_newCN/'
parID_list = pickle.load( open( results_path + '/' + 'parID_list_5716gaussian_L10J150T120N1200.pkl', "rb" ) )
parID_entropy_dict_red = {}
parID_entropy_dict_green = {}
parID_entropy_dict_redgreen= {}

for parID in tqdm(parID_list):

    filename = '2Dfinal_circuit2_variant5716gaussian_ca_nodeAdeleID%s_L10_J150_T120_N1200.pkl'%parID
    final_concentration = pickle.load( open( results_path + '/' + filename, "rb" ) )

    i, a,b = entropy_fourier(final_concentration,'red')
    kSI_red = kSI(i,a,b)
    parID_entropy_dict_red[parID]=kSI_red

    i, a,b = entropy_fourier(final_concentration,'green')
    kSI_green = kSI(i,a,b)
    parID_entropy_dict_green[parID]=kSI_green

    parID_entropy_dict_redgreen[parID]=kSI_red + kSI_green

parID_entropy_dict_redgreen = dict(sorted(parID_entropy_dict_redgreen.items(), key=lambda item: item[1]))
parID_entropy_dict_red= dict(sorted(parID_entropy_dict_red.items(), key=lambda item: item[1]))
parID_entropy_dict_green = dict(sorted(parID_entropy_dict_green.items(), key=lambda item: item[1]))

print(parID_entropy_dict_redgreen)
print(parID_entropy_dict_green)
print(parID_entropy_dict_red)

pickle.dump( parID_entropy_dict_redgreen, open( "EntropyDicts/EntropyDictdictRedGreen_v1_5716gaussian_L10J150T120N1200.pkl", "wb" ) )
pickle.dump( parID_entropy_dict_red, open( "EntropyDicts/EntropyDictdictRed_v1_5716gaussian_L10J150T120N1200.pkl", "wb" ) )
pickle.dump( parID_entropy_dict_green, open( "EntropyDicts/EntropyDictdictGreen_v1_5716gaussian_L10J150T120N1200.pkl", "wb" ) )
