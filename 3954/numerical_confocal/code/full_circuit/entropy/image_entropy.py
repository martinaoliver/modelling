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
import matplotlib as mpl
import matplotlib.pyplot as plt
# mpl.use('tkagg')

#######
#FUNCTIONS#
#######
def plot(parID,filename,results_path,L=10,mechanism='general',shape='ca',savefig_path='',x_gridpoints=8,save_figure=False):
    final_concentration = pickle.load( open( results_path + '/' + filename, "rb" ) )
    plot_redgreen_contrast(final_concentration,L,mechanism,shape,filename,savefig_path,parID=parID,scale_factor=x_gridpoints,save_figure=save_figure)

def entropy(parID,filename,results_path,show_fig=False,L=10,mechanism='general',shape='ca',savefig_path='',x_gridpoints=8,save_figure=False):
    final_concentration = pickle.load( open( results_path + '/' + filename, "rb" ) )

    if show_fig==True:
        plot(parID,filename,results_path)

    pixels_red,pixels_green=final_concentration[-2],final_concentration[-1]
    pixels_red[pixels_red == 0] = 0.0000000001
    pixels_green[pixels_green == 0] = 0.0000000001

    pixels_red/=np.sum(pixels_red)
    pixels_green/=np.sum(pixels_green)

    H_red=0
    H_green=0

    for counti,i in enumerate(pixels_red):
        for countj,p in enumerate(i):
            pixelH = p*np.log2(p)
            H_red -= pixelH

    for counti,i in enumerate(pixels_green):
        for countj,p in enumerate(i):
            pixelH = p*np.log2(p)
            H_green -= pixelH
    H = np.mean([H_red,H_green])
    return H_red,H_green,H


#######
#EXECUTE CODE#
#######

# parID = 11428
results_path = modelling_ephemeral + '/3954/numerical_confocal/results/simulation/1M_colony_ca/2D/full_circuit'
# parID_list  = pickle.load( open( results_path + "/parID_list_general_8x10T120.pkl", "rb" ) )
# a = entropy(parID,filename,results_path,show_fig=True)


from tqdm import tqdm
parID_list = pickle.load( open( results_path + '/' + 'parID_list_general_8x10T120.pkl', "rb" ) )
print(len(parID_list))
entropy_dict = {}
for parID in parID_list:
    filename = '2Dfinal_circuit2_variant0_boundary1_ca_generalID%s_L8_J80_T120_N23880.pkl'%parID
    entropy_dict[parID]=entropy(parID,filename,results_path,show_fig=False)
entropy_dict = sorted(entropy_dict.items(), key=lambda x: x[1][2],reverse=True)
pickle.dump( entropy_dict, open( "entropy_ordered_dict_fulldataset.pkl", "wb" ) )
print(entropy_dict[:10])
