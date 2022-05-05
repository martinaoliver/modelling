
#IMPORTS#
#############################
import sys
import os

from joblib import parallel_backend
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
L=5; x_gridpoints =10; J = L*x_gridpoints; dx = float(L)/float(J-1)
grid = np.array([j*dx for j in range(J)])
parID=0
def plot(parID,filename,results_path,grid):
    final_concentration = pickle.load( open( results_path + '/' + filename, "rb" ) )
    plt.pcolormesh(grid, grid, final_concentration[2], shading='auto')

def entropy(parID,filename,results_path, show_fig=True):
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
