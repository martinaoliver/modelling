import pickle
import numpy as np
import sys
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm
from scipy.signal import find_peaks, find_peaks_cwt
import os
pwd = os.getcwd()
root = pwd.split("home", 1)[0]
modelling_home = root + 'home/Documents/modelling'
modelling_ephemeral = root + 'ephemeral/Documents/modelling'
modulepath = modelling_home + '/6eq/modules'
sys.path.append(modulepath)


circuit_n=4
variant=0
parametersets_n = 10000

# open parameter dictionaries
general_df = pickle.load(open(modelling_home + '/3954/parameter_space_search/results/output_dataframes/lsa_df_circuit2_variant%r_%rparametersets.pkl'%(variant,parametersets_n), "rb"))

def plot_peaks(simulation,parID,print_stripped=True):

    green = simulation[-1]
    green = np.concatenate([np.zeros(1), green, np.zeros(1)])

    height = 1
    threshold = 0.1
    prominence = 0.1


    green_peaks, properties = find_peaks(green, height= height, threshold = threshold,prominence = prominence*np.amax(green))


    return len(green_peaks)

boundary_coef=int(sys.argv[1])
shape=str(sys.argv[2])
mechanism = 'lost_plasmid'
L,T,J,N = [8,100,64,12600]

stripped_parID = []
nonstripped_parID = []

for parID in range(1,10001):
    filename = 'final_circuit%r_variant%r_boundary%r_%s_%sID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundary_coef, shape,mechanism,parID,L,J,T,N)
    simulation = pickle.load( open( modelling_ephemeral + '/3954/numerical_confocal/results/simulation/1M_df_%s/1D%s.pkl'%(mechanism,filename), "rb" ) )
    peaks_n = plot_peaks(simulation,parID,print_stripped=True)
    if np.amax(peaks_n) > 2:
        stripped_parID.append(parID)
    else:
        nonstripped_parID.append(parID)
print(nonstripped_parID, stripped_parID)

pickle.dump( nonstripped_parID, open( "%s_boundary%r_nonstripped_parID_%r.pkl"%(shape,boundary_coef,parametersets_n), "wb" ) )
pickle.dump( stripped_parID, open( "%s_boundary%r_stripped_parID_%r.pkl"%(shape,boundary_coef,parametersets_n), "wb" ) )
