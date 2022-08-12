#############
###paths#####
#############
import sys
import os

from importlib_metadata import distribution
pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')
#############

from numerical.cn_plot import plot1D, surfpattern
from numerical.fourierAnalysisFunctions import powerspectrumFunction, entropyFunction, fourierAnalysisFunction
import pickle
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

circuit_n='turinghill'
variant= 0
n_species=2
mechanism='edgegrowth1'
L=50; x_gridpoints=5; J=L*x_gridpoints;I=J 
T=2000; t_gridpoints = 30; N=T*t_gridpoints #Number of timepoints
filename= lambda parID: '%s_variant%s_%s_ID%s_L%r_J%r_T%r_N%r'%(circuit_n,variant,mechanism,parID,L,J,T,N)

parID_list = pickle.load(open( modellingpath + '/growth/out/numerical/%s/%s/data/parID_list_%s.pkl'%(circuit_n,mechanism,filename('x')), "rb" ) )
start=0
stop=len(parID_list)
parID_list = [int(i) for i in parID_list[start:stop]] #turn string list into integer list
parID_list.sort() #sort from lower to higher values
parIDHpsDict = {}
test=False
for count,parID in enumerate(tqdm(parID_list, disable=True)):
    U = pickle.load( open(modellingpath + '/growth/out/numerical/%s/%s/data/2Dfinal_%s.pkl'%(circuit_n,mechanism,filename(parID)), 'rb'))
    # if test==True:
    #     plot1D(U)
    fourierAnalysisOut = [fourierAnalysisFunction(Ux) for Ux in U]
    meanH = np.mean([H for fft_U,ps,H in fourierAnalysisOut])
    parIDHpsDict[parID] = meanH
    if count%10000==0:
        pickle.dump(parIDHpsDict, open( modellingpath + '/growth/out/patternAnalysis/%s/%s/parIDHpsDict%s_batch%r.pkl'%(circuit_n,mechanism,filename('x'), count), 'wb'))

if test==False:
    pickle.dump(parIDHpsDict, open( modellingpath + '/growth/out/patternAnalysis/%s/%s/parIDHpsDict%s.pkl'%(circuit_n,mechanism,filename('x')), 'wb'))
    print('saved')

