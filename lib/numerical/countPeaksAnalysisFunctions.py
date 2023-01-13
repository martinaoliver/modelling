#############
###paths#####
#############
import sys
import os

pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')
#############

from numerical.cn_plot import plot1D
from scipy.signal import find_peaks
import numpy as np

# def countPeaks(U, height = 1, threshold=0.1, prominence=0.1, plot=True):
def countPeaks(U, showPlot1D=False):
    # height = 0.5
    # threshold = 0.01
    # prominence = 0.01

    peaks = [0,0]
    # peaks[0], _= find_peaks(U[0], height = height,threshold = threshold,prominence = prominence*np.amax(U[0]))
    # peaks[1], _ = find_peaks(U[1], height= height, threshold = threshold,prominence = prominence*np.amax(U[1]))


    peaks[0], _ = find_peaks(U[0],prominence=0.1)
    peaks[1], _ = find_peaks(U[1],prominence=0.1)

    
    if showPlot1D==True:
        plot1D(U, round=False, peaks=peaks, pad=0.2)
    return peaks

def varPeakDistFunction(U, showPlot1D = False, printVar=False):
    peaks = countPeaks(U, showPlot1D=showPlot1D)#,height = height, threshold=threshold, prominence=prominence )
    #calculate distance between peaks
    peak0 = peaks[0]
    #calculate distance between peaks in peak0
    var=[0,0]
    for count,peak in enumerate(peaks):
        if len(peak)>2:
            peak_dist = [np.linalg.norm(peak[i]-peak[i+1]) for i in range(len(peak)-1)]
            peak_dist = peak_dist/np.sum(peak_dist)
            var[count] = np.var(peak_dist)
        else:
            var[count] = 1
    if printVar==True:
        print(var)

    return var