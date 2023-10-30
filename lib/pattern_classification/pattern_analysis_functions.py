#############
###paths#####
#############
import sys
import os

pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')
#############

from numerical.cn_plot import plot1D, surfpattern
from scipy.signal import find_peaks
import numpy as np




def find_wavelenght(U,x_grid,showplot1D=True):
    peaks = [0, 0]
    prominence=0.05

    peaks[0], _ = find_peaks(U[0], prominence=prominence)
    peaks[1], _ = find_peaks(U[1], prominence=prominence)

    # Calculate the wavelength
    wavelength_x = np.mean(np.diff(x_grid[peaks[0]]))
    wavelength_y = np.mean(np.diff(x_grid[peaks[1]]))
    list_of_wavelength = np.array([wavelength_x, wavelength_y])
    avg_wavelength = np.mean([wavelength_x, wavelength_y])
    
    # Plot the 1D signal and peaks
    if showplot1D:
        plot1D(U, peaks=peaks)
    # if  one wavelength is not found, the other one will be picked up instead of the average
    if math.isnan(avg_wavelength):
        list_of_wavelength = list_of_wavelength[~np.isnan(list_of_wavelength)]
        if len(list_of_wavelength)>0:
            return list_of_wavelength[0]
        else: 
            return np.nan
    else:
        return  avg_wavelength


def find_convergence(U_record):
    #check if converged
    relRangeConverged=[0,0]
    for time in np.arange(199,0,-1):
        for count,Ux_record in enumerate(U_record):
            relRangeConverged[count] = [(np.amax(x) - np.amin(x))/(np.amax(x)+1e-8) for x in np.transpose(Ux_record[time:time+3])]
        # if np.amax(relRangeConverged[0])>0.001 or np.amax(relRangeConverged[1])>0.001:
        if np.amax(relRangeConverged[0])>0.05 or np.amax(relRangeConverged[1])>0.05:
            converged=False
            return time*10
        else:
            converged=True
