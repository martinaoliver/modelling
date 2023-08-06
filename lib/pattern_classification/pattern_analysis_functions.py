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
    peaks[0], _ = find_peaks(U[0], prominence=0.1)
    peaks[1], _ = find_peaks(U[1], prominence=0.1)

    # Calculate the wavelength
    wavelength_x = np.mean(np.diff(x_grid[peaks[0]]))
    wavelength_y = np.mean(np.diff(x_grid[peaks[1]]))
    avg_wavelength = np.mean([wavelength_x, wavelength_y])
    # Plot the 1D signal and peaks
    if showplot1D:
        plot1D(U, peaks=peaks)

    return  avg_wavelength


