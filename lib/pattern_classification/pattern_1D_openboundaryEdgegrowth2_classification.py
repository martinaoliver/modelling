import numpy as np
from scipy.signal import find_peaks
from numerical.cn_plot import plot1D, surfpattern


def NormalizeData(data):
    return (data - np.min(data)) / (np.max(data) - np.min(data)+1e-8)


def countPeaks(U, showPlot1D=True):
    peaks = [0,0]
    prominence=0.05
    peaks[0], _ = find_peaks(U[0], prominence=prominence)
    peaks[1], _ = find_peaks(U[1], prominence=prominence)
    if showPlot1D == True:
        plot1D(U,plotPeaks=True, peaks=peaks)

    return peaks

def patternClassification_openboundaryEdgegrowth2(U_final, showPlot1D=False):
    #check if regular
    # U_final_cropped = [U[50:-50] for U in U_final]
    U_final_cropped_norm = [NormalizeData(U) for U in U_final]
    peaks = countPeaks(U_final_cropped_norm, showPlot1D=showPlot1D)

    max_n_peaks = np.amax([len(morphogen_peaks) for morphogen_peaks in peaks])

    if max_n_peaks <= 1:
        pattern = 'no pattern, homogeneous'
    if max_n_peaks == 2:
        pattern = 'no pattern, boundary effect'
    if max_n_peaks == 3:
        pattern = 'weak pattern'
    if max_n_peaks == 4:
        pattern = 'intermediate pattern'
    if max_n_peaks >=5:
        pattern = 'strong pattern'
    
    return pattern, max_n_peaks