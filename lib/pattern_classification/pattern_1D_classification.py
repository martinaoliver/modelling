import numpy as np
from scipy.signal import find_peaks
from numerical.cn_plot import plot1D, surfpattern



def NormalizeData(data):
    return (data - np.min(data)) / (np.max(data) - np.min(data)+1e-8)

def countPeaks(U, showPlot1D=True):
    peaks = [0,0]
    peaks[0], _ = find_peaks(U[0], prominence=0.1)
    peaks[1], _ = find_peaks(U[1], prominence=0.1)
    if showPlot1D == True:
        plot1D(U,plotPeaks=True, peaks=peaks)

    return peaks

def patternClassification(U_final, U_record, normalize=True):
    #check if flat
    relRangeFlat = [(np.amax(U) - np.amin(U))/(np.amax(U)+1e-8) for U in U_final]
    if any(i<0.001 for i in relRangeFlat):
        flat=True
    else:
        flat=False

    #check if converged
    relRangeConverged=[0,0]
    for count,Ux_record in enumerate(U_record):
        relRangeConverged[count] = [(np.amax(x) - np.amin(x))/(np.amax(x)+1e-8) for x in np.transpose(Ux_record[-3:])]
    # if np.amax(relRangeConverged[0])>0.001 or np.amax(relRangeConverged[1])>0.001:
    if np.amax(relRangeConverged[0])>0.01 or np.amax(relRangeConverged[1])>0.01:
        converged=False
    else:
        converged=True

    #check if regular
    U_final_norm = [NormalizeData(U) for U in U_final]
    peaks = countPeaks(U_final_norm, showPlot1D=False)

    #calculate distance between peaks in peak0
    std=[0,0]
    for count,singleUpeak in enumerate(peaks):
        if len(singleUpeak)>2:
            peak_dist = [np.linalg.norm(singleUpeak[i]-singleUpeak[i+1]) for i in range(len(singleUpeak)-1)]
            peak_dist = peak_dist/np.sum(peak_dist)
            std[count] = np.std(peak_dist)
        else:
            std[count] = 1
    if std[0]<0.01 and std[1]<0.01:
        regular=True
    else:
        regular=False
        


    if flat==True and converged==True:
        pattern='Homogeneous'
    if flat==True and converged==False:
        pattern='Temporal Oscillator'
    if flat==False and converged==True and regular==True:
        pattern = 'Stationary regular pattern'
    if flat==False and converged==True and regular==False:
        pattern = 'Stationary irregular pattern'

    if flat==False and converged==False and regular==True:
        pattern = 'Non-Stationary regular pattern'
    if flat==False and converged==False and regular==False:
        pattern = 'Non-Stationary irregular pattern'
    # if var[0]>0.1 and var[1]>0.1:

    return pattern, converged, flat, regular