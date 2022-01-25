import pickle
import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
ephemeral_numerical_results =  '/rds/general/user/mo2016/ephemeral/Documents/modelling/6eq/numerical_confocal/results'


file = open(ephemeral_numerical_results + '/simulation/redgreen/8x20_T300/parID_list.txt')
parID_list =file.read().splitlines() # will create the list from the file which should contain only names without '\n'
file.close()
parID_list = [int(i) for i in parID_list] #turn string list into integer list



def count_peaks(timeseries_unstacked):
    #fig,ax = plt.subplots(1, 3,figsize=(20,4))
    #ax[0].imshow(timeseries_unstacked[-1].astype('uint8'))

    red_diameter,green_diameter = (timeseries_unstacked[150][int(160/2),:,i] for i in range(2))

    #findpeaks parameters
    height = 50
    threshold = 0
    prominence = 20
    red_peaks, _ = find_peaks(red_diameter, height=height, threshold = threshold, prominence = prominence)
    green_peaks, _ = find_peaks(green_diameter, height=height, threshold = threshold, prominence = prominence)

    return (len(red_peaks),len(green_peaks))

parID_peak_dict = {}
striped_parID = []
nonstriped_parID = []
for parID in parID_list:
    print(parID)
    filename = 'redgreen_ADI_circuit2boundary1_growing_colony_turingID%r_L8_J160_T300_N241800.pkl'%parID
    timeseries_unstacked = pickle.load( open(ephemeral_numerical_results + '/simulation/redgreen/8x20_T300/%s'%filename, 'rb' ) )
    peaks_n = count_peaks(timeseries_unstacked)
    parID_peak_dict[parID] = peaks_n
    if any(n>2 for n in peaks_n):
        striped_parID.append(parID)
    else:
        nonstriped_parID.append(parID)


print(len(striped_parID),len(nonstriped_parID))
pickle.dump( striped_parID, open( "peakresults/striped_parID_turingfulldataset.pkl", "wb" ) )
pickle.dump( nonstriped_parID, open( "peakresults/nonstriped_parID_turingfulldataset.pkl", "wb" ) )
