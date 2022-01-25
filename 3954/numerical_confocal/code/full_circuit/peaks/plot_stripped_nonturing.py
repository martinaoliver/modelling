import pickle
import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
ephemeral_numerical =  '/rds/general/user/mo2016/ephemeral/Documents/modelling/6eq/numerical_confocal'

striped_parID = pickle.load(open( "peakresults/striped_parID_nonturingfulldataset.pkl", "rb" ) )
nonstriped_parID = pickle.load(open( "peakresults/nonstriped_parID_nonturingfulldataset.pkl", "rb" ) )


n_col = 46
n_row = np.floor(len(striped_parID)/n_col)+1    # number of rows in the figure of the cluster
fig1 = plt.figure(figsize=(50,50))

for count,parID in enumerate(striped_parID):
    filename = 'redgreen_ADI_circuit2boundary1_growing_colony_nonturingID%r_df0_L8_J64_T150_N18900.pkl'%parID
    timeseries_unstacked = pickle.load( open( ephemeral_numerical + '/results/simulation/redgreen/nonturing/8x8_T150/%s'%filename, 'rb' ) )
    ax=plt.subplot(n_row,n_col, count)
    ax.imshow(timeseries_unstacked[149].astype('uint8'), origin= 'lower')
    ax.axis('off')
    ax.set_title(str(parID),fontsize=5)
fig.subplots_adjust(hspace=-0.6,wspace = -0.1)
plt.tight_layout()
plt.savefig('peakresults/stripped_nonturing_sample')
plt.close()
plt.clf()

n_col = 46
n_row = np.floor(len(nonstriped_parID)/n_col)+1    # number of rows in the figure of the cluster
fig2 = plt.figure(figsize=(50,50))

for count,parID in enumerate(nonstriped_parID):
    filename = 'redgreen_ADI_circuit2boundary1_growing_colony_nonturingID%r_df0_L8_J64_T150_N18900.pkl'%parID
    timeseries_unstacked = pickle.load( open( ephemeral_numerical + '/results/simulation/redgreen/nonturing/8x8_T150/%s'%filename, 'rb' ) )
    ax=plt.subplot(n_row,n_col, count)
    ax.imshow(timeseries_unstacked[149].astype('uint8'), origin= 'lower')
    ax.axis('off')
    ax.set_title(str(parID),fontsize=5)
fig.subplots_adjust(hspace=-0.6,wspace = -0.1)
plt.tight_layout()
plt.savefig('peakresults/nonstripped_nonturing_sample')
plt.close()

# for i in range(0,k):
#     n_col = 10
#     row = np.where(fit==i)[0]
#     print(row) # row in Z for elements of cluster i
#     num = row.shape[0]       #  number of elements for each cluster
#     n_row = np.floor(num/n_col)+1    # number of rows in the figure of the cluster
#     print("cluster "+str(i))
#     print(str(num)+" elements")
#     fig = plt.figure(figsize=(10,10))
#     for n in range(0, num):
#
#         ax=plt.subplot(n_row,n_col, n+1)
#         rgb_timeseries=timeseries_unstacked_list[row[n]] # Read the numpy matrix with images in the rows
#         ax.imshow(rgb_timeseries[160].astype('uint8'), origin= 'lower')
#         ax.axis('off')
#     plt.show()
