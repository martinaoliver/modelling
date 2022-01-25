
from tslearn.generators import random_walks
from sklearn.cluster import KMeans
import pickle
from tslearn.metrics import dtw
from tslearn.metrics import dtw_path
import matplotlib.pyplot as plt
import numpy as np
from tslearn.utils import to_time_series, to_time_series_dataset

numerical_confocal_path =  '/rds/general/user/mo2016/home/Documents/modelling/6eq/numerical_confocal'
ephemeral_numerical_results =  '/rds/general/user/mo2016/ephemeral/Documents/modelling/6eq/numerical_confocal/results/simulation/redgreen/8x20_T300'

import time
start_time = time.time()
#Opening list with parID's
file = open(ephemeral_numerical_results + '/parID_list.txt')
parID_list =file.read().splitlines() # will create the list from the file which should contain only names without '\n'
file.close()
parID_list = [int(i) for i in parID_list] #turn string list into integer list

timeseries_unstacked_list= []
timeseries_list = []

for parID in parID_list:
    # print(parID)
    filename = 'redgreen_ADI_circuit2boundary1_growing_colony_turingID%r_L8_J160_T300_N241800'%parID
    timeseries_unstacked = pickle.load( open( ephemeral_numerical_results + '/%s.pkl'%filename, 'rb' ) )
    timeseries_unstacked_list.append(timeseries_unstacked)

#     #stack timepoints
    timeseries= np.stack( timeseries_unstacked, axis=0 )
#
#     #obtain red and green and concatenate them
    timeseries_r, timeseries_g = (timeseries[160,80,:80,i] for i in range(2))
    timeseries_rg = np.concatenate((timeseries_g,timeseries_r),axis=-1)
    timeseries_list.append(timeseries_rg)
#
# # stack parID's
data = np.stack(timeseries_list, axis=0)

distortions = []
K = range(1,21)
for k in K:
    print(k)
    model = KMeans(n_clusters=k)
    fit = model.fit_predict(data)
    pickle.dump(model, open('../trained_model/cluster%r_KMeans_1.pkl'%k, "wb"))


# plt.figure(figsize=(16,8))
# plt.plot(K, distortions, 'bx-')
# plt.xlabel('k')
# plt.ylabel('Distortion')
# plt.title('The Elbow Method showing the optimal k')
# plt.show()
# plt.savefig('Elbowmethod k10-20')
# n=5
# model = KMeans(n_clusters=n)
# Z = model.fit_predict(data)


# print(model.labels_)
# print(fit)
# k=1
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



print("--- %s seconds ---" % (time.time() - start_time))
