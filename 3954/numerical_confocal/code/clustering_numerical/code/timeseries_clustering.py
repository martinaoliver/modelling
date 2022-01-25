# from tslearn.clustering import TimeSeriesKMeans
# model = TimeSeriesKMeans(n_clusters=3, metric="dtw", max_iter=10)
# model.fit(data)
from tslearn.generators import random_walks
from tslearn.clustering import TimeSeriesKMeans
import pickle
from tslearn.metrics import dtw
from tslearn.metrics import dtw_path
import matplotlib.pyplot as plt
import numpy as np
numerical_confocal_path =  '/rds/general/user/mo2016/home/Documents/modelling/6eq/numerical_confocal'
ephemeral_numerical_results =  '/rds/general/user/mo2016/ephemeral/Documents/modelling/6eq/numerical_confocal/results/simulation/redgreen/8x20_T300'

#Opening list with parID's
file = open(ephemeral_numerical_results + '/parID_list.txt')
parID_list =file.read().splitlines() # will create the list from the file which should contain only names without '\n'
file.close()
parID_list = [int(i) for i in parID_list] #turn string list into integer list


timeseries_unstacked_list= []
timeseries_list = []
for parID in parID_list:
    filename = 'redgreen_ADI_circuit2boundary1_growing_colony_turingID%r_L8_J160_T300_N241800'%parID
    timeseries_unstacked = pickle.load( open( ephemeral_numerical_results + '/%s.pkl'%filename, 'rb' ) )
    timeseries_unstacked_list.append(timeseries_unstacked)

    #stack timepoints
    timeseries= np.stack( timeseries_unstacked, axis=0 )


    #obtain red and green and concatenate them
    timeseries_r, timeseries_g = (timeseries[::10,80,:80,i] for i in range(2))
    # timeseries_r, timeseries_g = (timeseries[:10,80,:80,i] for i in range(2))

    timeseries_rg = np.concatenate((timeseries_g,timeseries_r),axis=-1)
    timeseries_list.append(timeseries_rg)

data = np.stack(timeseries_list, axis=0)

distortions = []
K = range(10,21)
for k in K:
    model = TimeSeriesKMeans(n_clusters=k, metric="softdtw",n_jobs=-1)#, max_iter=10)
    fit = model.fit_predict(data)
    distortions.append(model.inertia_)
    model.to_hdf5('../trained_model/cluster%r_trainedfull_..10h.hdf5'%k)
    print(k)

#plt.figure(figsize=(16,8))
#plt.plot(K, distortions, 'bx-')
#plt.xlabel('k')
#plt.ylabel('Distortion')
#plt.title('The Elbow Method showing the optimal k')
#plt.savefig('elbowmethod3.png')




# model = TimeSeriesKMeans(n_clusters=3, metric="dtw", max_iter=10)
# fit = model.fit_predict(data)
# inertia = model.inertia_


# def show_rgbvideo(n):
#     rgb_timeseries=timeseries_unstacked_list[n] # Read the numpy matrix with images in the rows
#     im=plt.imshow(rgb_timeseries[0].astype('uint8'), origin= 'lower')
#     for time in range(150):
#         im.set_data(rgb_timeseries[time].astype('uint8'))
#         plt.pause(0.01)
#     plt.show()

# show_rgbvideo(0)
# show_rgbvideo(3)
# show_rgbvideo(-1)
