from tslearn.clustering import TimeSeriesKMeans
import pickle
import numpy as np
import matplotlib.pyplot as plt
ephemeral_numerical =  '/rds/general/user/mo2016/ephemeral/Documents/modelling/6eq/numerical_confocal'

#Opening list with parID's
file = open(ephemeral_numerical + '/results/simulation/redgreen/8x20_T300/parID_list.txt')
parID_list =file.read().splitlines() # will create the list from the file which should contain only names without '\n'
file.close()
parID_list = [int(i) for i in parID_list] #turn string list into integer list

#timeseries_unstacked_list= []


#for parID in parID_list:
    # print(parID)
    #filename = 'redgreen_ADI_circuit2boundary1_growing_colony_turingID%r_L8_J160_T300_N241800'%parID
    #timeseries_unstacked = pickle.load( open( ephemeral_numerical + '/results/simulation/redgreen/8x20_T300/%s.pkl'%filename, 'rb' ) )
    #timeseries_unstacked_list.append(timeseries_unstacked)


k=20
model = pickle.load(open('../trained_model/cluster%r_KMeans_1.pkl'%k, 'rb'))
fit = model.labels_


for i in range(0,k):
    n_col = 30
    row = np.where(fit==i)[0]
    print(row) # row in Z for elements of cluster i
    num = row.shape[0]       #  number of elements for each cluster
    n_row = np.floor(num/n_col)+1    # number of rows in the figure of the cluster
    print("cluster "+str(i))
    print(str(num)+" elements")
    fig = plt.figure(figsize=(30,30))
    for n in range(0, num):

        ax=plt.subplot(n_row,n_col, n+1)
        #rgb_timeseries=timeseries_unstacked_list[row[n]] # Read the numpy matrix with images in the rows
        par_ID = parID_list[row[n]]
        filename = 'redgreen_ADI_circuit2boundary1_growing_colony_turingID%r_L8_J160_T300_N241800'%par_ID
        timeseries_unstacked = pickle.load( open( ephemeral_numerical + '/results/simulation/redgreen/8x20_T300/%s.pkl'%filename, 'rb' ) )
        rgb_timeseries=timeseries_unstacked # Read the numpy matrix with images in the rows
        ax.imshow(rgb_timeseries[160].astype('uint8'), origin= 'lower')
        ax.axis('off')
    plt.title('Cluster %r/%r'%(i,k))
    plt.tight_layout()
    plt.savefig('../figures/kmeans_singleframe_k%rcluster%r.png'%(k,i))
    plt.clf()
    plt.close(fig)
