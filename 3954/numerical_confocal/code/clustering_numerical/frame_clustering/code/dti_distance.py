from scipy import stats
import pickle
import matplotlib.pyplot as plt
import numpy as np

from dtaidistance import dtw
from dtaidistance import dtw_visualisation as dtwvis
from dtaidistance import ed

# numerical_confocal_path =  '/rds/general/user/mo2016/home/Documents/modelling/6eq/numerical_confocal'
ephemeral_numerical =  '/rds/general/user/mo2016/ephemeral/Documents/modelling/6eq/numerical_confocal'

import sys
import time
start_time = time.time()


#########
#QUERY
#########
#Opening list with parID's
file = open(ephemeral_numerical + '/results/simulation/redgreen/8x20_T300/parID_list.txt')
parID_list =file.read().splitlines() # will create the list from the file which should contain only names without '\n'
file.close()
parID_list = [int(i) for i in parID_list] #turn string list into integer list

#parID_list = parID_list[:10]




#########
#TEMPLATE
#########
template_singleframe_list = []
times = [65,77,102,110,124,135,149]
time_index = {}
for count, time in enumerate(times):
    single_frame = pickle.load( open( ephemeral_numerical + '/toy/fake/redgreen_mushroom_t%r.pkl'%time, 'rb' ) )
    template_singleframe_list.append(single_frame)
    time_index[time]=count


# query_parID = 3
# time = 149
window=int(sys.argv[1])
# index = time_index[time]

plot = False
if plot == True:
    plt.imshow(template_singleframe_list[index].astype('uint8'), origin= 'lower')
    plt.show()
    # plt.imshow(query_timeseries_list[query_parID][time].astype('uint8'), origin= 'lower')
    # plt.show()
def create_distance_list(time,template_singleframe_list):


    index = time_index[time]
    center_position = int(len(template_singleframe_list[0])/2)

    red_template, green_template = (template_singleframe_list[index][center_position,:center_position,i]  for i in range(2))

    distance_list = []

    parID_index = {}
    for index,parID in enumerate(parID_list):
     
        filename = 'redgreen_ADI_circuit2boundary1_growing_colony_turingID%r_L8_J160_T300_N241800'%parID
        timeseries_unstacked = pickle.load( open( ephemeral_numerical + '/results/simulation/redgreen/8x20_T300/%s.pkl'%filename, 'rb' ) )
        timeseries= np.stack( timeseries_unstacked, axis=0 )
        parID_index[parID]=index
        
        
        
        red_query, green_query= (timeseries[time,center_position,:center_position,i] for i in range(2))
        #sometimes array contains only zeros which gives a problem when normalising.
        if np.count_nonzero(red_query)==0:
            red_query = np.full_like(red_query,0.0001)
            print('zero_red')
        if np.count_nonzero(green_query)==0:
            green_query = np.full_like(green_query,0.0001)
            print('zero green')
        #z-normalisation
        red_template_z,red_query_z,green_template_z,green_query_z = (stats.zscore(a) for a in [red_template,red_query,green_template,green_query])
        #print('red')
        #print(red_template_z,red_query_z)
        #red_path,green_path = (dtw.warping_path(red_template, red_query,window=window), dtw.warping_path(green_template, green_query,window=window))
        red_distance, green_distance = (dtw.distance(red_template_z, red_query_z,window=window), dtw.distance(green_template_z, green_query_z,window=window))
        
        distance_list.append(int(np.mean([red_distance,green_distance])))

    return distance_list,parID_index

distance_list_list = []
for time in times:
    print('Comparing time %r' %time)
    distance_list,parID_index = create_distance_list(time,template_singleframe_list)
    distance_list_list.append(distance_list)

distance_list = np.mean(distance_list_list,axis=0)
print(distance_list)
distance_ordered_parID_list = [x for _,x in sorted(zip(distance_list,parID_list))]
pickle.dump( distance_ordered_parID_list, open( "distance_ordered_parID_list_window%r_normalised.pkl"%window, "wb" ) )


distance_ordered_parID_list = pickle.load( open('distance_ordered_parID_list_window%r_normalised.pkl'%window, 'rb' ) )

for time in times:
    print('plotting time %r'%time)
    index = time_index[time]
    n_col = 46
    num = len(distance_ordered_parID_list)      #  number of elements for each cluster
    n_row = np.floor(num/n_col)+1    # number of rows in the figure of the cluster
    fig = plt.figure(figsize=(50,50))

    ax=plt.subplot(n_row,n_col, 1)
    ax.imshow(template_singleframe_list[index].astype('uint8'), origin= 'lower')
    ax.title.set_text('T%rh'%time)
    ax.axis('off')


    for n in range(0, num):
        parID =distance_ordered_parID_list[n]
        filename = 'redgreen_ADI_circuit2boundary1_growing_colony_turingID%r_L8_J160_T300_N241800'%parID
        timeseries_unstacked = pickle.load( open( ephemeral_numerical + '/results/simulation/redgreen/8x20_T300/%s.pkl'%filename, 'rb' ) )
        #filename = 'redgreen_ADI_circuit2boundary1_growing_colony_turingID%r_L8_J160_T300_N241800'%parID
        #timeseries_unstacked = pickle.load( open( ephemeral_numerical + '/results/simulation/redgreen/8x20_T300/%s.pkl'%filename, 'rb' ) )
        timeseries= np.stack( timeseries_unstacked, axis=0 )


        ax=plt.subplot(n_row,n_col, n+2)
        ax.imshow(timeseries[time].astype('uint8'), origin= 'lower')
        ax.set_title('%r'%parID, fontsize=5)

        ax.axis('off')

    fig.subplots_adjust(hspace=-0.6,wspace = -0.1)
    plt.tight_layout()
    plt.savefig('../figures/dtw_ranking/dtw_window%r_t%r_fulldataset_normalised.png'%(window,time))
    plt.close()
    plt.clf()


    #
#
#


#
# #green
# fig,(ax00,ax2,ax3,ax4) = plt.subplots(1, 4,figsize=(20,4))
# path = dtw.warping_path(green_template, green_query,window=10)
# fig1, b  = dtwvis.plot_warping(green_template, green_query, path)#, filename="warp_green.png")
# # fig.suptitle("Vicky")
# # ax00 = fig.add_subplot(gs00[i, j])
# # fig[0].title.set_text('%rh'%1)
# distance = dtw.distance(green_template, green_query,window=10)
# print(distance)
#
#
# #red
# path = dtw.warping_path(red_template, red_query,window=10)
# dtwvis.plot_warping(red_template, red_query, path)#, filename="warp_red.png")
# distance = dtw.distance(red_template, red_query,window=10)
# print(distance)
# plt.show()
