#############################
#IMPORTS#
#############################
import sys
import os
pwd = os.getcwd()
root = pwd.rpartition("mo2016")[0] + pwd.rpartition("mo2016")[1] #/Volumes/mo2016/ or '/Users/mo2016/' or '/rds/general/mo2016/'
print(root)
if root == '/Users/mo2016':
    modelling_ephemeral = '/Volumes/mo2016/ephemeral/Documents/modelling'
    modelling_home = '/Volumes/mo2016/home/Documents/modelling'
    modelling_local = root + '/Documents/modelling'
    import matplotlib as mpl
    mpl.use('tkagg')

if root == '/Volumes/mo2016' or root=='/rds/general/user/mo2016': #'/rds/general' or root=='/Volumes':
        modelling_ephemeral = root + '/ephemeral/Documents/modelling'
        modelling_home = root  + '/home/Documents/modelling'
        modelling_local = modelling_home

modulepath = modelling_local + '/3954/modules/new_CN'
print(root)
sys.path.append(modulepath)


from plotting_numerical import plot_redgreen_contrast
from tqdm import tqdm
import pickle
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib.animation as animation
import matplotlib as mpl

#############################
#Opening list with parID's

lhs=True

if lhs==True:
    folder = 'fullcircuit_5716gaussian/1M'
    variant=0
    details = '1M_turingI'
else:
    var=float(sys.argv[1])
    # var=0.23
    folder = 'fullcircuit_5716gaussian/var%s'%var
    variant='5716gaussian'
    details = 'var%r'%var


circuit_n=2
shape='ca'
mechanism = 'fullcircuit'
dimension='2D'

L=10; x_gridpoints =15; J = L*x_gridpoints
T =120; t_gridpoints = 10; N = T*t_gridpoints


data_path = modelling_home + '/3954/numerical_confocal/results/simulation/ca/%s'%(folder)

if lhs==True:    
    if folder == 'fullcircuit_5716gaussian/1M_turingI':
        print('yes')
        general_filename = 'circuit%r_variant%s_%s_%s_L%r_J%r_T%r_N%r'%(circuit_n,'1MturingI', shape,mechanism,L,J,T,N)
    else: 
        general_filename = 'circuit%r_variant%s_%s_%s_L%r_J%r_T%r_N%r'%(circuit_n,variant, shape,mechanism,L,J,T,N)
else:
    general_filename = 'circuit%r_variant%svar%s_%s_%s_L%r_J%r_T%r_N%r'%(circuit_n,variant,var, shape,mechanism,L,J,T,N)



# metric=str(sys.argv[1])
metric='ps_min'
parID_dict = pickle.load( open(modelling_home + "/3954/numerical_confocal/results/entropy/EntropyDicts/%s_dict_%s.pkl"%(metric,general_filename), "rb" ) )
print(modelling_home + "/3954/numerical_confocal/results/entropy/EntropyDicts/%s_dict_%s.pkl"%(metric,general_filename))
len(parID_dict)
parID_list = []
entropy_list = []
parID_dict = dict(sorted(parID_dict.items(), key=lambda item: item[1])) #important to sort out dictionary by HKS values
for key in parID_dict:
    parID_list.append(key)
    entropy_list.append(parID_dict[key])
# print(parID_list, entropy_list)
parID_list = [int(i) for i in parID_list] #turn string list into integer list

print(parID_list)
# k=20
num=len(parID_list)
n_col = int(np.sqrt(num))
# for i in range(0,k):
# n_col = 10
#     row = np.where(fit==i)[0]
#     print(row) # row in Z for elements of cluster i
#     num = row.shape[0]       #  number of elements for each cluster
n_row = np.floor(num/n_col)+1    # number of rows in the figure of the cluster
#     print("cluster "+str(i))
#     print(str(num)+" elements")
fig = plt.figure(figsize=(n_col/10+2, n_row/10+2))
dx = float(L)/float(J-1)
grid = np.array([j*dx for j in range(J)])


for count,parID in tqdm(enumerate(parID_list),disable=False):
    ax=plt.subplot(n_row,n_col, count+1)
    #rgb_timeseries=timeseries_unstacked_list[row[n]] # Read the numpy matrix with images in the rows
    # par_ID = parID_list[row[n]]
    if lhs==True:
        filename = 'circuit%r_variant%s_%s_%sID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant, shape,mechanism,int(parID),L,J,T,N)
    else:
        filename = 'circuit%r_variant%svar%r_%s_%sID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant,var, shape,mechanism,int(parID),L,J,T,N)

    final_concentration = pickle.load( open(data_path + '/2Dfinal_%s.pkl'%filename, 'rb' ) )
    rgb = plot_redgreen_contrast(final_concentration,L,mechanism,shape,filename,modelling_ephemeral,parID=parID,scale_factor=x_gridpoints,save_figure='LargeImage')
    # rgb_timeseries=timeseries_unstacked # Read the numpy matrix with images in the rows
    # ax.set_title(parID,size=0.1)
    ax.set(yticks=[],xticks=[],yticklabels=[],xticklabels=[])
    ax.imshow(rgb.astype('uint8'), origin= 'lower')
    ax.set_ylabel('%s-%s'%(parID,entropy_list[count]),size= 1,c='y', labelpad=0.35)




# plt.title('1M numerical search 0-%r'%num)
plt.savefig(modelling_home + '/3954/numerical_confocal/results/entropy/LargeImages/%s_%s.png'%(metric,general_filename), dpi=2000)
# plt.show()
# plt.savefig('h.png',dpi=2000)
print('gh')
# plt.clf()
# plt.close(fig)

