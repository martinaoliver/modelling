import os
# files = [f for f in os.listdir('.') if os.path.isfile(f)]
import glob
import pickle
files = glob.glob('./2Dfinal*L5_J50_T2000_N20000*')
print(len(files))
parID_list = []
for f in files:

    f0 = f.rpartition('ID')[2]
    f1 = f0.rpartition('_L')[0]
    parID_list.append(f1)
print(parID_list[:10])
print(len(parID_list))

pickle.dump( parID_list, open( "parID_list_L5_J50_T2000_N20000.pkl", "wb" ) )
print('------')
print(parID_list.count('1'))
