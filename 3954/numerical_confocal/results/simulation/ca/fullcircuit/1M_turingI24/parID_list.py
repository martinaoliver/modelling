import os
# files = [f for f in os.listdir('.') if os.path.isfile(f)]
import glob
import pickle
files = glob.glob('./2Dfinal_circuit2_variant9_bc1.5_ca_fullcircuitID*_L2_J100_T24_N2400*')
print(len(files))
parID_list = []
for f in files:

    f0 = f.rpartition('ID')[2]
    f1 = f0.rpartition('_L')[0]
    parID_list.append(f1)
print(parID_list[:10])
print(len(parID_list))

pickle.dump( parID_list, open( "parID_list_circuit2_variant9_bc1.5_ca_fullcircuit_L2_J100_T24_N2400.pkl", "wb" ) )
print('------')
print(parID_list.count('1'))
