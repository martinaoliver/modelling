#############
###paths#####
#############
import sys
import os

pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
modellingephemeral = '/rds/general/ephemeral/user/mo2016/ephemeral/Documents/modelling'
sys.path.append(modellingpath + '/lib')
#############



import glob
import pickle


# Specify name of circuit and variant investigated
# circuit_n=14;variant='fitted1';n_species=6
circuit_n=14;variant='2nd';n_species=6
# Specifiy number of parameter sets in parameterset file to be loaded
n_samples = 1000000
Kce=100
# n_param_sets = 2000000
# balance = 'balanced'
# folder = 'circuit14variantfitted1'
folder = 'circuit14variant2ndBalancedKce100'
modelArgs = [circuit_n,variant,n_species,folder]

# Specifiy number of parameter sets in parameterset file to be loaded

# specify dimensions of system




# slow
L=20; dx =0.1; J = int(L/dx)
T =100; dt = 0.02; N = int(T/dt)
boundarycoeff = 1
divisionTimeHours=0.5
p_division=0.38;seed=1


# medium
L=20; dx =0.1; J = int(L/dx)
T =50; dt = 0.02; N = int(T/dt)
boundarycoeff = 1
divisionTimeHours=0.5
p_division=1;seed=1



# fast
L=20; dx =0.1; J = int(L/dx)
T =25; dt = 0.02; N = int(T/dt)
boundarycoeff = 1
divisionTimeHours=0.2
p_division=0.7;seed=1



shape = 'ca'
# filename= lambda parID: 'circuit%r_variant%s_bc%s_%s_ID%s_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundarycoeff, shape,parID,L,J,T,N)
# filename= lambda parID: 'circuit%r_variant%snsr%s_bc%s_%s_ID%s_L%r_J%r_T%r_N%r'%(circuit_n,variant,nsr,boundarycoeff, shape,parID,L,J,T,N)
# filename= lambda parID: 'circuit%r_variant%s_bc%s_%s_ID%s_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundarycoeff, shape,parID,L,J,T,N)
filename= lambda parID: 'circuit%r_variant%s_%sparametersets_balanced_Kce%s_bc%s_%s_ID%s_L%r_J%r_T%r_N%r'%(circuit_n,variant,n_samples,Kce,boundarycoeff, shape,parID,L,J,T,N)

# print(f'filename: {filename(1)}')
datafile= modellingpath + '/3954/paper/out/numerical/colonies/simulation/%s/2Dfinal_%s.pkl'%(folder,filename('*'))
print(datafile)
files = glob.glob(datafile)
print(len(files))
parID_list = []
for f in files:
    f0 = f.rpartition('ID')[2]
    f1 = f0.rpartition('_L')[0]
    parID_list.append(f1)
print(parID_list[:10])
print(len(parID_list))

pickle.dump( parID_list, open( modellingpath + '/3954/paper/out/numerical/colonies/simulation/%s/parID_list_%s.pkl'%(folder,filename('x')), "wb" ) )
print('------')
print(parID_list.count('1'))
