##########################
#########README##########
##########################
# Generate parameter sets using latin hypercube sampling in a loguniform distribution.
# run in commandline ' python parameterfiles_creator.py '. 64 dataframes will be generated with a specific number of samples.
# the number of samples is defined below in the 'numbercombinations' variable.
# $1 number of parameter combinations


##########################
#########IMPORTS##########
##########################

# import module folder containing general functions used frequently
import os.path
import sys

# CHANGES FROM CIRCUIT1 TO CIRCUIT2 IN PARAMETERS:
# - CONSTITUTIVE NODE A, BUT THIS DOESNT AFFECT PARAMETER
# - IN THIS CASE WE WILL TEST LOWER VALUES OF D_B AS POTENTIALLY DIFFUSION IS LOWER EXPERIMENTALLY
# - IN PC CIRCUIT DIFFUSION OF PC AND HSL IS TUNED SIMULTANEOUSLY AND THEREFORE MUA = MUB

# modulepath = os.path.expanduser(
#     '~/Documents/modelling/6eq/modules')  # os.path.expanduser(path) : return the argument with an initial component of ~ or ~user replaced by that user’s home directory.
# sys.path.append(modulepath)

# other imports
import pickle
import pandas as pd
import numpy as np
from datetime import date
from tqdm import tqdm
import time
start_time = time.time()
# date = str(date.today())
# path = os.path.expanduser(
#     '~/Documents/modelling/6eq/parameter_space_search')  # path of project folder: folder where code, results and parameterfiles are found.
circuit_n=11
variant=8 #variant1 is variant0 but with varying kce as well.

#######################
#########CODE##########
#######################


def loguniform(low=-3, high=3, size=None):
    return (10) ** (np.random.uniform(low, high, size))


def uniform(low=-3, high=3, size=None):
    return np.random.uniform(low, high, size)


def lhs(data, nsample):
    m, nvar = len(data),1
    ran = np.random.uniform(size=(nsample, nvar))
    s = np.zeros((nsample, nvar))
    for j in range(0, nvar):
        idx = np.random.permutation(nsample) + 1
        P = ((idx - ran[:, j]) / nsample) * 100
        s[:, j] = np.percentile(data[:, j], P)
    return s

def parameterfile_creator_function(numbercombinations):
    loguniformdist = loguniform(size=1000000)
    parameters_lhs = ['Vm','km','mu','R','d_H','a','kD']

    Vm_range = (10, 1000)
    km_range = (0.1, 250)
    mu_range = (0.25, 2)
    R_range = (10, 200)
    d_H_range = (0.1, 10)
    a_range = (0.1,100)
    kD_range = (28,1000)

    parameter_range_list = [Vm_range,km_range,mu_range,R_range,d_H_range,a_range,kD_range]

    parameter_distribution_list = []
    for parameter_range in parameter_range_list:
        distribution = [x for x in loguniformdist if parameter_range[0] <= x <= parameter_range[1]]
        parameter_distribution_list.append(distribution)

    minimumlenghtdistribution = np.amin([len(x) for x in parameter_distribution_list])
    for count,parameter_distribution in enumerate(parameter_distribution_list):
        globals()[f'{parameters_lhs[count]}_distribution'] =  np.column_stack( (parameter_distribution[:minimumlenghtdistribution])).transpose()

    VmB,VmD,VmE,VmF = [lhs(Vm_distribution, numbercombinations) for i in range(4)]
    kemb,khmd,kfme,keme = [lhs(km_distribution, numbercombinations) for i in range(4)]
    aB,aD,aE,aF = [lhs(a_distribution, numbercombinations) for i in range(4)]
    muH = lhs(mu_distribution,numbercombinations)
    d_H = lhs(d_H_distribution,numbercombinations)
    R = lhs(R_distribution,numbercombinations)
    kD = lhs(kD_distribution,numbercombinations)

    # VmF = lhs(Vm_distribution,numbercombinations)
    # VmF = lhs(Vm_distribution)
    bmX = np.full((numbercombinations, 1), 0.1)#nM
    d_A = np.full((numbercombinations, 1), 1)#mm/h
    muLVA = np.full((numbercombinations, 1), 1.08)#h^-1
    muRNA= np.full((numbercombinations, 1), 12)#h^-1
    nH = np.full((numbercombinations, 1), 1)
    nE = np.full((numbercombinations, 1), 3)
    nF = np.full((numbercombinations, 1), 3)
    alphaH = np.full((numbercombinations, 1), 3000)
    # d_B = np.full((numbercombinations, 1), 1.8)
    index = np.arange(1, numbercombinations + 1, dtype=np.int).reshape(numbercombinations, 1)
    # parameternames = ['index','bmB','bmD','bmE','bmF','VmB','VmD','VmE','VmF','kemb','khmd','kfme','keme','kD','aB','aD','aE','aF','alphaH','muRNA','muLVA','muH','R','d_H','n']
    parameternames = ['index','bmB','bmD','bmE','bmF','VmB','VmD','VmE','VmF','kemb','khmd','kfme','keme','kD','aB','aD','aE','aF','alphaH','muRNA','muLVA','muH','R','d_H','nH','nE','nF']
    points = np.concatenate((index,bmX,bmX,bmX,bmX,VmB,VmD,VmE,VmF,kemb,khmd,kfme,keme,kD,aB,aD,aE,aF,alphaH,muRNA,muLVA,muH,R,d_H,nH,nE,nF), 1)

    df = pd.DataFrame(data=points, columns=parameternames)
    df['index'] = df['index'].astype(int)
    # df['n'] = df['n'].astype(int)
    df = df.set_index('index')

    return df

# number_of_dataframes=64

# n_param_sets = 10  # this number needs to be a multiple of 64
n_param_sets = int(sys.argv[1])  # this number needs to be a multiple of 64
# n_param_sets = int(numbercombinations / number_of_dataframes)
# count = 0
# n_batches= 64
# batch_size = int(n_param_sets/n_batches)

# for n in tqdm(range(n_batches)):
for n in range(1):
    # df = parameterfile_creator_function(batch_size)
    df = parameterfile_creator_function(n_param_sets)
    # df['kce'] = 10 ** 4
    # pickle.dump(df, open('df_circuit%r_variant%r_%rparametersets_batch%r.pkl'%(circuit_n,variant,n_param_sets,n), 'wb'))
    pickle.dump(df, open('df_circuit%r_variant%r_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), 'wb'))
    print ('Parameterfile a %r created' %n)
    print(df.head)
#
# my_data = {}
#
# # Load all batch dataframes
# for batch in range(n_batches):
#     my_data[batch] = pickle.load(open('df_circuit%r_variant%r_%rparametersets_batch%r.pkl'%(circuit_n,variant,n_param_sets,batch), "rb" ) )
#
# # Join all batch results to large results dataframe
# concat_df = pd.concat(my_data.values(), ignore_index=True)
# pickle.dump(concat_df, open('df_circuit%r_variant%r_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), 'wb'))
print("--- %s seconds ---" % (time.time() - start_time))
