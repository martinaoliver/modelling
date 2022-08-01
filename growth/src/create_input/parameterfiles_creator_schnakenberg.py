##########################
#########README##########
##########################
# Generate parameter sets using latin hypercube sampling in a loguniform distribution.
# run in commandline ' python parameterfiles_creator.py '
# the number of samples is defined below in the 'numbercombinations' variable.
# $1 number of parameter combinations

#############
###paths#####
#############
import sys
import os

from importlib_metadata import distribution
pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')
#############

#############
###imports###
#############
import pickle
import pandas as pd
import numpy as np
from tqdm import tqdm



#######################
#########CODE##########
#######################

circuit_n='schnakenberg'
variant=0 #variant1 is variant0 but with varying kce as well.
np.random.seed(1)


def loguniform(low=-3, high=3, size=None):
    return (10) ** (np.random.uniform(low, high, size))


def uniform(low=-3, high=3, size=None):
    return np.random.uniform(low, high, size)


def lhs(data, nsample):
    m, nvar = data.shape
    ran = np.random.uniform(size=(nsample, nvar))
    s = np.zeros((nsample, nvar))
    for j in tqdm(range(0, nvar)):
        idx = np.random.permutation(nsample) + 1
        P = ((idx - ran[:, j]) / nsample) * 100
        s[:, j] = np.percentile(data[:, j], P)

    return s


def parameterfile_creator_function(numbercombinations):
    loguniformdist = loguniform(size=1000000)

    c_range = (0.001,100)
    d_B_range = (0.001,10)

    c_distribution = [x for x in loguniformdist if c_range[0] <= x <= c_range[1]]
    d_B_distribution = [x for x in loguniformdist if d_B_range[0] <= x <= d_B_range[1]]


    lenghtsdistributions = ( len(c_distribution), len(d_B_distribution))
    minimumlenghtdistribution = np.amin(lenghtsdistributions)
    c_distribution = c_distribution[:minimumlenghtdistribution]
    d_B_distribution = d_B_distribution[:minimumlenghtdistribution]

    # A general matrix is generated with the distributions for each parameter
    c_matrix = np.column_stack((c_distribution, c_distribution, c_distribution, c_distribution,))
    d_B_matrix = np.column_stack( (d_B_distribution)).transpose()  # needs transposing as only one element leads to np.array the other way around.

    # par_distribution_matrix = np.concatenate((Vm_matrix, km_matrix, mu_matrix), 1)
    par_distribution_matrix = np.concatenate(( c_matrix,d_B_matrix), 1)
    #
    points = lhs(par_distribution_matrix, numbercombinations)

    d_A = np.full((numbercombinations, 1), 1)
    parameterindex = np.arange(1, numbercombinations + 1, dtype=int).reshape(numbercombinations, 1)
    points = np.concatenate((parameterindex, points, d_A), 1)
    parameternames = (
    'index', 'c1', 'c2', 'c3', 'c4', 'd_B', 'd_A')
    print(np.shape(points), np.shape(parameternames))
    df = pd.DataFrame(data=points, columns=parameternames)
    df['index'] = df['index'].astype(int)    
    df = df.set_index('index')

    return df

# number_of_dataframes=64

n_param_sets = int(sys.argv[1]) # this number needs to be a multiple of 64
# n_param_sets = int(numbercombinations / number_of_dataframes)
# count = 0
# n_batches= 64
# batch_size = int(n_param_sets/n_batches)

# for n in tqdm(range(n_batches)):
for n in range(1):
    # df = parameterfile_creator_function(batch_size)
    df = parameterfile_creator_function(n_param_sets)
    pickle.dump(df, open(modellingpath + '/growth/input/parameterfiles/df_circuit%s_variant%r_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), 'wb'))
    print ('Parameterfile a %r created' %n)
    print(df)
