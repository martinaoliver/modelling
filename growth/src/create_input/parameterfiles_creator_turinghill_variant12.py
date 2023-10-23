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
from tkinter import N
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

circuit_n='turinghill'
variant=12
n_samples =int(sys.argv[1])
seed = int(sys.argv[2])
np.random.seed(seed)


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

    # b_range = (0.1,100)
    # Vm_range = (0.1, 100)
    # km_range = (0.1, 100)
    # mu_range = (0.01, 1)


    b_range = (10,10000)
    Vm_range = (10,10000)
    km_range = (0.1,100)
    mu_range = (1,100)

    # d_B_range = (10**(-3), 10**3)
    n_range=(2,4)

    b_distribution = [x for x in loguniformdist if b_range[0] <= x <= b_range[1]]
    Vm_distribution = [x for x in loguniformdist if Vm_range[0] <= x <= Vm_range[1]]
    km_distribution = [x for x in loguniformdist if km_range[0] <= x <= km_range[1]]
    mu_distribution = [x for x in loguniformdist if mu_range[0] <= x <= mu_range[1]]
    # d_B_distribution = [x for x in loguniformdist if d_B_range[0] <= x <= d_B_range[1]]
    n_distribution = [x for x in loguniformdist if n_range[0] <= x <= n_range[1]]

    # lenghtsdistributions = ( len(Vm_distribution), len(km_distribution), len(mu_distribution))
    # lenghtsdistributions = ( len(b_distribution), len(Vm_distribution), len(km_distribution), len(mu_distribution), len(d_B_distribution))
    lenghtsdistributions = ( len(b_distribution), len(Vm_distribution), len(km_distribution), len(mu_distribution), len(n_distribution))
    minimumlenghtdistribution = np.amin(lenghtsdistributions)
    b_distribution = b_distribution[:minimumlenghtdistribution]
    Vm_distribution = Vm_distribution[:minimumlenghtdistribution]
    km_distribution = km_distribution[:minimumlenghtdistribution]
    mu_distribution = mu_distribution[:minimumlenghtdistribution]
    # d_B_distribution = d_B_distribution[:minimumlenghtdistribution] 
    n_distribution = n_distribution[:minimumlenghtdistribution]

    # A general matrix is generated with the distributions for each parameter
    b_matrix = np.column_stack(( b_distribution, b_distribution))
    Vm_matrix = np.column_stack((Vm_distribution, Vm_distribution))
    km_matrix = np.column_stack((km_distribution, km_distribution, km_distribution, km_distribution,))
    mu_matrix = np.column_stack((mu_distribution, mu_distribution))
    # d_B_matrix = np.column_stack( (d_B_distribution)).transpose()  # needs transposing as only one element leads to np.array the other way around.
    n_matrix = np.column_stack( (n_distribution)).transpose() 

    # par_distribution_matrix = np.concatenate((Vm_matrix, km_matrix, mu_matrix), 1)
    par_distribution_matrix = np.concatenate((b_matrix, Vm_matrix, km_matrix, mu_matrix,n_matrix), 1)
    #
    points = lhs(par_distribution_matrix, numbercombinations)

    # bx = np.full((numbercombinations, 1), 0.01)
    # cooperativity = np.full((numbercombinations, 1), 3)
    d_A = np.full((numbercombinations, 1), 10)
    d_B= np.full((numbercombinations, 1), 0.01)
    # d_B = np.full((numbercombinations, 1), 0.8)
    parameterindex = np.arange(0, numbercombinations , dtype=int).reshape(numbercombinations, 1)
    points = np.concatenate((parameterindex, points, d_A, d_B), 1)
    parameternames = (
    'index', 'ba', 'bb', 'Va', 'Vb', 'kaa', 'kba', 'kab', 'kbb', 'mua', 'mub','n','d_A', 'd_B')
    print(np.shape(points), np.shape(parameternames))
    df = pd.DataFrame(data=points, columns=parameternames)
    df['index'] = df['index'].astype(int)
    # df['n'] = df['n'].astype(int)
    df = df.set_index('index')

    return df

for n in range(1):
    # df = parameterfile_creator_function(batch_size)
    df = parameterfile_creator_function(n_samples)
    df.index += seed*n_samples

    pickle.dump(df, open(modellingpath + '/growth/input/parameterfiles/df_%s_variant%r_%rparametersets_seed%s.pkl'%(circuit_n,variant,n_samples,seed), 'wb'))
    print ('Parameterfile a %r created' %n)
    print(df)



# add_to_db=True
# if add_to_db==True:
#     from database.databaseFunctions import modelParam_df_to_sql
#     modelParam_df_to_sql(df, circuit_n, str(variant), n_samples)
