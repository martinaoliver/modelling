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
circuit_n=2
variant=9 #variant1 is variant0 but with varying kce as well.

#######################
#########CODE##########
#######################


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

    # b_range = (0.1,1)
    Vm_range = (10, 1000)
    km_range = (0.1, 250)
    mu_range = (0.001, 50)
    # d_B_range = (0.001, 10)

    # b_distribution = [x for x in loguniformdist if b_range[0] <= x <= b_range[1]]
    Vm_distribution = [x for x in loguniformdist if Vm_range[0] <= x <= Vm_range[1]]
    km_distribution = [x for x in loguniformdist if km_range[0] <= x <= km_range[1]]
    mu_distribution = [x for x in loguniformdist if mu_range[0] <= x <= mu_range[1]]
    # d_B_distribution = [x for x in loguniformdist if d_B_range[0] <= x <= d_B_range[1]]


    lenghtsdistributions = ( len(Vm_distribution), len(km_distribution), len(mu_distribution))
    # lenghtsdistributions = ( len(Vm_distribution), len(km_distribution), len(mu_distribution), len(d_B_distribution))
    minimumlenghtdistribution = np.amin(lenghtsdistributions)
    # b_distribution = b_distribution[:minimumlenghtdistribution]
    Vm_distribution = Vm_distribution[:minimumlenghtdistribution]
    km_distribution = km_distribution[:minimumlenghtdistribution]
    mu_distribution = mu_distribution[:minimumlenghtdistribution]
    # d_B_distribution = d_B_distribution[:minimumlenghtdistribution]

    # A general matrix is generated with the distributions for each parameter

    # b_matrix = np.column_stack(
    #     (b_distribution, b_distribution, b_distribution, b_distribution, b_distribution, b_distribution))
    Vm_matrix = np.column_stack(
        (Vm_distribution, Vm_distribution, Vm_distribution, Vm_distribution, Vm_distribution, Vm_distribution))
    mu_matrix = np.column_stack((mu_distribution, mu_distribution))
    km_matrix = np.column_stack((km_distribution, km_distribution, km_distribution, km_distribution, km_distribution,
                                 km_distribution, km_distribution))
    # d_B_matrix = np.column_stack( (d_B_distribution)).transpose()  # needs transposing as only one element leads to np.array the other way around.

    par_distribution_matrix = np.concatenate((Vm_matrix, km_matrix, mu_matrix), 1)
    # par_distribution_matrix = np.concatenate(( Vm_matrix, km_matrix, mu_matrix,d_B_matrix), 1)
    #
    points = lhs(par_distribution_matrix, numbercombinations)

    bx = np.full((numbercombinations, 1), 0.01)
    cooperativity = np.full((numbercombinations, 1), 2)
    d_A = np.full((numbercombinations, 1), 2)
    d_B = np.full((numbercombinations, 1), 0.8)
    parameterindex = np.arange(1, numbercombinations + 1, dtype=np.int).reshape(numbercombinations, 1)
    points = np.concatenate((parameterindex, bx, bx, bx, bx, bx, bx, points, d_A, d_B, cooperativity), 1)

    parameternames = (
    'index', 'ba', 'bb', 'bc', 'bd', 'be', 'bf', 'Va', 'Vb', 'Vc', 'Vd', 'Ve', 'Vf', 'kaa', 'kda', 'keb', 'kbd', 'kce',
    'kfe', 'kee', 'mua', 'mulva', 'd_A', 'd_B', 'n')
    df = pd.DataFrame(data=points, columns=parameternames)
    df['index'] = df['index'].astype(int)
    df['n'] = df['n'].astype(int)
    df = df.set_index('index')

    return df

# number_of_dataframes=64

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
    # print(df)
    print(df)
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
