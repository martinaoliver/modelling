import sys
import os
pwd = os.getcwd()
root = pwd.split("home", 1)[0]
modelling_home = root + 'home/Documents/modelling'
modelling_ephemeral = root + 'ephemeral/Documents/modelling'
modulepath = modelling_home + '/3954/modules'
sys.path.append(modulepath)



#importing functions from that module folder
from parametercombination_analysis import *
from randomfunctions import wavelenght_from_dispersion, plot_highest_dispersion
from class_circuit_eq import *
from lhs import *
from findsteadystates_functions import *
from dispersionrelation_functions import *
# from lsafunctions import jacobianlsa


#other imports
import time
import datetime
import scipy.io
import numpy as np
from numpy import linalg as LA
import pickle
import matplotlib.pyplot as plt
import re
import pandas as pd
#path of folder containing code, results, parameter files
path = modelling_home + '/3954/parameter_space_search'

from datetime import date, timedelta


circuit=2
variant=0
n_parametersets = 1000000
df = pd.read_pickle(path + "/results/output_dataframes/lsa_df_circuit%r_variant%r_%rparametersets.pkl"%(circuit,variant,n_parametersets))
print(df[0:2])
df_singleindex = df.xs(0, level=1)
print(len(df_singleindex))
print(len(df))


print(df['system_class'].value_counts())
print(df_singleindex['system_class'].value_counts())
print(df_singleindex['ss_n'].value_counts())
