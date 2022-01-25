import sys
import os
pwd = os.getcwd()
root = pwd.split("home", 1)[0]
modelling_home = root + 'home/Documents/modelling'
modelling_ephemeral = root + 'ephemeral/Documents/modelling'
modulepath = modelling_home + '/3954/modules'
sys.path.append(modulepath)

from numerical_solvers_variableboundary import *
import pickle
import pandas as pd
import seaborn as sns
# general_df = pickle.load(open('bullseye_df.pkl', "rb"))
general_df = pickle.load(open('/Volumes/mo2016/home/Documents/modelling/3954/parameter_space_search/results/output_dataframes/interesting_variations/ATC_24240.pkl', "rb"))
# general_df = pickle.load(open('interesting_parIDs.pkl', "rb"))
# parID = int(sys.argv[1])
# print(general_df.loc[parID])
print(general_df['kce'])
