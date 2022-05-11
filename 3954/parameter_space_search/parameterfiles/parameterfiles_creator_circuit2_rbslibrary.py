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



# other imports
import pickle
import pandas as pd
import numpy as np
from datetime import date
from tqdm import tqdm
import time
from itertools import product

df = pickle.load( open( "df_circuit2_variant1_100parametersets.pkl", "rb" ) )
df['original_index'] =df.index
library_df = pd.DataFrame(columns=df.columns)
rbs_strenght = np.logspace(1, 3, num=8)
rbs_combinations = list(product(rbs_strenght, repeat=3))

for index, row in df.iterrows():
    for combination in rbs_combinations:
        # library_df = library_df.append(row)
        row['Va','Vb','Ve'] = combination
        library_df  = library_df.append(row)
library_df = library_df.reset_index(drop=True)

par_dict = library_df.loc[1].to_dict()
print(par_dict)
pickle.dump( library_df, open( "df_circuit2_variant1_100parametersets_rbslibrary0.pkl", "wb" ) )
print(len(library_df))
print(library_df)