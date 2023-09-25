#%%
import sys
import os

pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')

from database.databaseFunctions import *
import pickle

query = '''SELECT system_class, COUNT(*) AS count
FROM analytical_output ao
join model_param mp on ao.model_param_id = mp.model_param_id
where mp.variant='9'
    GROUP BY system_class'''


# query = '''SELECT system_class, COUNT(*) AS count
# FROM analytical_output 
#     GROUP BY system_class'''

output = dict(general_query(query))
#%%

threshold = int(sys.argv[1])

if any([value < threshold for value in output.values()]):
    print('true')
else:
    print('false')



# %%
