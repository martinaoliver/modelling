#%%
import sys
import os

pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')

from database.databaseFunctions import *
import pickle

#############

# nsr_list=[0.01,0.05,0.1,0.2]
# for nsr in nsr_list:




    #%%
    # Specify name of circuit and variant investigated
    # circuit_n='14'
    # parID=4187715
    # # nsr=0.01
    # old_variant='fitted7'
    # variant=f'{old_variant}_gaussian{parID}_nsr{nsr}' #variant='fitted7'
    # Specifiy number of parameter sets in parameterset file to be loaded
    # n_samples = 2000 #n_samples = 13700000


circuit_n='\'14\'';variant='\'2nd\'';n_species=6;n_samples = 1000000;balance='\'Balanced\''



model_param_dict = {'circuit_n': circuit_n, 'variant': variant, 'n_samples': n_samples, 'balance': balance}



    # %%






# #%%
# import psycopg2 
# import pandas as pd
# credentials=f"postgresql://moliver:moliver@ld-rendres07.bc.ic.ac.uk/moliver"

# def query_analyticalOutput_df_from_sql(model_param_dict):
#     with psycopg2.connect(credentials) as conn:
#         with conn.cursor() as cursor:
#         # Build the SQL query dynamically based on the provided dictionaries
#             query = """
#             SELECT *
#             FROM analytical_output ao
#             JOIN model_param mp on mp.model_param_id = ao.model_param_id
#             WHERE 1=1
#             """
            
#             # Add filters for model_params
#             for key, value in model_param_dict.items():
#                 query += f"AND mp.{key} = {value}\n"
            

#             print('query',query)

#             # # Execute the query
#             # cursor.execute(query)
            
#             # # Fetch all the rows
#             # rows = cursor.fetchmany(2)
            
#             # # Print or process the rows as per your requirement
#             # for row in rows:
#             #     print(row)

#             df = pd.read_sql_query(query, conn)

#             # Close the cursor and the connection



#             return df




result_df = query_analyticalOutput_df_from_sql(model_param_dict)
print(result_df)
result_df
# %%
