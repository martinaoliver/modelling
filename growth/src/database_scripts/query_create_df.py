#%%
import sys
import os

pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')

from database.databaseFunctions import *
import pickle
import pandas as pd

#############

#%%


# query = '''select ao.ss_n, mp."parID", ao."ssID", mp.variant,simulation_param_uuid,  ao.model_param_id,ao.system_class,pco.pattern_class_nogrowth from pattern_class_output pco
#     join analytical_output ao on (pco.model_param_id, pco."ssID") = (ao.model_param_id, ao."ssID")
#     join model_param mp on mp.model_param_id = ao.model_param_id
#     where (simulation_param_uuid = '132323a4-3f93-4287-aca9-d18e84848e37'
#     and ( mp.variant='11' or mp.variant='12')
#     and mp.n_samples=1000000
#     and ss_n=1)

#     or (simulation_param_uuid = '132323a4-3f93-4287-aca9-d18e84848e37'
#     and mp.variant='0' 
#     and mp.n_samples=2000000
#     and ss_n=1);
#     '''




query = '''select mp.* from pattern_class_output pco
    join analytical_output ao on (pco.model_param_id, pco."ssID") = (ao.model_param_id, ao."ssID")
    join model_param mp on mp.model_param_id = ao.model_param_id
    where (simulation_param_uuid = '132323a4-3f93-4287-aca9-d18e84848e37'
    and ( mp.variant='11' or mp.variant='12')
    and mp.n_samples=1000000
    and ss_n=1)

    or (simulation_param_uuid = '132323a4-3f93-4287-aca9-d18e84848e37'
    and mp.variant='0' 
    and mp.n_samples=2000000
    and ss_n=1);
    '''
def df_from_query(query):

    credentials=f"postgresql://moliver:moliver@ld-rendres07.bc.ic.ac.uk/moliver"
    with psycopg2.connect(credentials) as conn:
        with conn.cursor() as cursor:
            df = pd.read_sql_query(query, conn)
    df = df.dropna(axis=1, how='all')
    df = df.drop(['model_param_id'], axis=1)
    # df = df.drop(model_param_dict.keys(), axis=1)
    df.set_index('parID', inplace=True)
    return df


df = df_from_query(query)
df

    #%%
    # Specify name of circuit and variant investigated
    # circuit_n='14'
    # parID=4187715
    # # nsr=0.01
    # old_variant='fitted7'
    # variant=f'{old_variant}_gaussian{parID}_nsr{nsr}' #variant='fitted7'
    # Specifiy number of parameter sets in parameterset file to be loaded
    # n_samples = 2000 #n_samples = 13700000


circuit_n=14;variant='2nd';n_species=6;n_samples = 1000000



model_param_dict = {'circuit_n': circuit_n, 'variant': variant, 'n_samples': n_samples}
#%%
result_df = query_modelParam_df_from_sql(model_param_dict)


# pickle.dump(result_df, open('file.txt', 'wb' ) )
# pickle.dump(result_df, open(modellingpath + '/3954/paper/input/gaussian_parameterfiles/df_circuit%s_variant%s_%rparametersets.pkl'%(circuit_n,variant,n_samples), 'wb' ) )
# pickle.dump(result_df, open(modellingpath + '/3954/paper/gaussian_parameterfiles/df_circuit%s_variant%s_%rparametersets.pkl'%(circuit_n,variant,n_samples), 'wb' ) )
# with open(modellingpath + '/3954/paper/input/gaussian_parameterfiles/df_circuit%r_variant%s_%rparametersets.pkl'%(circuit_n,variant,n_samples), 'wb' ) as f:
    # pickle.dump(object, result_df)


import psycopg2
credentials=f"postgresql://moliver:moliver@ld-rendres07.bc.ic.ac.uk/moliver"
with psycopg2.connect(credentials) as conn:
    with conn.cursor() as cursor:
        query = '''select ao.ss_n, mp."parID", ao."ssID", mp.variant,simulation_param_uuid,  ao.model_param_id,ao.system_class,pco.pattern_class_nogrowth from pattern_class_output pco
join analytical_output ao on (pco.model_param_id, pco."ssID") = (ao.model_param_id, ao."ssID")
join model_param mp on mp.model_param_id = ao.model_param_id
where (simulation_param_uuid = '132323a4-3f93-4287-aca9-d18e84848e37'
and ( mp.variant='11' or mp.variant='12')
and mp.n_samples=1000000
and ss_n=1)

or (simulation_param_uuid = '132323a4-3f93-4287-aca9-d18e84848e37'
and mp.variant='0' 
and mp.n_samples=2000000
and ss_n=1);
'''
        lsa_vs_numerical_df = pd.read_sql_query(query, conn)

lsa_vs_numerical_df
        