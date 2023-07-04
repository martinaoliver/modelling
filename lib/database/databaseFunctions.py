#############
###paths#####
#############
import sys
import os


import psycopg2
import os
import time
import numpy as np
import random
from datetime import datetime
import string
import pandas as pd
# credentials=f"postgresql://moliver:{password}@ld-rendres07.bc.ic.ac.uk/moliver"
credentials=f"postgresql://moliver:moliver@ld-rendres07.bc.ic.ac.uk/moliver"

#user=moliver
#server=ld-rendres07.bc.ic.ac.uk
#database=moliver

#this function allows us to "on conflict - update"
def postgres_upsert(table, conn, keys, data_iter):
    from sqlalchemy.dialects.postgresql import insert

    data = [dict(zip(keys, row)) for row in data_iter]

    insert_statement = insert(table.table).values(data)
    upsert_statement = insert_statement.on_conflict_do_update(
        constraint=f"{table.table.name}_pkey",
        set_={c.key: c for c in insert_statement.excluded},
    )
    conn.execute(upsert_statement)

def generate_insert_uuid(id_name,table_name,cursor,conn):# Generate a unique key for simID
    #generate random
    random.seed(datetime.now().timestamp())
    id = ''.join(random.choices(string.ascii_lowercase, k=20))

    # Check if the simID already exists in Table3
    check_unique_query = f'SELECT COUNT(*) FROM {table_name} WHERE {table_name}."{id_name}" = \'{id}\''
    cursor.execute(check_unique_query)
    existing_count = cursor.fetchone()[0]
    # If the simID already exists, regenerate the key
    while existing_count > 0:
        id = ''.join(random.choices(string.ascii_lowercase, k=10))
        cursor.execute(check_unique_query)
        existing_count = cursor.fetchone()[0]

    return id


# df has:
# - columns: parID, circuit_n, var, ...parameters
def modelParam_df_to_sql(df, circuit_n, variant, n_samples):

    df['circuit_n'] = circuit_n
    df['variant'] = variant
    df['n_samples'] = n_samples
    conn = psycopg2.connect(credentials)
    cursor= conn.cursor()
    for column in df.columns:
        # cur.execute("SELECT f_add_col('public.%s', '%s', 'numeric');"%('model_param',column))
        # cur.execute('ALTER TABLE %s ADD COLUMN "%s" numeric' % ('model_param', column))
        try:
            cursor.execute('ALTER TABLE %s ADD COLUMN "%s" numeric' % ('model_param', column))

            conn.commit()
            print(f'New columns added: {column}')
        except psycopg2.errors.DuplicateColumn:
            conn.rollback()
            # print('failed to commit: ')
            
    conn.close()
    df.index.name = 'parID'
    print('Preparing to insert data into database')
    st = time.time()
    rows = df.to_sql('model_param', con=credentials, if_exists='append', index=True, index_label='parID', method=postgres_upsert)
    et = time.time()
    # get the execution time
    elapsed_time = et - st
    print('Execution time:', elapsed_time, 'seconds')   
    print('Inserted data into database')
    return


# df has:
# - columns: parID, circuit_n, var, ...parameters, analyticalResult
def analyticalOutput_df_to_sql(lsa_df, circuit_n, variant, n_samples):
    # TODO remove modelParameters from df
    # TODO will insert 


    lsa_df['circuit_n'] = circuit_n
    lsa_df['variant'] = variant
    lsa_df['n_samples'] = n_samples

    lsa_df['maxeig'] = np.real(lsa_df['maxeig'] )
    lsa_df['ss_list'] = lsa_df['ss_list'].apply(lambda x: x.tolist())

    conn = psycopg2.connect(credentials)
    cursor= conn.cursor()
    cursor.execute("SELECT column_name FROM information_schema.columns WHERE table_schema = 'public' AND table_name = 'model_param';")
    column_names = [row[0] for row in cur]
    [column_names.remove(x) for x in ['parID','circuit_n', 'variant', 'n_samples']]
    lsa_df = lsa_df.rename_axis(['parID','ssID']).reset_index()
    lsa_df = lsa_df.drop(columns = column_names)
    print(lsa_df)
    type_dict = {'parID':'numeric','circuit_n':'text','n_samples':'numeric','variant':'text','ssID': 'numeric', 'ss_n': 'numeric',  'ss_list': 'numeric[]', 'ss_class':'text', 'system_class':'text', 'maxeig':'numeric', 'estimated_wvl': 'numeric', 'complex_dispersion': 'bool'} 
# serialized_data = pickle.dumps(my_object)
    for column in lsa_df.columns:

        try:
            cursor.execute('ALTER TABLE %s ADD COLUMN "%s" %s' % ('analytical_output', column, type_dict[column]))
            conn.commit()
            print(f'New columns added: {column}')
        except psycopg2.errors.DuplicateColumn:
            conn.rollback()
            # print('failed to commit: ')
            
    conn.close()
    
    lsa_df.index.name = 'parID'
    print(lsa_df)
    print('Preparing to insert data into database')
    st = time.time()
    rows = lsa_df.to_sql('analytical_output', con=credentials, if_exists='append',  method=postgres_upsert,  index=False)
    et = time.time()
    elapsed_time = et - st
    print('Execution time:', elapsed_time, 'seconds')   
    print('Inserted data into database')

    return






def simulationParam_to_sql(sim_dict):
    print(sim_dict)
    conn = psycopg2.connect(credentials)
    cursor= conn.cursor()
    table_name = 'simulation_param'
    id = generate_insert_uuid('simID',table_name,cursor,conn)
    sim_dict['simID'] = id
    print(sim_dict)

    df = pd.DataFrame(sim_dict, index=[sim_dict['simID']])
    df.index.name = 'simID'
    df = df.drop(['simID'], axis=1)
    
    print('Preparing to insert data into database')
    st = time.time()
    df.to_sql(table_name, con=credentials, if_exists='append', index=True, index_label='simID', method=postgres_upsert)

    et = time.time()
    # get the execution time
    elapsed_time = et - st
    print('Execution time:', elapsed_time, 'seconds')   
    print('Inserted data into database')

    conn.commit()

    cursor.close()
    conn.close()
    









def saveNumericalResult(parID, simID, var, circuit_n, snapshot, timeseries):
    # TODO will insert 
    return



def queryAnalyticalResult(parID, circuit_n, var):
    # TODO will query 
    return df


def queryNumericalResult(parID, simID, circuit_n, var):
    # TODO will query 
    return snapshot, timeseries



# %%
