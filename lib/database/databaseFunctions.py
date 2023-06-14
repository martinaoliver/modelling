#############
###paths#####
#############
import sys
import os


import psycopg2
import os
import pickle
from tqdm import tqdm
import time
import numpy as np

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


# df has:
# - columns: parID, circuit_n, var, ...parameters
def modelParam_df_to_sql(df, circuit_n, variant, n_samples):

    df['circuit_n'] = circuit_n
    df['variant'] = variant
    df['n_samples'] = n_samples
    conn = psycopg2.connect(credentials)
    cur = conn.cursor()
    for column in df.columns:
        # cur.execute("SELECT f_add_col('public.%s', '%s', 'numeric');"%('model_param',column))
        # cur.execute('ALTER TABLE %s ADD COLUMN "%s" numeric' % ('model_param', column))
        try:
            cur.execute('ALTER TABLE %s ADD COLUMN "%s" numeric' % ('model_param', column))

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
    cur = conn.cursor()
    cur.execute("SELECT column_name FROM information_schema.columns WHERE table_schema = 'public' AND table_name = 'model_param';")
    column_names = [row[0] for row in cur]
    [column_names.remove(x) for x in ['parID','circuit_n', 'variant', 'n_samples']]
    lsa_df = lsa_df.rename_axis(['parID','ssID']).reset_index()
    lsa_df = lsa_df.drop(columns = column_names)
    print(lsa_df)
    type_dict = {'parID':'numeric','circuit_n':'text','n_samples':'numeric','variant':'text','ssID': 'numeric', 'ss_n': 'numeric',  'ss_list': 'numeric[]', 'ss_class':'text', 'system_class':'text', 'maxeig':'numeric', 'estimated_wvl': 'numeric', 'complex_dispersion': 'bool'} 
# serialized_data = pickle.dumps(my_object)
    for column in lsa_df.columns:

        try:
            cur.execute('ALTER TABLE %s ADD COLUMN "%s" %s' % ('analytical_output', column, type_dict[column]))
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



# df has:
# - columns: parID, circuit_n, var, ...parameters
def saveSimParameters(sim_dict):
    simID = 'simdictstring'
    # TODO will save simParameters to database
    return




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
