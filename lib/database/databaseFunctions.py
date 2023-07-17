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
import math as math
#user=moliver
#server=ld-rendres07.bc.ic.ac.uk
#database=moliver

credentials=f"postgresql://moliver:moliver@ld-rendres07.bc.ic.ac.uk/moliver"



    

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

    with psycopg2.connect(credentials) as conn:
        with conn.cursor() as cursor:
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
                    
            df.index.name = 'parID'
            print('Preparing to insert data into database')
            st = time.time()
            rows = df.to_sql('model_param', con=credentials, if_exists='append', index=True, index_label='parID', method=postgres_upsert)
            et = time.time()
            # get the execution time
            elapsed_time = et - st
            print('Execution time:', elapsed_time, 'seconds')   
            print('Inserted data into database')


def analyticalOutput_df_to_sql(lsa_df, circuit_n, variant, n_samples):

    with psycopg2.connect(credentials) as conn:
        with conn.cursor() as cursor:
            lsa_df['circuit_n'] = circuit_n
            lsa_df['variant'] = variant
            lsa_df['n_samples'] = n_samples
            lsa_df = lsa_df.rename_axis(['parID','ssID']).reset_index()

            lsa_df['model_param_id'] = lsa_df.apply(lambda row: f"{row['parID']}_circuit:{row['circuit_n']}_variant:{row['variant']}_samples:{row['n_samples']}", axis=1)

            lsa_df['maxeig'] = np.real(lsa_df['maxeig'] )
            print('prelambda')
            lsa_df['ss_list'] = lsa_df['ss_list'].apply(lambda x: x.tolist() if type(x) == 'numpy.ndarray' else [None])


            cursor.execute("SELECT column_name FROM information_schema.columns WHERE table_schema = 'public' AND table_name = 'analytical_output';")
            column_names = [row[0] for row in cursor]
            lsa_df = lsa_df[column_names]
            print(f'Preparing to insert analytical output data into database {circuit_n}, {variant}, {n_samples}')
            st = time.time()
            rows = lsa_df.to_sql('analytical_output', con=credentials, if_exists='append',  method=postgres_upsert,  index=False)
            et = time.time()
            elapsed_time = et - st
            print('Execution time:', elapsed_time, 'seconds')   
            print('Inserted data into database')

            return lsa_df




def simulationParam_to_sql(sim_dict):
    print(sim_dict)
    with psycopg2.connect(credentials) as conn:
        with conn.cursor() as cursor:
            table_name = 'simulation_param'
            id = generate_insert_uuid('simulation_param_id',table_name,cursor,conn)
            sim_dict['simulation_param_id'] = id
            print(sim_dict)

            df = pd.DataFrame(sim_dict, index=[sim_dict['simulation_param_id']])
            df.index.name = 'simulation_param_id'
            df = df.drop(['simulation_param_id'], axis=1)
            
            print('Preparing to insert data into database')
            st = time.time()
            df.to_sql(table_name, con=credentials, if_exists='append', index=True, index_label='simulation_param_id', method=postgres_upsert)

            et = time.time()
            # get the execution time
            elapsed_time = et - st
            print('Execution time:', elapsed_time, 'seconds')   
            print('Inserted data into database')

            conn.commit()

            

def simulationOutput_to_sql(sim_param_dict,model_param_dict,U_final_1D,U_record_1D, ssID=0):
    U_final_1D_list = np.array(U_final_1D).tolist()
    U_record_1D_list = np.array(U_record_1D).tolist()
    with psycopg2.connect(credentials) as conn:
        with conn.cursor() as cursor:    
            
            table_name = "simulation_param"
                        
                        # Build the query dynamically
            query = "SELECT simulation_param_id FROM {} WHERE ".format(table_name)
            conditions = []
            values = []
            for key, value in sim_param_dict.items():
                conditions.append('"{0}" = %s'.format(key))
                values.append(value)
            cursor.execute(query + ' AND '.join(conditions), values)

            result = cursor.fetchall()
            if len(result)>1:
                print('error!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            else:
                simulation_param_id = result[0]

            print(f"simulation_param_id:{simulation_param_id}")
            model_param_id =  f"{model_param_dict['parID']}_circuit:{model_param_dict['circuit_n']}_variant:{model_param_dict['variant']}_samples:{model_param_dict['n_samples']}"
            print(f"model_param_id:{model_param_id}")
            insert_query = 'INSERT INTO simulation_output ("simulation_param_id", "model_param_id", "ssID", "U_final_1D","U_record_1D") VALUES (%s, %s, %s,%s,%s) ON CONFLICT DO NOTHING'
            values = (simulation_param_id, model_param_id, ssID, U_final_1D_list,U_record_1D_list)
            cursor.execute(insert_query, values)
            conn.commit()

            print('simulation_output inserted')
            return model_param_id




def query_simulationOutput_from_sql(sim_param_dict,model_param_dict,query_column, ssID=0):
    with psycopg2.connect(credentials) as conn:
        with conn.cursor() as cursor:
            cursor.execute('SELECT "simulation_param_id" from simulation_param WHERE "L" = 50 and "dx" = 0.1;')
            conn.commit()
            simulation_param_id = cursor.fetchall()[0]
            print(f"simulation_param_id:{simulation_param_id}")
            model_param_id =  f"{model_param_dict['parID']}_circuit:{model_param_dict['circuit_n']}_variant:{model_param_dict['variant']}_samples:{model_param_dict['n_samples']}"
            print(f"model_param_id:{model_param_id}")
            
            insert_query = 'SELECT "U_final_1D" from simulation_output where "model_param_id"=(%s) and "simulation_param_id"=(%s) and "ssID"=(%s)'
            values = (model_param_id, simulation_param_id, ssID)
            cursor.execute(insert_query, values)
            simulationOutput = np.array(cursor.fetchall()[0][0],dtype=float)

            conn.commit()




            return simulationOutput




# %%
