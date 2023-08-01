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

def general_query(query):
    with psycopg2.connect(credentials) as conn:
        with conn.cursor() as cursor:
            cursor.execute(query)
            return cursor

    

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




def insert_simulationParam_to_sql(model_param_dict, table_name='simulation_param'):
    with psycopg2.connect(credentials) as conn:
        with conn.cursor() as cursor:
            try: # Prepare the SQL query to insert data into the table
                columns = ", ".join([f"\"{key}\"" for key in model_param_dict.keys()])
                values = ", ".join([f"'{val}'" for val in model_param_dict.values()])
                query = f"INSERT INTO {table_name} ({columns}) VALUES ({values}) ;"
                print(query)
        
                cursor.execute(query)
                conn.commit()
                print("Data inserted successfully!")
                
            except (Exception, psycopg2.Error) as error:
                print("Error while inserting data:", error)

            


def insert_simulationOutput_to_sql(sim_param_dict,model_param_dict,U_final,U_record, ssID=0, dimensions='2D'):
    U_final_1D_list = np.array(U_final).tolist()
    U_record_1D_list = np.array(U_record).tolist()
    with psycopg2.connect(credentials) as conn:
        with conn.cursor() as cursor:    
            
            table_name = "simulation_param"
                        
                        # Build the query dynamically
            query = "SELECT simulation_param_uuid FROM {} WHERE ".format(table_name)
            conditions = []
            values = []
            for key, value in sim_param_dict.items():
                conditions.append('"{0}" = %s'.format(key))
                values.append(value)
            cursor.execute(query + ' AND '.join(conditions), values)

            result = cursor.fetchall()
            if len(result)>1:
                print('error!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
                return 'error'
            else:
                simulation_param_uuid = result[0]

            print(f"simulation_param_uuid:{simulation_param_uuid}")
            model_param_id =  f"{model_param_dict['parID']}_circuit:{model_param_dict['circuit_n']}_variant:{model_param_dict['variant']}_samples:{model_param_dict['n_samples']}"
            print(f"model_param_id:{model_param_id}")
            if dimensions=='1D':
                print('1D')
                insert_query = 'INSERT INTO simulation_output ("simulation_param_uuid", "model_param_id", "ssID", "U_final_1D","U_record_1D") VALUES (%s, %s, %s,%s,%s)  ON CONFLICT ("simulation_param_uuid", "model_param_id","ssID")  DO UPDATE SET "U_final_1D" = EXCLUDED."U_final_1D", "U_record_1D" = EXCLUDED."U_record_1D";'
            elif dimensions=='2D':
                print('2D')
                insert_query = 'INSERT INTO simulation_output ("simulation_param_uuid", "model_param_id", "ssID", "U_final_2D","U_record_2D") VALUES (%s, %s, %s,%s,%s)'
            values = (simulation_param_uuid, model_param_id, ssID, U_final_1D_list,U_record_1D_list)
            cursor.execute(insert_query, values)
            conn.commit()

            print('simulation_output inserted')
            return model_param_id
        
        
def query_modelParam_df_from_sql( model_param_dict):
    with psycopg2.connect(credentials) as conn:
        with conn.cursor() as cursor:

            # Build the SQL query
            columns = ", ".join(model_param_dict.keys())
            values = ", ".join([f"'{val}'" for val in model_param_dict.values()])
            query = f"SELECT * FROM model_param WHERE ({columns}) = ({values});"

            # Fetch the data into a pandas DataFrame
            df = pd.read_sql_query(query, conn)
    df = df.dropna(axis=1, how='all')
    df = df.drop(['model_param_id'], axis=1)
    df = df.drop(model_param_dict.keys(), axis=1)
    df.set_index('parID', inplace=True)
    return df



def query_analyticalOutput_df_from_sql(model_param_dict,limit='None'):
    with psycopg2.connect(credentials) as conn:
        with conn.cursor() as cursor:
        # Build the SQL query dynamically based on the provided dictionaries
            #measure time
            query = """
            SELECT *
            FROM analytical_output ao
            JOIN model_param mp on mp.model_param_id = ao.model_param_id
            WHERE 1=1
            """
            
            # Add filters for model_params
            for key, value in model_param_dict.items():
                query += f"AND mp.{key} = {value}\n"
            if limit != 'None':
                licensequery += f'LIMIT {limit}'

            print('query',query)



            df = pd.read_sql_query(query, conn)

            df = df.dropna(axis=1, how='all')
            df = df.drop(['model_param_id'], axis=1)
            df = df.drop(model_param_dict.keys(), axis=1)
            df.set_index('parID', inplace=True)


            return df


def query_simulationOutput_from_sql(sim_param_dict,model_param_dict,query_column, ssID=0):
    with psycopg2.connect(credentials) as conn:
        with conn.cursor() as cursor:
            table_name = "simulation_param"
                        
                        # Build the query dynamically
            query = "SELECT simulation_param_uuid FROM {} WHERE ".format(table_name)
            conditions = []
            values = []
            for key, value in sim_param_dict.items():
                conditions.append('"{0}" = %s'.format(key))
                values.append(value)
            cursor.execute(query + ' AND '.join(conditions), values)

            result = cursor.fetchall()
            if len(result)>1:
                print('error!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
                simulation_param_uuid = result[0]

            else:
                simulation_param_uuid = result[0]
            print(f"simulation_param_uuid:{simulation_param_uuid}")
            model_param_id =  f"{model_param_dict['parID']}_circuit:{model_param_dict['circuit_n']}_variant:{model_param_dict['variant']}_samples:{model_param_dict['n_samples']}"
            print(f"model_param_id:{model_param_id}")
            
            # if query_column == '
            insert_query = f'SELECT "{query_column}" from simulation_output where "model_param_id"=(%s) and "simulation_param_uuid"=(%s) and "ssID"=(%s)'
            values = (model_param_id, simulation_param_uuid, ssID)
            cursor.execute(insert_query, values)
            simulationOutput = np.array(cursor.fetchall()[0][0],dtype=float)

            conn.commit()




            return simulationOutput

