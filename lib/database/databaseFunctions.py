#############
###paths#####
#############
import sys
import os


import psycopg2
import os
import pickle



# credentials=f"postgresql://moliver:{password}@ld-rendres07.bc.ic.ac.uk/moliver"
credentials=f"postgresql://moliver:moliver@ld-rendres07.bc.ic.ac.uk/moliver"

#user=moliver
#server=ld-rendres07.bc.ic.ac.uk
#database=moliver

# df has:
# - columns: parID, circuit_n, var, ...parameters
def saveModelParameters(df, circuit_n, variant):

    df['circuit_n'] = circuit_n
    df['variant'] = variant
    conn = psycopg2.connect(credentials)
    cur = conn.cursor()
    for column in df.columns:
        # cur.execute("SELECT f_add_col('public.%s', '%s', 'numeric');"%('model_param',column))
        # cur.execute('ALTER TABLE %s ADD COLUMN "%s" numeric' % ('model_param', column))
        try:
            cur.execute('ALTER TABLE %s ADD COLUMN "%s" numeric' % ('model_param', column))

            conn.commit()
        except psycopg2.errors.DuplicateColumn:
            conn.rollback()
            # print('failed to commit: ')
            
    conn.close()

    df.index.name = 'parID'
    # print(lhs_df)
    print('aaaa')
    rows = df.to_sql('model_param', con=credentials, if_exists='append', index=True, index_label='parID')
    # print('inserted_rows: ', rows)
    return


# df has:
# - columns: parID, circuit_n, var, ...parameters
def saveSimParameters(sim_dict):
    simID = 'simdictstring'
    # TODO will save simParameters to database
    return


# df has:
# - columns: parID, circuit_n, var, ...parameters, analyticalResult
def saveAnalyticalResult(analyticalResultDf):
    # TODO remove modelParameters from df
    # TODO will insert 
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
