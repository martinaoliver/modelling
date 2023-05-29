#############
###paths#####
#############
import sys
import os

pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')
#############




import psycopg2
import os
import pickle

password = os.environ['DB_PASSWORD']
credentials=f"postgresql://marti:{password}@vps.dcotta.eu:5432/marti_phd"



# df has:
# - columns: parID, circuit_n, var, ...parameters
def saveModelParameters(df, circuit_n, variant):

    df['circuit_n'] = circuit_n
    df['variant'] = variant
    conn = psycopg2.connect(credentials)
    cur = conn.cursor()
    for column in df.columns:
        cur.execute('ALTER TABLE %s ADD COLUMN if not exists "%s" numeric' % ('model_param', column))
        conn.commit()
    conn.close()

    lhs_df.index.name = 'parID'
    # print(lhs_df)
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


#%%
# Specify name of circuit and variant investigated
circuit_n='circuit14'
variant='2nd'
# Specifiy number of parameter sets in parameterset file to be loaded
n_param_sets = 1000000

print(f'Circuit:{circuit_n}, Variant:{variant}')
lhs_df = pickle.load( open(modellingpath + '/3954/paper/input/lhs_parameterfiles/df_%s_variant%s_%rparametersets.pkl'%(circuit_n,variant,100000), "rb"))

saveModelParameters(lhs_df, circuit_n, variant)

# %%
