import sys
import os

pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')

from database.databaseFunctions import *
import pickle
#############


# #%%
# Specify name of circuit and variant investigated

circuit_n='turinghill'
variant=9

# Specifiy number of parameter sets in parameterset file to be loaded
n_samples = 2000000

print(f'Circuit:{circuit_n}, Variant:{variant}')
lsa_df = pickle.load( open(modellingpath + '/growth/out/analytical/lsa_dataframes/lsa_df_%s_variant%s_%rparametersets.pkl'%(circuit_n,variant,n_samples), "rb"))
# lsa_df = lsa_df.iloc[:10000]
analyticalOutput_df_to_sql(lsa_df, circuit_n, variant, n_samples)
# credentials=f"postgresql://moliver:moliver@ld-rendres07.bc.ic.ac.uk/moliver"
# conn = psycopg2.connect(credentials)
# cur = conn.cursor()
# cur.execute("SELECT column_name FROM information_schema.columns WHERE table_schema = 'public' AND table_name = 'model_param';")
# column_names = [row[0] for row in cur]
# [column_names.remove(x) for x in ['parID','circuit_n', 'variant', 'n_samples']]
# lsa_df = lsa_df.rename_axis(['parID','ssID']).reset_index()
# lsa_df.drop(columns = column_names)