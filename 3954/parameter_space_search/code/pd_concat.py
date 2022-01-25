import pandas as pd
import pickle


circuit_n=2
variant= 0

# Specifiy number of parameter sets in parameterset file to be loaded
n_param_sets = 1000000


batch_indices = list(range(0, 1000000,10000 ))

my_data = {}

# Load all batch dataframes
for start_batch_index in batch_indices:
    print(start_batch_index)
    my_data[start_batch_index] = pickle.load(open('../results/output_dataframes/lsa_df_circuit%r_variant%r_%rparametersets_batch%r.pkl'%(circuit_n,variant,n_param_sets,start_batch_index), "rb" ) )

# Join all batch results to large results dataframe
results_df = pd.concat(my_data.values(), ignore_index=False)

# Pickle and save results dataframe
tupled_index =  [tuple(l) for l in results_df.index]
multi_index = pd.MultiIndex.from_tuples(tupled_index)
results_df = results_df.set_index(multi_index)
pickle.dump(results_df, open('../results/output_dataframes/lsa_df_circuit%r_variant%r_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), 'wb'))
