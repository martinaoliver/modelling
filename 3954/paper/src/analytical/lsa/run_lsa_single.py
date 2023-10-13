#%%
#############

###paths#####
#############
import sys
import os

pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')
#############

from analytical.linear_stability_analysis import big_turing_analysis_df, detailed_turing_analysis_dict

from randomfunctions import *

import pickle

#######################
#########CODE##########
#######################







#%%
circuit_n='circuit14'
# variant='fitted7'
# variant='2nd'
# balance = 'balanced'#'balancedSemiBalanced'
parID = 41 #takes the first parameter set of the dataframe... can choose any
n_species=6 #number of molecular species in circuit_n (#Circuit2 has 6 molecular species)
# #turing 995206
# n_samples = 1000000#13700000
system_class='Stable'
# obtain a dictionary with some parameters to use in our analysis
# lsa_df= pickle.load( open(modellingpath + "/3954/paper/out/analytical/lsa_dataframes/all_dataframes/lsa_df_%s_variant%s_%rparametersets_%s.pkl"%(circuit_n,variant,n_samples, balance), "rb"))
# df= pickle.load( open(modellingpath + "/3954/paper/input/balanced_parameterfiles/df_%s_variant%s_%rparametersets_%s.pkl"%(circuit_n,variant,n_samples, balance), "rb"))
df= pickle.load( open(modellingpath + "/3954/paper/out/analytical/lsa_dataframes/all_dataframes/lsa_df_circuit14_variantfitted7_gaussian4187715_nsr0.01_2000parametersets.pkl", "rb"))
# df= pickle.load( open(modellingpath + "/3954/paper/input/fitted_parameterfiles/df_%s_variant%s_%rparametersets.pkl"%(circuit_n,variant,n_samples), "rb"))
# df= pickle.load( open(modellingpath + "/3954/paper/input/fitted_parameterfiles/df_%s_variant%s_%rparametersets.pkl"%(circuit_n,variant,n_param_sets), "rb"))
# n_analysed_samples = 10
# df = df.iloc[parID-1:parID+1]

par_dict = df.loc[parID,0].to_dict()
#Run analysis on 1M parameter sets
# output_df = big_turing_analysis_df(df,circuit_n,variant,n_samples, n_species,print_parID=False, tqdm_disable=False)
out = detailed_turing_analysis_dict(par_dict, circuit_n,n_species,top_dispersion=1000,calculate_unstable=False,steadystate=False)
# plot_all_dispersion(out[4][0],n_species, crop=100, top=300)
# plt.show()
#%%


# plot_highest_dispersion_hopf(out[4][0],crop = 300, top = 2000)
plot_highest_dispersion(out[4][0],crop = 100, top = 2000)
plt.title(f'{system_class}')
plt.tight_layout()
plt.savefig(modellingpath + "/3954/paper/out/analytical/pyPlots/dispersion_relation/%s_dispersion_circuit14_variantfitted7_gaussian4187715_nsr0.01_parID41.pdf"%system_class)
# plt.show()

# print(out[3])



# %%
