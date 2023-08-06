#%%#############

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








circuit_n='circuit14'
# variant='fitted7'
variant='2nd'
balance = 'balanced'#'balancedSemiBalanced'
parID = 999961 #takes the first parameter set of the dataframe... can choose any
n_species=6 #number of molecular species in circuit_n (#Circuit2 has 6 molecular species)

n_samples = 1000000#13700000

# obtain a dictionary with some parameters to use in our analysis
lsa_df= pickle.load( open(modellingpath + "/3954/paper/out/analytical/lsa_dataframes/all_dataframes/lsa_df_%s_variant%s_%rparametersets_%s.pkl"%(circuit_n,variant,n_samples, balance), "rb"))
df= pickle.load( open(modellingpath + "/3954/paper/input/balanced_parameterfiles/df_%s_variant%s_%rparametersets_%s.pkl"%(circuit_n,variant,n_samples, balance), "rb"))
# df= pickle.load( open(modellingpath + "/3954/paper/input/fitted_parameterfiles/df_%s_variant%s_%rparametersets.pkl"%(circuit_n,variant,n_samples), "rb"))
# df= pickle.load( open(modellingpath + "/3954/paper/input/fitted_parameterfiles/df_%s_variant%s_%rparametersets.pkl"%(circuit_n,variant,n_param_sets), "rb"))
# n_analysed_samples = 10
# df = df.iloc[parID-1:parID+1]
par_dict = df.loc[parID].to_dict()
#Run analysis on 1M parameter sets
# output_df = big_turing_analysis_df(df,circuit_n,variant,n_samples, n_species,print_parID=False, tqdm_disable=False)
out = detailed_turing_analysis_dict(par_dict, circuit_n,n_species,top_dispersion=1000,calculate_unstable=False,steadystate=False)
# plot_all_dispersion(out[4][0],n_species, crop=100, top=300)
# plt.show()
#%%
def plot_highest_dispersion_hopf(eigenvalues,crop = 1000, top = 5000, L=100):
    wvn_list = np.array(list(range(0,top+1)))*np.pi/L
    # wvn_list = np.array(list(range(0,5000+1)))*np.pi/100

    plt.plot(wvn_list[:crop], eigenvalues.real[:crop,[-1]], label='Real highest eigenvalue', c='k')
    plt.plot(wvn_list[1:crop], eigenvalues.imag[1:crop,[-1]], linestyle = '--', label = 'Imaginary highest eigenvalue', c='k')
    plt.plot(wvn_list[1:crop], eigenvalues.imag[1:crop,[-2]], linestyle = '--', label = 'Imaginary highest eigenvalue', c='k')

    plt.legend()
    plt.xlabel('Wavenumber')
    plt.ylabel('Highest eigenvalue')
    plt.axhline(y=0, color='k', linestyle='-', linewidth = 0.1)
    plt.grid()
    plt.tight_layout()

plot_highest_dispersion(out[4][3],crop = 300, top = 2000)
plt.title(f'Turing instability: {circuit_n}, {variant}, parID {parID}')
plt.tight_layout()
plt.savefig(modellingpath + "/3954/paper/out/analytical/pyPlots/dispersion_relation/unstable_dispersion.pdf")
plt.show()

# print(out[3])



# %%
