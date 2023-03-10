#############################
#IMPORTS#
#############################
import sys
import os
pwd = os.getcwd()
root = pwd.rpartition("mo2016")[0] + pwd.rpartition("mo2016")[1] #/Volumes/mo2016/ or '/Users/mo2016/' or '/rds/general/mo2016/'

# print(root)
if root == '/Users/mo2016':
    modelling_ephemeral = '/Volumes/mo2016/ephemeral/Documents/modelling'
    modelling_home = '/Volumes/mo2016/home/Documents/modelling'
    modelling_local = root + '/Documents/modelling'
    import matplotlib as mpl
    mpl.use('tkagg')

if root == '/Volumes/mo2016' or  root=='/rds/general': #'/rds/general' or root=='/Volumes':
        modelling_ephemeral = root + '/ephemeral/Documents/modelling'
        modelling_home = root  + '/home/Documents/modelling'
        modelling_local = modelling_home

modulepath = modelling_local + '/3954/modules'
sys.path.append(modulepath)

import sys
import time
import pickle
from scipy.integrate import odeint
import matplotlib.pyplot as plt

from linear_stability_analysis import *
from class_subcircuit_eq import *

#######################
#########CODE##########
#######################

circuit_n = 2 #ID of circuit we want to analyse
#(parameter sets provided correspond to circuit2 which is the one currently being implemented experimentally)
parID = 0 #takes the first parameter set of the dataframe... can choose any
n_species=6 #number of molecular species in circuit_n (#Circuit2 has 6 molecular species)
variant=0
n_param_sets=1000000
#obtain a dictionary with some parameters to use in our analysis
df= pickle.load( open(modelling_home + '/3954/parameter_space_search/results/output_dataframes/lsa_df_circuit%r_variant%r_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), "rb" ) )
# df= pickle.load( open(modelling_home + '/3954/parameter_space_search/parameterfiles/df_circuit%r_variant%r_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), "rb" ) )
# turing_df = turing_df.droplevel(level=1)
turing_df = df.loc[df['system_class'] == 'turing I']
turing_df = turing_df.droplevel(level=1)

stable = df.loc[df['system_class'] == 'simple stable']
stable = stable.droplevel(level=1)
turing_df = turing_df.reset_index()
stable = stable.reset_index()

print(turing_df)
print(stable)
# print(df.loc[653])
# print(par_dict)
 #Run analysis on 1M parameter sets
# output_df = detailed_turing_analysis_dict(par_dict, circuit_n,n_species,top_dispersion=5000,calculate_unstable=False)
# print(output_df)
#Runs ODE

%matplotlib inline

def loguniform(low=-3, high=3, size=None):
    return (10) ** (np.random.uniform(low, high, size))

n_parIDs = 50
responses_turing= {}
responses_prime_turing= {}

responses_nonturing= {}
responses_prime_nonturing= {}

responses_list_nonturing = []
responses_prime_list_nonturing = []
km_list =np.logspace(-3,3,50)
# fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(18, 6))
for parID in turing_df.index[:n_parIDs]:
    par_dict = turing_df.loc[parID].to_dict()
    U_ss = []
    for km in km_list:
        par_dict['kce']=km
        T=300
        t = np.linspace(0, T, T*1000)
        initial_conditions = [0]
        subcircuit = subcircuitAC_1(par_dict)
        sol = odeint(subcircuit.ddt, initial_conditions, t)

        Ustar_ODE=sol[-1]
        U_ss.append(sol[-1,0])
    U_ss_prime = np.gradient(U_ss, km_list)
    ax[1].plot(km_list, U_ss_prime, color='k', label='Turing')
    ax[0].plot(km_list, U_ss,color='k', label='Turing')
    responses_turing[parID]=U_ss
    responses_prime_turing[parID] = U_ss_prime

print('turing_done')
for parID in stable.index[:n_parIDs]:
    par_dict = stable.loc[parID].to_dict()
    U_ss = []
    for km in km_list:
        par_dict['kce']=km
        T=300
        t = np.linspace(0, T, T*1000)
        initial_conditions = [0]
        subcircuit = subcircuitAC_1(par_dict)
        sol = odeint(subcircuit.ddt, initial_conditions, t)

        Ustar_ODE=sol[-1]
        U_ss.append(sol[-1,0])
    U_ss_prime = np.gradient(U_ss, km_list)
    ax[1].plot(km_list, U_ss_prime, color='r', label='Non Turing')
    ax[0].plot(km_list, U_ss,color='r', label='Non Turing')
    responses_nonturing[parID]=U_ss
    responses_prime_nonturing[parID] = U_ss_prime




###Dose responses 
fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(18, 6))
mylabel = ['Turing', 'Non Turing']
ax[0].plot(km_list, responses_turing[1], color='k',alpha=1)
ax[1].plot(km_list, responses_prime_turing[1],color='k', label=mylabel[0],alpha=1)

mylabel = ["_nolegend_","_nolegend_"]
for i in range(2):
    ax[i].set_xscale("log")
    # ax[i].set_yscale("log")
    ax[i].set_xlabel('Relative [ATC]')
    ax[i].legend()

ax[0].set_ylabel('Relative RFP')
ax[1].set_ylabel('Derivative (Relative RFP)')



alpha=0.7

TC, NTC = 'mediumaquamarine', 'lightcoral'
###Dose responses 
fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(18, 6))
mylabel = ['Turing', 'Non Turing']
for parID in stable.index[:n_parIDs]:
    ax[0].plot(km_list, responses_turing[parID], color=TC, label=mylabel[0],alpha=alpha)
    ax[1].plot(km_list, responses_prime_turing[parID],color=TC, label=mylabel[0],alpha=alpha)

    ax[0].plot(km_list, responses_nonturing[parID], color=NTC, label=mylabel[1],alpha=alpha)
    ax[1].plot(km_list, responses_prime_nonturing[parID],color=NTC, label=mylabel[1],alpha=alpha)
    mylabel = ["_nolegend_","_nolegend_"]
for i in range(2):
    ax[i].set_xscale("log")
    # ax[i].set_yscale("log")
    ax[i].set_xlabel('Relative [ATC]')
    ax[i].legend()

ax[0].set_ylabel('Relative RFP')
ax[1].set_ylabel('Derivative (Relative RFP)')


####Dose responses averaged
fig,ax = plt.subplots(nrows=1, ncols=2, figsize=(18, 6))

sum_responses_turing = np.sum([responses_turing[item] for item in responses_turing], axis=0)/n_parIDs
sum_responses_nonturing = np.sum([responses_nonturing[item] for item in responses_nonturing], axis=0)/n_parIDs
ax[0].plot(km_list,sum_responses_turing,color=TC, label='Turing')
ax[0].plot(km_list,sum_responses_nonturing,color=NTC, label='Non Turing')

sum_responses_prime_turing = np.sum([responses_prime_turing[item] for item in responses_prime_turing], axis=0)/n_parIDs
sum_responses_prime_nonturing = np.sum([responses_prime_nonturing[item] for item in responses_prime_nonturing], axis=0)/n_parIDs
ax[1].plot(km_list,sum_responses_prime_turing,color=TC, label='Turing')
ax[1].plot(km_list,sum_responses_prime_nonturing,color=NTC, label='Non Turing')

for i in range(2):
    ax[i].set_xscale("log")
    ax[i].legend()
    ax[i].set_xlabel('Relative [ATC]')

ax[0].set_ylabel('Relative RFP')
ax[1].set_ylabel('Derivative (Relative RFP)')


####Dynamic ranges of dose response

def dynamic_range(U_ss):
    return U_ss[-1]-U_ss[0]

turing_dr = {}
non_turing_dr = {}
for parID in stable.index[:n_parIDs]:
    U_ss = responses_turing[parID]
    turing_dr[parID] = dynamic_range(U_ss)
    U_ss = responses_nonturing[parID]
    non_turing_dr[parID] = dynamic_range(U_ss)

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 4))
nbins=50
alpha=0.5
ax.hist(turing_dr.values(), bins=nbins, color=TC, label='Turing', alpha=alpha)
ax.hist(non_turing_dr.values(), bins=nbins, color=NTC,  label='Non Turing', alpha=alpha)
ax.legend()
ax.set_xlabel('Dynamic range')
ax.set_ylabel('Frequency')
