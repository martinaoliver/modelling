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
import pickle
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
print('heehehe')






circuit_n='turinghill'
variant= 0
n_species=2
# Specifiy number of parameter sets in parameterset file to be loaded
n_param_sets = 2000000





df= pickle.load( open(modellingpath + '/growth/out/analytical/lsa_dataframes/lsa_df_%s_variant%r_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), "rb"))
df = df.loc[df['ss_n']==1]# df= pickle.load( open(modellingpath + f'/growth/out/analytical/lsa_dataframes/lsa_df_turinghill_variant{variant}_{n_param_sets}parametersets_seed0.pkl', "rb"))
print(f'Variant {variant}')
print(df['system_class'].value_counts())
#%%
#cut df
cutDf = True
if cutDf == True:
    cropValue = 2000
    cuttedDf = df.iloc[:cropValue]
    pickle.dump( cuttedDf, open(modellingpath + '/growth/out/analytical/lsa_dataframes/lsa_df_%s_variant%r_%rparametersets.pkl'%(circuit_n,variant,cropValue), "wb" ) )
#%%
#plot pieChart
def pieChart_lsa(valueCounts_dict,title,log=True):
    colors=['grey','grey','grey','grey','peachpuff','coral','coral','coral','coral','coral','grey','coral','coral']
    labels = []
    sizes = []
    
    for x, y in valueCounts_dict.items():
        labels.append(x)
        sizes.append(y)
    if log==True:
        sizes = np.log(sizes)
    plt.pie(sizes,colors=colors, labels=labels)
    plt.axis('equal')
    plt.title(title)
    plt.show()

valueCounts_dict = dict(df['system_class'].value_counts())
title = f'{circuit_n} Variant {variant}'
pieChart_lsa(valueCounts_dict,title)
#%%


#values for which complex dispersion = true
complex_df = df[df['complex_dispersion']==True]
complex_df.index  = complex_df.index.droplevel(-1)

#values that have instabilities
saveInstabilitiesComplex = True
if saveInstabilitiesComplex ==True:
    instabilities = ['turing I', 'turing II', 'turing I hopf', 'turing I oscillatory', 'turing II hopf','hopf', 'turing semi-hopf']  
    instabilities_df = df.loc[df['system_class'].isin(instabilities)]
    # instabilities_df.index  = instabilities_df.index.droplevel(-1)
    print(instabilities_df['system_class'].value_counts())
    pickle.dump( instabilities_df, open(modellingpath + '/growth/out/analytical/instability/instability_df_%s_variant%r_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), "wb" ) )

    instabilityComplex_df = pd.concat([complex_df, instabilities_df])
    pickle.dump( instabilityComplex_df, open(modellingpath + '/growth/out/analytical/instabilityComplex/instabilityComplex_df_%s_variant%r_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), "wb" ) )

    print(instabilityComplex_df)

# #values that have turing
saveTuring = True
if saveTuring == True:
    turingStates = ['turing I','turing I oscillatory']  
    turing_df = df.loc[df['system_class'].isin(turingStates)]
    # turing_df = turing_df.xs(0, level=1)
    # turing_df.index  = turing_df.index.droplevel(-1)
    # turing_df = turing_df.loc[turing_df['ss_n']==1]
    # turing_df = turing_df.sort_values(by=['maxeig'],  ascending=False)
    pickle.dump( turing_df, open(modellingpath + '/growth/out/analytical/turing/turing_df_%s_variant%r_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), "wb" ) )
    print(turing_df)




# %%

