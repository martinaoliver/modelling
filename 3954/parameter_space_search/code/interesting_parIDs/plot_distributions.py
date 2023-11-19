import sys
import os
pwd = os.getcwd()
root = pwd.split("home", 1)[0]
modelling_home = root + 'home/Documents/modelling'
modelling_ephemeral = root + 'ephemeral/Documents/modelling'
modulepath = modelling_home + '/3954/modules'
sys.path.append(modulepath)

from numerical_solvers_variableboundary import *
import pickle
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import numpy as np
# general_df = pickle.load(open('bullseye_df.pkl', "rb"))
general_df = pickle.load(open('interesting_parIDs.pkl', "rb"))

df_concat = general_df

#Distribution of Turing I Monostable kinetic parameters
# turing_monostable_df, loguniform_turing_monostable_df, joint_turing_monostable_df = open_dfs('circuit2_turingI_monostable_df')
# lenght_turing = len(turing_monostable_df)
print(df_concat)

limits_ba = (0,10)
limits_Vm = (10, 1000)
limits_km = (0.1, 250)
limits_mu = (0.001, 50)
limits_db = (0.001, 10)
ylimits= [limits_ba,limits_ba,limits_ba,limits_ba,limits_ba,limits_ba,limits_Vm,limits_Vm,limits_Vm,limits_Vm,limits_Vm,limits_Vm,limits_km,limits_km,limits_km,limits_km,limits_km,limits_km,limits_km,limits_mu,limits_mu,limits_db,limits_db,limits_db]
par_dict = df_concat.iloc[0].to_dict()
parameter_list = [key for key, value in par_dict.items()]
sns.set(style="white", palette="muted", color_codes=False)

# fig, axs = plt.subplots(nrows=4, ncols=1, figsize=(35,20))
fig, axs = plt.subplots(nrows=4, ncols=6, figsize=(24,16))
axs = axs.flatten()

palette = sns.diverging_palette(10, 220, n=2)
for count, parameter in enumerate(parameter_list[:23]):
    LogMin, LogMax = np.log10(df_concat.iloc[:,count].min()),np.log10(df_concat.iloc[:,count].max())
    newBins = np.logspace(LogMin, LogMax,100)
#     sns.histplot(df_concat.iloc[lenght_turing:,count].values, bins=newBins, kde=False, color = 'grey',alpha=0.001, ax = axs[count])
#     sns.histplot(df_concat.iloc[:lenght_turing,count].values, bins=newBins, kde=False, color = 'mediumseagreen', alpha=0.001, ax = axs[count])
#     lower_limit_distribution = np.amin(big_loguniform.iloc[:,count])
#     upper_limit_distribution = np.amax(big_loguniform.iloc[:,count])
    # sns.kdeplot(df_concat.iloc[:,count].values, fill=True,log_scale=True,cut=1,color='teal', linewidth = 4, alpha = 0.1, ax = axs[count])
    axs[count].hist(df_concat.iloc[:,count].values, bins = len(df_concat))
    print(df_concat.iloc[:,count].values)
    axs[count].set_xscale('log')
#     axs[count].set_yticklabels(axs[count].get_yticks(), size = 20)
    axs[count].set_ylabel('f', size = 0.001)
    # axs[count].set_xticklabels(axs[count].get_xticks(), size = 20)
    axs[count].set_title(str(parameter),fontsize=10)
    # axs[count].set(ylim=(0,1))
    axs[count].set(xlim=ylimits[count])




fig.tight_layout()
# fig.text(0.3,1.01,'Parameter Distribution of Turing I (Monostable)', size=50)
# fig.text(0.5,-0.03, "Parameter Value (log-scale)", ha="center", va="center",  fontsize=40)
# fig.text(-0.03,0.5, "Relative robustness", ha="center", va="center", rotation=90 ,fontsize=40)
# fig.savefig('interesting_parIDs_distributions.png', bbox_inches="tight")
plt.show()
