#############
###paths#####
#############
import sys
import os



pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')
#############
from equations.parameterCreation_functions import *
#############
import numpy as np
import pandas as pd

# %matplotlib inline



# Gaussiansmall = {'name':'Gaussiansmall','distribution':'lognormal', 'mean':1, 'sigma':0.1}
# Gaussianbig = {'name':'Gaussianbig','distribution':'lognormal', 'mean':10, 'sigma':0.1}
# LogUniform = {'name':'LogUniform','distribution':'loguniform', 'min':1, 'max':100}
# Fixed = {'name':'Fixed','distribution':'fixed','value':1}

# nsamples=1000

# # parameterDictList = [muLVA,muASV,muU,muV]
# parameterDictList = [Gaussiansmall,Gaussianbig, LogUniform,Fixed]
# stackedDistributions = preLhs(parameterDictList)

# lhsDist = lhs(stackedDistributions,nsamples)

# lhsDist_df = pd.DataFrame(data = lhsDist, columns=[parameter['name'] for parameter in parameterDictList])

# plotDist(parameterDictList,lhsDist_df)

import seaborn as sns
import matplotlib.pyplot as plt

# fig, ax = plt.subplots(1)

# ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
vx=100
mean = np.log(vx)
noisetosignal=0.3
sigma = (noisetosignal**2)*(mean**2)
a = np.random.lognormal(mean=mean,sigma=sigma,size=1000)
sns.histplot(a)
plt.xscale('log')
# plt.xticklabels(tick_labels.astype(int))
plt.show()