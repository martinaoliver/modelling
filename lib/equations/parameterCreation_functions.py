from unicodedata import name
from sklearn.utils import column_or_1d

import numpy as np
import pandas as pd
from tqdm import tqdm
import seaborn as sns
from matplotlib import pyplot as plt


np.random.seed(1)

def lhs(data, nsample,seed=1, tqdm_disable=False):
    np.random.seed(seed)
    m, nvar = data.shape
    ran = np.random.uniform(size=(nsample, nvar))
    s = np.zeros((nsample, nvar))
    for j in tqdm(range(0, nvar), disable=tqdm_disable):
        idx = np.random.permutation(nsample) + 1
        P = ((idx - ran[:, j]) / nsample) * 100
        s[:, j] = np.percentile(data[:, j], P)

    if np.any(s<=0):
        print('WARNING: negative values in lhs')
        s[s<0] = 0.001
        
    return s

def loguniform(size, low=-3, high=3):
    return (10) ** (np.random.uniform(low, high, size))



def parameterGaussian(name, mean, noisetosignal, size):
    stdev = noisetosignal * mean
    gaussianDistribution = np.random.normal(mean, stdev, size)
    return gaussianDistribution

def parameterLogNormal(name, mean, noisetosignal, size):
    sigma = noisetosignal * mean
    normal_std = np.sqrt(np.log(1 + (sigma/mean)**2))
    normal_mean = np.log(mean) - normal_std**2 / 2
    lognormalDistribution = np.random.lognormal(normal_mean, normal_std, size)
    return lognormalDistribution

def parameterLogUniform(name, min, max, size):
    loguniformDistribution = loguniform(size)
    croppedLoguniformDistribution = np.array([x for x in loguniformDistribution if min <= x <= max])
    return croppedLoguniformDistribution

def parameterFixed(name, value, size):
    fixedDistribution = np.full((size), value)
    return fixedDistribution



def parameterDistribution(parameterDict,size):
    if parameterDict['distribution']=='gaussian':
        dist = parameterGaussian(parameterDict['name'],parameterDict['mean'], parameterDict['noisetosignal'],size)
    if parameterDict['distribution']=='lognormal':
        dist = parameterLogNormal(parameterDict['name'],parameterDict['mean'], parameterDict['noisetosignal'],size)
    if parameterDict['distribution']=='loguniform':
        dist =  parameterLogUniform(parameterDict['name'],parameterDict['min'], parameterDict['max'],size)
    if parameterDict['distribution']=='fixed':
        dist =  parameterFixed(parameterDict['name'],parameterDict['value'],size)
    # sns.histplot(dist, bins=100)
    # plt.show()
    return dist

def preLhs(parameterDictList):
    parameterDistributionList = [parameterDistribution(parameterDict,100000) for parameterDict in parameterDictList]
    distributionMinimumLenght = np.amin([len(x) for x in parameterDistributionList])
    croppedParameterDistributionList = [x[:distributionMinimumLenght] for x in parameterDistributionList]
    # for dist in croppedParameterDistributionList:
    #     sns.histplot(dist, bins=100)
    #     plt.show()
    stackedDistributions = np.column_stack((croppedParameterDistributionList))
    return stackedDistributions

def plotDist(parameterDictList,lhsDist_df):
    nvar = len(parameterDictList)
  
    fig,axs = plt.subplots(nrows=1,ncols=nvar,figsize=(nvar*5,5))
    for count,parameter in enumerate(parameterDictList):
        name = parameter['name']
        lhsDistColumn = lhsDist_df[name]
        sns.histplot(lhsDistColumn, ax=axs[count], bins=100)
        axs[count].set(ylabel ='',yticks=[],yticklabels=[])
        axs[count].set_xlabel(name, fontsize=15)
        # axs[count].set_xscale('log')
    plt.show()

# def plotDist(parameterDictList,lhsDist_df):

#     nvar=len(parameterDictList)
#     n_col = int(np.sqrt(nvar))
#     n_col = 7
#     n_row = int(np.floor(nvar/n_col)+1)    # number of rows in the figure of the cluster


#     fig = plt.figure(figsize=(n_col/44, n_row/4))
#     for count,parameter in enumerate(parameterDictList):
#         ax=plt.subplot(n_row,n_col, count+1)
#         name = parameter['name']
#         lhsDistColumn = lhsDist_df[name]
#         sns.histplot(lhsDistColumn, ax=ax, bins=100)
#         ax.set(ylabel='',yticks=[],yticklabels=[])

#     plt.show()
#     plt.close
