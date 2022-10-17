from sklearn.utils import column_or_1d

import numpy as np
import pandas as pd
from tqdm import tqdm

def lhs(data, nsample):
    nvar,m = data.shape
    ran = np.random.uniform(size=(nsample, nvar))
    s = np.zeros((nsample, nvar))
    for j in tqdm(range(0, nvar)):
        idx = np.random.permutation(nsample) + 1
        P = ((idx - ran[:, j]) / nsample) * 100
        s[:, j] = np.percentile(data[:, j], P)

    return s





def loguniform(size, low=-3, high=3):
    return (10) ** (np.random.uniform(low, high, size))



def parameterGaussian(name, mean, stdev, size):
    gaussianDistribution = np.random.normal(mean, stdev, size)
    # pd_column = pd.DataFrame({name:distribution})
    return gaussianDistribution

def parameterLoguniform(name, min, max, size):
    loguniformDistribution = loguniform(size)
    croppedLoguniformDistribution = np.array([x for x in loguniformDistribution if min <= x <= max])
    return croppedLoguniformDistribution

def parameterFixed(name, value, size):
    fixedDistribution = np.full((size), value)
    return fixedDistribution



def parameterDistribution(parameterDict,size):
    if parameterDict['distribution']=='gaussian':
        dist = parameterGaussian(parameterDict['name'],parameterDict['mean'], parameterDict['stdev'],size)
    if parameterDict['distribution']=='loguniform':
        dist =  parameterLoguniform(parameterDict['name'],parameterDict['min'], parameterDict['max'],size)
    if parameterDict['distribution']=='fixed':
        dist =  parameterLoguniform(parameterDict['name'],parameterDict['value'],size)
    return dist

def preLhs(parameterDictList):
    parameterDistributionList = [parameterDistribution(parameterDict,100000) for parameterDict in parameterDictList]
    distributionMinimumLenght = np.amin([len(x) for x in parameterDistributionList])
    croppedParameterDistributionList = [x[:distributionMinimumLenght] for x in parameterDistributionList]
    stackedDistributions = np.stack(croppedParameterDistributionList,axis=0)
    return stackedDistributions
