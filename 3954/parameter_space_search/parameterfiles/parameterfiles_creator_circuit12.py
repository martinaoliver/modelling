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


muLVA = {'name':'muLVA','distribution':'gaussian', 'mean':2, 'stdev':3}
muASV = {'name':'muASV','distribution':'gaussian', 'mean':2, 'stdev':3}
muV = {'name':'muV','distribution':'loguniform', 'min':2, 'max':10}
muU = {'name':'muU','distribution':'loguniform', 'min':2, 'max':3}

nsamples=10

parameterDictList = [muLVA,muASV,muU,muV]
stackedDistributions = preLhs(parameterDictList)
lhsDist = lhs(stackedDistributions,nsamples)
print(lhsDist)