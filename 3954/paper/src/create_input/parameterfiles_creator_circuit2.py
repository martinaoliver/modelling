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
import pickle as pkl
# %matplotlib inline
circuit_n=2
variant=0
#diffusion parameters
# DU = {'name':'DU','distribution':'gaussian', 'mean':1, 'noisetosignal':0.05}
# DV= {'name':'DV','distribution':'gaussian', 'mean':1, 'noisetosignal':0.05}
DA = {'name':'DA','distribution':'fixed', 'value':1}
DB = {'name':'DB','distribution':'loguniform', 'min':0.001, 'max':10}

D_parameters = [DA,DB]

#background parameters
bA = {'name':'bA','distribution':'fixed', 'value':0.01}
bB = {'name':'bB','distribution':'fixed', 'value':0.01}
bC = {'name':'bC','distribution':'fixed', 'value':0.01}
bD = {'name':'bD','distribution':'fixed', 'value':0.01}
bE = {'name':'bE','distribution':'fixed', 'value':0.01}
bF = {'name':'bF','distribution':'fixed', 'value':0.01}

b_parameters = [bA,bB,bC,bD,bE,bF]

#maximum production parameters (V)


VA = {'name':'VA','distribution':'loguniform', 'min':10, 'max':1000}
VB = {'name':'VB','distribution':'loguniform', 'min':10, 'max':1000}
VC = {'name':'VC','distribution':'loguniform', 'min':10, 'max':1000}
VD = {'name':'VD','distribution':'loguniform', 'min':10, 'max':1000}
VE = {'name':'VE','distribution':'loguniform', 'min':10, 'max':1000}
VF = {'name':'VF','distribution':'loguniform', 'min':10, 'max':1000}
V_parameters = [VA,VB,VC,VD,VE,VF]

#[] at half activation parameters (K)
Kda = {'name':'Kda','distribution':'loguniform', 'min':10, 'max':1000}
Kab = {'name':'Kab','distribution':'loguniform', 'min':10, 'max':1000}
Keb = {'name':'Keb','distribution':'loguniform', 'min':10, 'max':1000}
Kbd = {'name':'Kbd','distribution':'loguniform', 'min':10, 'max':1000}
Kfe = {'name':'Kfe','distribution':'loguniform', 'min':10, 'max':1000}
Kee = {'name':'Kee','distribution':'loguniform', 'min':10, 'max':1000}
Kce = {'name':'Kce','distribution':'loguniform', 'min':10, 'max':1000}
K_parameters = [Kda, Kab, Keb, Kbd, Kfe, Kee, Kce]

#degradation parameters (mu)
muLVA = {'name':'muLVA','distribution':'loguniform', 'min':0.001, 'max':50}
muASV = {'name':'muASV','distribution':'loguniform', 'min':0.001, 'max':50}
mu_parameters = [muLVA,muASV]

#cooperativity parameters (n)
nbd = {'name':'nbd','distribution':'fixed', 'value':2}
nab = {'name':'nab','distribution':'fixed', 'value':1}
nda = {'name':'nda','distribution':'fixed', 'value':1}
nfe = {'name':'nfe','distribution':'fixed', 'value':4}
nee = {'name':'nee','distribution':'fixed', 'value':4}
neb = {'name':'neb','distribution':'fixed', 'value':4}
nce = {'name':'nce','distribution':'fixed', 'value':1}
n_parameters = [nbd,nab,nda,nfe,nee,neb,nce]



nsamples=int(sys.argv[1])
plotDistributions=False
if plotDistributions == True:
    parameterTypeList = [D_parameters,b_parameters,V_parameters,K_parameters,mu_parameters,n_parameters,k_parameters]
    for parameterType in parameterTypeList:
        stackedDistributions = preLhs(parameterType)
        lhsDist = lhs(stackedDistributions,nsamples)
        lhsDist_df = pd.DataFrame(data = lhsDist, columns=[parameter['name'] for parameter in parameterType])
        plotDist(parameterType,lhsDist_df)

parameterDictList = D_parameters + b_parameters + V_parameters + K_parameters + mu_parameters + n_parameters
# parameterDictList = [DU, DV, bA, bB, bC, bD, bE, bF, VA, VB, VC, VD, VE, VF, Kbd, Kab, Kda, Kfe, Kee, Keb, Kce, KaTc, Kiptg, muLVA, muAAV, muASV, muUb, muVb, muaTc, muU, muV, nbd, nab, nda, nfe, nee, neb, nce, naTc, niptg, k1, k2, iptg]
stackedDistributions = preLhs(parameterDictList)
lhsDist = lhs(stackedDistributions,nsamples)
lhsDist_df = pd.DataFrame(data = lhsDist, columns=[parameter['name'] for parameter in parameterDictList])
# plotDist(parameterDictList,lhsDist_df)
pkl.dump(lhsDist_df, open(modellingpath + '/3954/paper/input/lhs_parameterfiles/df_circuit%r_variant%r_%rparametersets.pkl'%(circuit_n,variant,nsamples), 'wb'))

print(lhsDist_df)