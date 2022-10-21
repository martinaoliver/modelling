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
circuit_n=12
variant=1
#diffusion parameters
DU = {'name':'DU','distribution':'gaussian', 'mean':1, 'noisetosignal':0.05}
DV= {'name':'DV','distribution':'gaussian', 'mean':1, 'noisetosignal':0.05}
D_parameters = [DU,DV]

#background parameters
bA = {'name':'bA','distribution':'gaussian', 'mean':0.1, 'noisetosignal':0.05}
bB = {'name':'bB','distribution':'gaussian', 'mean':0.1, 'noisetosignal':0.05}
bC = {'name':'bC','distribution':'gaussian', 'mean':0.1, 'noisetosignal':0.05}
bD = {'name':'bD','distribution':'gaussian', 'mean':0.1, 'noisetosignal':0.05}
bE = {'name':'bE','distribution':'gaussian', 'mean':0.1, 'noisetosignal':0.05}
bF = {'name':'bF','distribution':'gaussian', 'mean':0.1, 'noisetosignal':0.05}
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
Kbd = {'name':'Kbd','distribution':'gaussian', 'mean':430, 'noisetosignal':0.05} #lit 430, exp 870
Kab = {'name':'Kab','distribution':'gaussian', 'mean':0.4, 'noisetosignal':0.05} 
Kda = {'name':'Kda','distribution':'gaussian', 'mean':0.6, 'noisetosignal':0.05} 
Kfe = {'name':'Kfe','distribution':'gaussian', 'mean':20, 'noisetosignal':0.05}
Kee = {'name':'Kee','distribution':'gaussian', 'mean':20, 'noisetosignal':0.05}
Keb = {'name':'Keb','distribution':'gaussian', 'mean':20, 'noisetosignal':0.05}
Kce = {'name':'Kce','distribution':'gaussian', 'mean':5.6, 'noisetosignal':0.05}
KaTc = {'name':'KaTc','distribution':'gaussian', 'mean':1.26*10**3, 'noisetosignal':0.05}
Kiptg = {'name':'Kiptg','distribution':'gaussian', 'mean':1.5*10**3, 'noisetosignal':0.05}
K_parameters = [Kbd,Kab,Kda,Kfe,Kee,Keb,Kce,KaTc,Kiptg]

#degradation parameters (mu)
muLVA = {'name':'muLVA','distribution':'gaussian', 'mean':1.143, 'noisetosignal':0.05}
muAAV = {'name':'muAAV','distribution':'gaussian', 'mean':0.633, 'noisetosignal':0.05}
muASV = {'name':'muASV','distribution':'gaussian', 'mean':0.3, 'noisetosignal':0.05}
muUb = {'name':'muUb','distribution':'gaussian', 'mean':0.0225, 'noisetosignal':0.05}
muVb = {'name':'muVb','distribution':'gaussian', 'mean':0.0225, 'noisetosignal':0.05}
muaTc = {'name':'muaTc','distribution':'gaussian', 'mean':0.0101, 'noisetosignal':0.05}
muU = {'name':'muU','distribution':'loguniform', 'min':0.0225, 'max':3}
muV = {'name':'muV','distribution':'loguniform', 'min':0.0225, 'max':3}
mu_parameters = [muLVA,muAAV,muASV,muUb,muVb,muaTc,muU,muV]

#cooperativity parameters (n)
nbd = {'name':'nbd','distribution':'fixed', 'value':2}
nab = {'name':'nab','distribution':'fixed', 'value':1}
nda = {'name':'nda','distribution':'fixed', 'value':1}
nfe = {'name':'nfe','distribution':'fixed', 'value':4}
nee = {'name':'nee','distribution':'fixed', 'value':4}
neb = {'name':'neb','distribution':'fixed', 'value':4}
nce = {'name':'nce','distribution':'fixed', 'value':1}
naTc = {'name':'naTc','distribution':'fixed', 'value':2}
niptg = {'name':'niptg','distribution':'fixed', 'value':1}
n_parameters = [nbd,nab,nda,nfe,nee,neb,nce,naTc,niptg]

#kinetic rates parameters (k)
k1 = {'name':'k1','distribution':'gaussian', 'mean':0.0183, 'noisetosignal':0.05}
k2 = {'name':'k2','distribution':'gaussian', 'mean':0.0183, 'noisetosignal':0.05}
k_parameters = [k1,k2]

#molecule concentrations parameters (C)
iptg= {'name':'iptg','distribution':'gaussian', 'mean':1000, 'noisetosignal':0.05}
#atc concentration comes as an initial condition


nsamples=int(sys.argv[1])
plotDistributions=False
if plotDistributions == True:
    parameterTypeList = [D_parameters,b_parameters,V_parameters,K_parameters,mu_parameters,n_parameters,k_parameters]
    for parameterType in parameterTypeList:
        stackedDistributions = preLhs(parameterType)
        lhsDist = lhs(stackedDistributions,nsamples)
        lhsDist_df = pd.DataFrame(data = lhsDist, columns=[parameter['name'] for parameter in parameterType])
        plotDist(parameterType,lhsDist_df)

parameterDictList = [DU, DV, bA, bB, bC, bD, bE, bF, VA, VB, VC, VD, VE, VF, Kbd, Kab, Kda, Kfe, Kee, Keb, Kce, KaTc, Kiptg, muLVA, muAAV, muASV, muUb, muVb, muaTc, muU, muV, nbd, nab, nda, nfe, nee, neb, nce, naTc, niptg, k1, k2, iptg]
stackedDistributions = preLhs(parameterDictList)
lhsDist = lhs(stackedDistributions,nsamples)
lhsDist_df = pd.DataFrame(data = lhsDist, columns=[parameter['name'] for parameter in parameterDictList])
# plotDist(parameterDictList,lhsDist_df)
pkl.dump(lhsDist_df, open(modellingpath + '/3954/paper/input/parameterfiles/df_circuit%r_variant%r_%rparametersets.pkl'%(circuit_n,variant,nsamples), 'wb'))

print(lhsDist_df)