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
circuit_n=14
variant='1nd'
#diffusion parameters


#maximum production parameters (V*)
# V* = V/b
# V = 10-1000
# b=0.1-1
minV = 10;maxV=1000;minb=0.1;maxb=1
VA = {'name':'VA','distribution':'loguniform', 'min':minV/maxb, 'max':maxV/minb}
VB = {'name':'VB','distribution':'loguniform', 'min':minV/maxb, 'max':maxV/minb}
VC = {'name':'VC','distribution':'loguniform', 'min':minV/maxb, 'max':maxV/minb}
VD = {'name':'VD','distribution':'loguniform', 'min':minV/maxb, 'max':maxV/minb}
VE = {'name':'VE','distribution':'loguniform', 'min':minV/maxb, 'max':maxV/minb}
VF = {'name':'VF','distribution':'loguniform', 'min':minV/maxb, 'max':maxV/minb}
V_parameters = [VA,VB,VC,VD,VE,VF]



K1=0.0183; K2=0.0183
DUmin=0.1; DUmax=10; DVmin=0.1; DVmax=10
muU=0.0225; muV=0.0225
KdiffpromMin=0.1;KdiffpromMax=250

def Dr_(K1,K2,muU,muV,DU_,DV_):
    return (K2*DV_*muU)/(K1*DU_*muV)


def Kdiffstar(mudiff,Kdiffprom,kdiff):
    return mudiff*Kdiffprom/kdiff





Dr = {'name':'Dr','distribution':'loguniform', 'min':Dr_(K1,K2,muU,muV,DUmax,DVmin), 'max':Dr_(K1,K2,muU,muV,DUmin,DVmax)}
D_parameters = [Dr, Dr]

#[] at half activation parameters (K)
Kda = {'name':'Kda','distribution':'loguniform', 'min':0.1, 'max':250}
Keb = {'name':'Keb','distribution':'loguniform', 'min':0.1, 'max':250}
Kfe = {'name':'Kfe','distribution':'loguniform', 'min':0.1, 'max':250}
Kee = {'name':'Kee','distribution':'loguniform', 'min':0.1, 'max':250}
Kce = {'name':'Kce','distribution':'loguniform', 'min':0.1, 'max':250}
Kbd = {'name':'Kbd','distribution':'loguniform', 'min':Kdiffstar(muV,KdiffpromMin,K1), 'max':Kdiffstar(muV,KdiffpromMax,K1)}
Kab = {'name':'Kab','distribution':'loguniform', 'min':Kdiffstar(muU,KdiffpromMin,K2), 'max':Kdiffstar(muU,KdiffpromMax,K2)}

# Kab = {'name':'Kab','distribution':'loguniform', 'min':muU_low*DU_low/k1, 'max':muU_high*DU_high/k1}
# Kbd = {'name':'Kbd','distribution':'loguniform', 'min':muV_low*DV_low/k2, 'max':muV_high*DV_high/k2}
K_parameters = [Kda, Kab, Keb, Kbd, Kfe, Kee, Kce]

#protein degradation parameters (mu)
# mu = mux/mua
# minmua=0.001;maxmua=50;minmux=0.001;maxmux=50
muLVA_estimate =1.143
muAAV_estimate =0.633
muASV_estimate=0.300 #this corresponds to mua

muLVA = {'name':'muLVA','distribution': 'gaussian','mean':muLVA_estimate/muASV_estimate, 'noisetosignal':0.1}
muASV = {'name':'muASV','distribution':'fixed', 'value':muASV_estimate/muASV_estimate}
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



plotDistributions=False
if plotDistributions == True:
    nsamples=10000
    parameterTypeList = [ D_parameters  , V_parameters , K_parameters , mu_parameters , n_parameters]

    for parameterType in parameterTypeList:
        stackedDistributions = preLhs(parameterType)
        lhsDist = lhs(stackedDistributions,nsamples)
        lhsDist_df = pd.DataFrame(data = lhsDist, columns=[parameter['name'] for parameter in parameterType])
        plotDist(parameterType,lhsDist_df)
createParams=True
if createParams == True:
    # nsamples=1000000
    nsamples=int(sys.argv[1])
    # nsamples=1000000
    parameterDictList = D_parameters  + V_parameters + K_parameters + mu_parameters + n_parameters
    # parameterDictList = [DU, DV, bA, bB, bC, bD, bE, bF, VA, VB, VC, VD, VE, VF, Kbd, Kab, Kda, Kfe, Kee, Keb, Kce, KaTc, Kiptg, muLVA, muAAV, muASV, muUb, muVb, muaTc, muU, muV, nbd, nab, nda, nfe, nee, neb, nce, naTc, niptg, k1, k2, iptg]
    stackedDistributions = preLhs(parameterDictList)
    lhsDist = lhs(stackedDistributions,nsamples)
    lhsDist_df = pd.DataFrame(data = lhsDist, columns=[parameter['name'] for parameter in parameterDictList])
    # plotDist(parameterDictList,lhsDist_df)
    pkl.dump(lhsDist_df, open(modellingpath + '/3954/paper/input/lhs_parameterfiles/df_circuit%r_variant%s_%rparametersets.pkl'%(circuit_n,variant,nsamples), 'wb'))

    print(lhsDist_df)

