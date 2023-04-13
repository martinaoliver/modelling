#############
###paths#####
#############
import sys
import os

pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')
#############
import pickle as pkl
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error
from scipy.optimize import curve_fit

from equations.parameterCreation_functions import *

#plotting functions
def plotData(inducer, rfpExp_list, gfpExp_list, semRed, semGreen,pad=0.01, inducerName='OC14'):
    fig,ax = plt.subplots()

    ax.plot(inducer,rfpExp_list,label='RFP', c='red')
    ax.scatter(inducer,rfpExp_list, c='red')
    ax.errorbar(inducer,rfpExp_list,yerr=semRed,c='red',fmt='o')
    ax.legend(loc='center left') #upper right
    ax.set_ylabel('RFP / ($A_{600}$ $RFP_{basal})$')
    ax.set_xscale('log')


    ax2=ax.twinx()
    ax2.plot(inducer,gfpExp_list,label='GFP', c='green')
    ax2.scatter(inducer,gfpExp_list, c='green')
    ax2.errorbar(inducer,gfpExp_list,yerr=semGreen,c='green',fmt='o')
    ax2.legend(loc='center right') #upper left
    ax2.set_ylabel('GFP / ($A_{600}$ $GFP_{basal})$')
    ax.set_xscale('log')

    ax.set_xlabel(f'{inducerName} concentration (µM)')

    plt.show()

def plotFitvsData(inducer,inducer_continuous, gfpExp_list, rfpExp_list, semGreen, semRed,doseResponseGreen,doseResponseRed,pad=0.01):
    fig,ax = plt.subplots()

    ax.plot(inducer_continuous,doseResponseRed,label='RFP', c='red')
    ax.scatter(inducer,rfpExp_list, c='red')
    ax.errorbar(inducer,rfpExp_list,yerr=semRed,c='red',fmt='o')
    ax.legend(loc='center left') #upper right
    ax.set_ylabel('RFP / ($A_{600}$ $RFP_{basal})$')
    ax.set_xscale('log')


    ax2=ax.twinx()
    ax2.plot(inducer_continuous,doseResponseGreen,label='GFP', c='green')
    ax2.scatter(inducer,gfpExp_list, c='green')
    ax2.errorbar(inducer,gfpExp_list,yerr=semGreen,c='green',fmt='o')
    ax2.legend(loc='center right') #upper left
    ax2.set_ylabel('GFP / ($A_{600}$ $GFP_{basal})$')
    ax.set_xscale('log')


    plt.show()

#load data
doseResponseExp1 = pkl.load(open('input/liquid_culture/curatedData/nonTransferred_subcircuit1_circuit14_doseResponseOC14_-0.5ATC.pkl','rb'))
OC14_list1= np.array(doseResponseExp1['OC14']); gfpExp_list1 = list(doseResponseExp1['mean_gfp']); rfpExp_list1 = list(doseResponseExp1['mean_rfp'])
semGreen1 = doseResponseExp1['std_gfp']; semRed1 = doseResponseExp1['std_rfp']
plotData(OC14_list1, rfpExp_list1, gfpExp_list1, semRed1, semGreen1)


doseResponseExp3 = pkl.load(open('input/liquid_culture/curatedData/nonTransferred_subcircuit3_circuit14_doseResponseOC14_-0.5ATC.pkl','rb'))
OC14_list3= np.array(doseResponseExp3['OC14']); gfpExp_list3 = list(doseResponseExp3['mean_gfp']); rfpExp_list3 = list(doseResponseExp3['mean_rfp'])
semGreen3 = doseResponseExp3['std_gfp']; semRed3 = doseResponseExp3['std_rfp']
plotData(OC14_list3, rfpExp_list3, gfpExp_list3, semRed3, semGreen3)

##fitting params

nvd = 2
nfe = 5
nda=2
nce=3

def gfp1_steadystate(OC14, Vf,Kvd): 
    muv = 0.0225 ; kv =  0.0183 ;
    F1 = 1 + Vf*(1/(1+((muv*Kvd)/(kv*OC14 + 1e-8))**nvd ))
    return F1

def rfp1_steadystate(OC14, Vf,Kvd,VeKfenfe): 
    muv = 0.0225 ; kv =  0.0183 ;
    E1 = 1 + VeKfenfe*(1/(1+((1 + Vf*(1/(1+((muv*Kvd)/(kv*OC14 + 1e-8))**nvd ))))**nfe))

    return E1


def gfp3_steadystate(OC14,  Vd,Kvd): 
    muv = 0.0225 ; kv =  0.0183 ;
    D3 = 1 + Vd*(1/(1+((muv*Kvd)/(kv*OC14 + 1e-8))**nvd ))
    return D3

def bfp3_steadystate(D,VcKdanda): 
    C3= 1 + VcKdanda*(1/(1+((D)**nda)))
    return C3

def rfp3_steadystate(OC14,  Vd,Kvd,VcKdanda,  VeKcence): 
    E3 = 1 + VeKcence*(1/(1+((bfp3_steadystate(gfp3_steadystate(OC14,  Vd,Kvd), VcKdanda))**nce)))
    return E3


OC14_continuous = np.logspace(-3,1, 100)


def steadystate(OC14,Vf,Vd, Kvd, VeKfenfe, VcKdanda, VeKcence):
  
  if len(OC14) == 22:
      gaps = [5,5,6,6]
  else:
        gaps = [int(len(OC14)/4)]*4
  F1 = gfp1_steadystate(OC14[:np.sum(gaps[:1])],  Vf,Kvd)
  E1 = rfp1_steadystate(OC14[np.sum(gaps[:1]):np.sum(gaps[:2])],Vf,Kvd,VeKfenfe)
  D3 = gfp3_steadystate(OC14[np.sum(gaps[:2]):np.sum(gaps[:3])], Vd,Kvd)
  E3 = rfp3_steadystate(OC14[np.sum(gaps[:3]):np.sum(gaps[:4])],  Vd,Kvd,VcKdanda,  VeKcence)
  FE = np.hstack([F1,E1, D3, E3])
  return FE


fluorescenceData = np.hstack([gfpExp_list1,rfpExp_list1, gfpExp_list3,rfpExp_list3])
OC14data_new = np.hstack([OC14_list1,OC14_list1, OC14_list3,OC14_list3])
OC14data_continuous= np.hstack([OC14_continuous,OC14_continuous, OC14_continuous,OC14_continuous])
semStacked= np.hstack([semGreen1,semRed1, semGreen3,semRed3])

popt, pcov = curve_fit(f=steadystate, xdata=OC14data_new, ydata=fluorescenceData ,sigma =semStacked, bounds = (0,100), maxfev = 100000000)

paramNames = ['Vf','Vd', 'Kvd', 'VeKfenfe', 'VcKdanda', 'VeKcence']
pfitDict = {}
for param in popt:
    pfitDict[paramNames[popt.tolist().index(param)]] = param

fluorescenceFit = steadystate(OC14data_new, *popt)
fluorescenceFit_continuous = steadystate(OC14data_continuous, *popt)
gfpFit1 = fluorescenceFit[:5]; rfpFit1 = fluorescenceFit[5:10]; gfpFit3 = fluorescenceFit[10:16]; rfpFit3 = fluorescenceFit[16:22]
gfpFit1_continuous = fluorescenceFit_continuous[:100]; rfpFit1_continuous = fluorescenceFit_continuous[100:200]; gfpFit3_continuous = fluorescenceFit_continuous[200:300]; rfpFit3_continuous = fluorescenceFit_continuous[300:400]

plotFitvsData(OC14_list1,OC14_continuous, gfpExp_list1, rfpExp_list1, semGreen1, semRed1, gfpFit1_continuous,rfpFit1_continuous)

plotFitvsData(OC14_list3,OC14_continuous, gfpExp_list3, rfpExp_list3, semGreen3, semRed3, gfpFit3_continuous,rfpFit3_continuous)
gfpFit1_continuous_copy,rfpFit1_continuous_copy, gfpFit3_continuous_copy,rfpFit3_continuous_copy = gfpFit1_continuous,rfpFit1_continuous, gfpFit3_continuous,rfpFit3_continuous 


#creating df
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
variant='fitted3'
#diffusion parameters


#maximum production parameters (V*)
# V* = V/b
# V = 10-1000
# b=0.1-1
minV = 10;maxV=1000;minb=0.01;maxb=1
Va = {'name':'Va','distribution':'loguniform', 'min':minV/maxb, 'max':maxV/minb}
Vb = {'name':'Vb','distribution':'loguniform', 'min':minV/maxb, 'max':maxV/minb}

# Vc = {'name':'Vc','distribution':'gaussian', 'mean':pfitDict['Vc'], 'noisetosignal':0.3}
# Vd = {'name':'Vd','distribution':'gaussian', 'mean':pfitDict['Vd'], 'noisetosignal':0.3}
# Ve = {'name':'Ve','distribution':'gaussian', 'mean':pfitDict['Ve'], 'noisetosignal':0.3}
# Vf = {'name':'Vf','distribution':'gaussian', 'mean':pfitDict['Vf'], 'noisetosignal':0.3}


Vc = {'name':'Vc','distribution':'lognormal', 'mean':pfitDict['Vc'], 'noisetosignal':1}
Vd = {'name':'Vd','distribution':'lognormal', 'mean':pfitDict['Vd'], 'noisetosignal':1}
Ve = {'name':'Ve','distribution':'lognormal', 'mean':pfitDict['Ve'], 'noisetosignal':1}
Vf = {'name':'Vf','distribution':'lognormal', 'mean':pfitDict['Vf'], 'noisetosignal':1}

V_parameters = [Va,Vb, Vc, Vd, Ve, Vf]



K1=0.0183; K2=0.0183
DUmin=1; DUmax=1; DVmin=0.001; DVmax=1
muU=0.0225; muV=0.0225
KdiffpromMin=0.1;KdiffpromMax=250
muLVA_estimate =1.143
muAAV_estimate =0.633
muASV_estimate=0.300 #this corresponds to mua

def Dr_(K1,K2,muU,muV,DU_,DV_):
    return (K2*DV_*muU)/(K1*DU_*muV)


def Kdiffstar(mudiff,Kdiffprom,kdiff):
    return mudiff*Kdiffprom/kdiff

def Kstar(mu,b,K):
    return mu/b*K

Dr = {'name':'Dr','distribution':'loguniform', 'min':Dr_(K1,K2,muU,muV,DUmax,DVmin), 'max':Dr_(K1,K2,muU,muV,DUmin,DVmax)}
D_parameters = [Dr]

#[] at half activation parameters (K)
minK=0.1;maxK=250

Keb = {'name':'Keb','distribution':'loguniform', 'min':Kstar(muLVA_estimate,maxb,minK), 'max':Kstar(muLVA_estimate,minb,maxK)}

Kee = {'name':'Kee','distribution':'fixed','value':0.001}

Kub = {'name':'Kub','distribution':'loguniform', 'min':Kdiffstar(muU,KdiffpromMin,K1), 'max':Kdiffstar(muU,KdiffpromMax,K1)}

# Kvd = {'name':'Kvd','distribution':'gaussian', 'mean':pfitDict['Kvd'], 'noisetosignal':0.3}
# Kda = {'name':'Kda','distribution':'gaussian', 'mean':pfitDict['Kda'], 'noisetosignal':0.3}
# Kce = {'name':'Kce','distribution':'gaussian', 'mean':pfitDict['Kce'], 'noisetosignal':0.3}
# Kfe = {'name':'Kfe','distribution':'gaussian', 'mean':pfitDict['Kfe'], 'noisetosignal':0.3}


Kvd = {'name':'Kvd','distribution':'lognormal', 'mean':pfitDict['Kvd'], 'noisetosignal':1}
Kda = {'name':'Kda','distribution':'lognormal', 'mean':pfitDict['Kda'], 'noisetosignal':1}
Kce = {'name':'Kce','distribution':'lognormal', 'mean':pfitDict['Kce'], 'noisetosignal':1}
Kfe = {'name':'Kfe','distribution':'lognormal', 'mean':pfitDict['Kfe'], 'noisetosignal':1}
print(pfitDict['Kvd'], pfitDict['Kda'], pfitDict['Kce'], pfitDict['Kfe'])

# Kab = {'name':'Kab','distribution':'loguniform', 'min':muU_low*DU_low/k1, 'max':muU_high*DU_high/k1}
# Kbd = {'name':'Kbd','distribution':'loguniform', 'min':muV_low*DV_low/k2, 'max':muV_high*DV_high/k2}
K_parameters = [ Kub, Keb, Kee, Kvd, Kda, Kce, Kfe]



#protein degradation parameters (mu)
# mu = mux/mua
muASV = {'name':'muASV','distribution':'fixed', 'value':muASV_estimate/muASV_estimate}
muLVA = {'name':'muLVA','distribution': 'lognormal','mean':muLVA_estimate /muASV_estimate, 'noisetosignal':1}
mu_parameters = [muLVA,muASV]


#cooperativity parameters (n)
nvd = {'name':'nvd','distribution':'fixed', 'value':2}
nfe = {'name':'nfe','distribution':'fixed', 'value':5} #ideally hihger but we keep lower to ensure numerics work
nda = {'name':'nda','distribution':'fixed', 'value':2}
nce = {'name':'nce','distribution':'fixed', 'value':3}
nub = {'name':'nub','distribution':'fixed', 'value':1}
nee = {'name':'nee','distribution':'fixed', 'value':4}
neb = {'name':'neb','distribution':'fixed', 'value':4}
nfe = {'name':'nfe','distribution':'fixed', 'value':8}
n_parameters = [nub,nee,neb,nvd,nda,nce,nfe]



plotDistributions=False
if plotDistributions == True:
    D_parameters = [Dr, Dr]
    nsamples=1000
    parameterTypeList = [ D_parameters  , V_parameters , K_parameters , mu_parameters , n_parameters]

    for parameterType in parameterTypeList:
        stackedDistributions = preLhs(parameterType)
        lhsDist = lhs(stackedDistributions,nsamples)
        lhsDist_df = pd.DataFrame(data = lhsDist, columns=[parameter['name'] for parameter in parameterType])
        plotDist(parameterType,lhsDist_df)

createParams=True
if createParams == True:
    nsamples=2000000
    # nsamples=int(sys.argv[1])
    # nsamples=14
    parameterDictList = D_parameters  + V_parameters + K_parameters + mu_parameters + n_parameters
    # parameterDictList = [DU, DV, bA, bB, bC, bD, bE, bF, VA, VB, VC, VD, VE, VF, Kbd, Kab, Kda, Kfe, Kee, Keb, Kce, KaTc, Kiptg, muLVA, muAAV, muASV, muUb, muVb, muaTc, muU, muV, nbd, nab, nda, nfe, nee, neb, nce, naTc, niptg, k1, k2, iptg]
    stackedDistributions = preLhs(parameterDictList)
    lhsDist = lhs(stackedDistributions,nsamples)
    lhsDist_df = pd.DataFrame(data = lhsDist, columns=[parameter['name'] for parameter in parameterDictList])
    # plotDist(parameterDictList,lhsDist_df)
    pkl.dump(lhsDist_df, open(modellingpath + '/3954/paper/input/fitted_parameterfiles/df_circuit%r_variant%s_%rparametersets.pkl'%(circuit_n,variant,nsamples), 'wb'))

    print(lhsDist_df)
