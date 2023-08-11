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
import pickle as pkl
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error
from scipy.optimize import curve_fit
from tqdm import tqdm
import pandas as pd
import seaborn as sns

from equations.parameterCreation_functions import *
from fitting.load_dose_response import *
from fitting.dose_response_plotting import *
#%%


test=True


circuit_n=14
variant='fitted7'
n_samples = 13700000
# n_samples = 10
if test==True:
    n_samples=1000
#############
###load dataset#####
#############


HSLtransform = 0.14*10**3
filename = 'Jure_subcircuit1_circuit14_doseResponseOC14_0.5ATC.pkl'
OC14_list1_red, gfpExp_list1, rfpExp_list1, semGreen1, semRed1 = load_dose_response(filename)
print(OC14_list1_red)
del gfpExp_list1[1:6]
del semGreen1[1:6]
OC14_list1_green = [OC14_list1_red[x] for x in [0,6,7,8,9]]
plotData_croppedGreen(OC14_list1_green, OC14_list1_red, rfpExp_list1, gfpExp_list1, semRed1, semGreen1)

filename = 'Jure_subcircuit3_circuit14_doseResponseOC14_0.5ATC.pkl'
OC14_list3_red, gfpExp_list3, rfpExp_list3, semGreen3, semRed3 = load_dose_response(filename)

del gfpExp_list3[1:6]
del semGreen3[1:6]
OC14_list3_green = [OC14_list3_red[x] for x in [0,6,7,8,9]]
plotData_croppedGreen(OC14_list3_green, OC14_list3_red, rfpExp_list3, gfpExp_list3, semRed3, semGreen3)





HSLtransform = 0.14*10**3
filename = 'Jure_subcircuit1_circuit14_doseResponseOC14_0.5ATC.pkl'
OC14_list1_red, gfpExp_list1, rfpExp_list1, semGreen1, semRed1 = load_dose_response(filename)
print(OC14_list1_red)
del gfpExp_list1[1:6]
del semGreen1[1:6]
OC14_list1_green = [OC14_list1_red[x] for x in [0,6,7,8,9]]
plotData_croppedGreen(OC14_list1_green, OC14_list1_red, rfpExp_list1, gfpExp_list1, semRed1, semGreen1)

filename = 'Jure_subcircuit3_circuit14_doseResponseOC14_0.5ATC.pkl'
OC14_list3_red, gfpExp_list3, rfpExp_list3, semGreen3, semRed3 = load_dose_response(filename)

del gfpExp_list3[1:6]
del semGreen3[1:6]
OC14_list3_green = [OC14_list3_red[x] for x in [0,6,7,8,9]]
plotData_croppedGreen(OC14_list3_green, OC14_list3_red, rfpExp_list3, gfpExp_list3, semRed3, semGreen3)




#############
###fit dataset#####
#############

nvd = 2
nfe = 5
nda=2
nce=3

OC14_continuous = np.logspace(-3,2, 100)*HSLtransform


###subcircuit1
def gfp1_steadystate(OC14, Vf,Kvd): 
    muv = 0.0225 ; kv =  0.0183 ;
    F1 = 1 + Vf*(1/(1+((muv*Kvd)/(kv*OC14 + 1e-8))**nvd ))
    return F1


def rfp1_steadystate(OC14, Ve, Vf,Kvd,Kfe): 
    muv = 0.0225 ; kv =  0.0183 ;
    E1 = 1 + Ve*(1/(1+((1 + Vf*(1/(1+((muv*Kvd)/(kv*OC14 + 1e-8))**nvd )))/(Kfe+1e-8))**nfe))
    
    return E1


def steadystate_subcircuit1(OC14, Ve, Vf, Kvd, Kfe):
  
  if len(OC14) == 15:
      gaps = [5,10]
  else:
     gaps = [int(len(OC14)/2)]*2
  F1 = gfp1_steadystate(OC14[:np.sum(gaps[:1])],  Vf,Kvd)
  E1 = rfp1_steadystate(OC14[np.sum(gaps[:1]):np.sum(gaps[:2])],Ve, Vf,Kvd,Kfe)
  FE = np.hstack([F1,E1])
  return FE


fluorescenceData = np.hstack([gfpExp_list1,rfpExp_list1])
OC14data_new = np.hstack([OC14_list1_green,OC14_list1_red])
OC14data_continuous= np.hstack([OC14_continuous,OC14_continuous])
semStacked= np.hstack([semGreen1,semRed1])

popt_subcircuit1, pcov_subcircuit1 = curve_fit(f=steadystate_subcircuit1, xdata=OC14data_new, ydata=fluorescenceData,  p0 = [2,3.5,0.1,3],maxfev=10000000, bounds=(0,100))

paramNames_subcircuit1 = ['Ve1','Vf', 'Kvd','Kfe']
pfitDict = {}
for param in popt_subcircuit1:
    pfitDict[paramNames_subcircuit1[popt_subcircuit1.tolist().index(param)]] = param


###subcircuit3
Kvd = pfitDict['Kvd']
Ve = np.amax(rfpExp_list3) - 1

def gfp3_steadystate(OC14,  Vd): 
    muv = 0.0225 ; kv =  0.0183 ;
    D3 = 1 + Vd*(1/(1+((muv*Kvd)/(kv*OC14 + 1e-8))**nvd ))

    return D3

def bfp3_steadystate(D,Vc,Kda):
    C3 = 1 + Vc*(1/(1+((D/(Kda+1e-8))**nda)))
    return C3

def rfp3_steadystate(D,Vc, Kda,Kce): 

    E3 = 1 + Ve*(1/(1+((bfp3_steadystate(D, Vc, Kda)/(Kce+1e-8))**nce)))
    return E3


def steadystate_subcircuit3(OC14, Vc, Vd, Kda, Kce):
  if len(OC14) == 15:
      gaps = [5,10]
  else:
     gaps = [int(len(OC14)/2)]*2
  D3 = gfp3_steadystate(OC14[:np.sum(gaps[:1])], Vd)
  E3 = rfp3_steadystate(OC14[np.sum(gaps[:1]):np.sum(gaps[:2])], Vc, Kda,Kce)
  DE = np.hstack([D3,E3])
  return DE


fluorescenceData = np.hstack([gfpExp_list3,rfpExp_list3])
OC14data_new = np.hstack([OC14_list3_green,OC14_list3_red])
OC14data_continuous= np.hstack([OC14_continuous,OC14_continuous])
semStacked= np.hstack([semGreen3,semRed3])

popt_subcircuit3, pcov_subcircuit3 = curve_fit(f=steadystate_subcircuit3, xdata=OC14data_new, ydata=fluorescenceData, maxfev=10000, bounds=(0,100))

paramNames_subcircuit3 = ['Vc', 'Vd',  'Kda', 'Kce']
for count, param in enumerate(popt_subcircuit3):
    pfitDict[paramNames_subcircuit3[count]] = param

if test==True:
    fluorescenceFit = steadystate_subcircuit1(OC14data_new, *popt_subcircuit1)
    fluorescenceFit_continuous = steadystate_subcircuit1(OC14data_continuous,*popt_subcircuit1)
    gfpFit1 = fluorescenceFit[:10]; rfpFit1 = fluorescenceFit[10:20]
    gfpFit1_continuous = fluorescenceFit_continuous[:100]; rfpFit1_continuous = fluorescenceFit_continuous[100:200]

    plotFitvsData_croppedGreen(OC14_list1_green, OC14_list1_red,OC14_continuous, gfpExp_list1, rfpExp_list1, semGreen1, semRed1, gfpFit1_continuous,rfpFit1_continuous)

    gfpFit1_continuous_copy,rfpFit1_continuous_copy = gfpFit1_continuous,rfpFit1_continuous



    fluorescenceFit = steadystate_subcircuit3(OC14data_new, *popt_subcircuit3)
    fluorescenceFit_continuous = steadystate_subcircuit3(OC14data_continuous,*popt_subcircuit3)
    gfpFit3 = fluorescenceFit[:10]; rfpFit3 = fluorescenceFit[10:20]
    gfpFit3_continuous = fluorescenceFit_continuous[:100]; rfpFit3_continuous = fluorescenceFit_continuous[100:200]

    plotFitvsData_croppedGreen(OC14_list3_green, OC14_list3_red,OC14_continuous, gfpExp_list3, rfpExp_list3, semGreen3, semRed3, gfpFit3_continuous,rfpFit3_continuous)

    gfpFit3_continuous_copy,rfpFit3_continuous_copy = gfpFit3_continuous,rfpFit3_continuous


#############
###generate fit distribution#####
#############

#sample from the two subcircuit fits
print('Creating subcircuit 1 distribution fit')
sampled_parameters_subcircuit1 = np.random.multivariate_normal(popt_subcircuit1,pcov_subcircuit1*10, size= n_samples*2, check_valid='warn')#
print('Creating subcircuit 3 distribution fit')
sampled_parameters_subcircuit3 = np.random.multivariate_normal(popt_subcircuit3,pcov_subcircuit3*10, size=n_samples*2, check_valid='warn')#


#stack two arrays horizontally
sampled_parameters = np.hstack([sampled_parameters_subcircuit1,sampled_parameters_subcircuit3])

#filter parameters for negative values
print('Filtering parameters for negative values')
filtered_parameters = [p for p in tqdm(sampled_parameters) if np.all(p>0)]
print('Filtered!')
if len(filtered_parameters) < n_samples:
    print('not enough fitted samples')

if test==True:
    def steadystate_combined_subcircuits(OC14,Ve, Vf, Kvd, Kfe, Vc, Vd, Kda, Kce):
        F1 = gfp1_steadystate(OC14, Vf,Kvd)
        E1 = rfp1_steadystate(OC14,Ve,Vf,Kvd,Kfe)
        D3 = gfp3_steadystate(OC14, Vd)
        E3 = rfp3_steadystate(OC14,  Vc, Kda,Kce)
        stacked_fluorescence = [F1,E1, D3, E3]

        return stacked_fluorescence

    fig,ax = plt.subplots()
    ax2=ax.twinx()
    OC14_continuous = np.logspace(-3,2, 100)*HSLtransform
    best_fit_p = np.hstack([popt_subcircuit1, popt_subcircuit3])
    fluorescenceFit_continuous_list = []
    for p in filtered_parameters:
        fluorescenceFit = steadystate_combined_subcircuits(OC14_continuous, *p)
        fluorescenceFit_continuous = steadystate_combined_subcircuits(OC14_continuous, *p)
        fluorescenceSingleFit_continuous = steadystate_combined_subcircuits(OC14_continuous, *best_fit_p)

        
        ax2.plot(OC14_continuous, fluorescenceFit_continuous[0], c='darkseagreen', alpha=0.08)
        ax.plot(OC14_continuous, fluorescenceFit_continuous[1], c='lightcoral', alpha=0.08)
        ax2.scatter(OC14_list1_green,gfpExp_list1 , label='data', c='green')
        ax2.errorbar(OC14_list1_green,gfpExp_list1,yerr=semGreen1,c='green',fmt='o')
        ax.scatter(OC14_list1_red,rfpExp_list1 , label='data', c='red')
        ax.errorbar(OC14_list1_red,rfpExp_list1,yerr=semRed1,c='red',fmt='o')
        plt.xscale('log')
        fluorescenceFit_continuous_list.append(fluorescenceFit_continuous)


    # ax.legend(loc='center left') #upper right
    ax.set_ylabel('RFP / ($A_{600}$ $RFP_{basal})$', fontsize=15)
    ax.set_xscale('log')
    # ax2.legend(loc='center right') #upper left
    ax2.set_ylabel('GFP / ($A_{600}$ $GFP_{basal})$',fontsize=15)
    ax.set_xscale('log')
    ax.set_xlabel(f'3OHC14-HSL concentration (µM)',fontsize=15)
    ax.plot(OC14_continuous, rfpFit1_continuous_copy, c='red', alpha=1)
    ax2.plot(OC14_continuous, gfpFit1_continuous_copy, c='green', alpha=1)
    plt.show()

    gfpFit1_continuous_copy
    fig,ax = plt.subplots()
    ax2=ax.twinx()
    for p in filtered_parameters:
        fluorescenceFit = steadystate_combined_subcircuits(OC14_continuous, *p)
        fluorescenceFit_continuous = steadystate_combined_subcircuits(OC14_continuous, *p)
        fluorescenceSingleFit_continuous = steadystate_combined_subcircuits(OC14_continuous, *best_fit_p)

        
        ax2.plot(OC14_continuous, fluorescenceFit_continuous[2], c='darkseagreen', alpha=0.08)
        ax.plot(OC14_continuous, fluorescenceFit_continuous[3], c='lightcoral', alpha=0.08)
        ax2.scatter(OC14_list1_green,gfpExp_list3 , label='data', c='green')
        ax2.errorbar(OC14_list1_green,gfpExp_list3,yerr=semGreen3,c='green',fmt='o')
        ax.scatter(OC14_list1_red,rfpExp_list3 , label='data', c='red')
        ax.errorbar(OC14_list1_red,rfpExp_list3,yerr=semRed3,c='red',fmt='o')
        plt.xscale('log')

    # ax.legend(loc='center left') #upper right
    ax.set_ylabel('RFP / ($A_{600}$ $RFP_{basal})$', fontsize=15)
    ax.set_xscale('log')
    # ax2.legend(loc='center right') #upper left
    ax2.set_ylabel('GFP / ($A_{600}$ $GFP_{basal})$', fontsize=15)
    ax.set_xscale('log')
    ax.set_xlabel(f'3OHC14-HSL concentration (µM)', fontsize=15)

    ax.plot(OC14_continuous, rfpFit3_continuous_copy, c='red', alpha=1)
    ax2.plot(OC14_continuous, gfpFit3_continuous_copy, c='green', alpha=1)
    plt.show()

df = pd.DataFrame(filtered_parameters, columns= np.hstack([ paramNames_subcircuit1 , paramNames_subcircuit3]))
path = '/Users/mo2016/Documents/modelling/3954/paper/out/dose_response/dose_response_paper/'
pkl.dump(fluorescenceFit_continuous_list, open(path+'fluorescenceFit_continuous_list.pkl', 'wb'))
pkl.dump(OC14_continuous, open(path+'OC14_continuous.pkl', 'wb'))
pkl.dump(OC14_list1_green, open(path+'OC14_green.pkl', 'wb'))
pkl.dump(OC14_list1_red, open(path+'OC14_red.pkl', 'wb'))
pkl.dump(gfpExp_list1, open(path+'gfpExp_list1.pkl', 'wb'))
pkl.dump(rfpExp_list1, open(path+'rfpExp_list1.pkl', 'wb'))
pkl.dump(semGreen1, open(path+'semGreen1.pkl', 'wb'))
pkl.dump(semRed1, open(path+'semRed1.pkl', 'wb'))
pkl.dump(gfpExp_list3, open(path+'gfpExp_list3.pkl', 'wb'))
pkl.dump(rfpExp_list3, open(path+'rfpExp_list3.pkl', 'wb'))
pkl.dump(semGreen3, open(path+'semGreen3.pkl', 'wb'))
pkl.dump(semRed3, open(path+'semRed3.pkl', 'wb'))

pkl.dump(gfpFit1_continuous_copy, open(path+'gfp_bestFit_subcircuit1_continous.pkl', 'wb'))
pkl.dump(rfpFit1_continuous_copy, open(path+'rfp_bestFit_subcircuit1_continous.pkl', 'wb'))
pkl.dump(gfpFit3_continuous_copy, open(path+'gfp_bestFit_subcircuit3_continous.pkl', 'wb'))
pkl.dump(rfpFit3_continuous_copy, open(path+'rfp_bestFit_subcircuit3_continous.pkl', 'wb'))

import numpy, scipy.io
scipy.io.savemat(path + 'fluorescenceFit_continuous_list.mat', mdict={'fluorescenceFit_continuous_list': fluorescenceFit_continuous_list})
scipy.io.savemat(path + 'OC14_continuous.mat', mdict={'OC14_continuous': OC14_continuous})
scipy.io.savemat(path + 'OC14_green.mat', mdict={'OC14_green': OC14_list1_green})
scipy.io.savemat(path + 'OC14_red.mat', mdict={'OC14_red': OC14_list1_red})
scipy.io.savemat(path + 'gfpExp_list1.mat', mdict={'gfpExp_list1': gfpExp_list1})
scipy.io.savemat(path + 'rfpExp_list1.mat', mdict={'rfpExp_list1': rfpExp_list1})
scipy.io.savemat(path + 'semGreen1.mat', mdict={'semGreen1': semGreen1})
scipy.io.savemat(path + 'semRed1.mat', mdict={'semRed1': semRed1})
scipy.io.savemat(path + 'gfpExp_list3.mat', mdict={'gfpExp_list3': gfpExp_list3})
scipy.io.savemat(path + 'rfpExp_list3.mat', mdict={'rfpExp_list3': rfpExp_list3})
scipy.io.savemat(path + 'semGreen3.mat', mdict={'semGreen3': semGreen3})
scipy.io.savemat(path + 'semRed3.mat', mdict={'semRed3': semRed3})

scipy.io.savemat(path + 'gfp_bestFit_subcircuit1_continous.mat', mdict={'gfp_bestFit_subcircuit1_continous': gfpFit1_continuous_copy})
scipy.io.savemat(path + 'rfp_bestFit_subcircuit1_continous.mat', mdict={'rfp_bestFit_subcircuit1_continous': rfpFit1_continuous_copy})
scipy.io.savemat(path + 'gfp_bestFit_subcircuit3_continous.mat', mdict={'gfp_bestFit_subcircuit3_continous': gfpFit3_continuous_copy})
scipy.io.savemat(path + 'rfp_bestFit_subcircuit3_continous.mat', mdict={'rfp_bestFit_subcircuit3_continous': rfpFit3_continuous_copy})


fluorescenceFit_continuous_list
OC14_continuous
OC14_list1_green
OC14_list1_red
gfpExp_list1
rfpExp_list1
semGreen1
semRed1
gfpExp_list3
rfpExp_list3
semGreen3
semRed3

gfp_bestFit_subcircuit1_continous
rfp_bestFit_subcircuit1_continous
gfp_bestFit_subcircuit3_continous
rfp_bestFit_subcircuit3_continous




if test==True:
    sns.set_context("paper", rc={"axes.labelsize":20, "xtick.labelsize":15, "ytick.labelsize":15, "legend.fontsize":15, "axes.titlesize":20})
    sns.pairplot(df)

#############
###generate lhs distribution#####
#############
#%%
# pfitDict['Ve'] = np.amax([rfpExp_list3, rfpExp_list1]) - 1
pfitDict['nvd'] = nvd
pfitDict['nfe'] = nfe
pfitDict['nda'] = nda
pfitDict['nce'] = nce
pfitDict['Ve'] = Ve



#maximum production parameters (V*)
minV = 10;maxV=1000;minb=0.1;maxb=1
Va = {'name':'Va','distribution':'loguniform', 'min':minV/maxb, 'max':maxV/minb}
Vb = {'name':'Vb','distribution':'loguniform', 'min':minV/maxb, 'max':maxV/minb}
Ve = {'name':'Ve','distribution':'gaussian', 'mean':pfitDict['Ve'], 'noisetosignal':0.3}
V_parameters = [Va,Vb,Ve]



K1=0.0183; K2=0.0183
DUmin=0.1; DUmax=10; DVmin=0.1; DVmax=10
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

minK=0.1;maxK=250

Keb = {'name':'Keb','distribution':'loguniform', 'min':Kstar(muLVA_estimate,maxb,minK), 'max':Kstar(muLVA_estimate,minb,maxK)}

Kee = {'name':'Kee','distribution':'fixed','value':0.001}

Kub = {'name':'Kub','distribution':'loguniform', 'min':Kdiffstar(muU,KdiffpromMin,K2), 'max':Kdiffstar(muU,KdiffpromMax,K2)}


K_parameters = [ Kub, Keb, Kee]



#protein degradation parameters (mu)
muASV = {'name':'muASV','distribution':'fixed', 'value':muASV_estimate/muASV_estimate}
muLVA = {'name':'muLVA','distribution': 'gaussian','mean':muLVA_estimate /muASV_estimate, 'noisetosignal':0.1}
mu_parameters = [muLVA,muASV]


#cooperativity parameters (n)
nvd = {'name':'nvd','distribution':'fixed', 'value':pfitDict['nvd']}
nfe = {'name':'nfe','distribution':'fixed', 'value': pfitDict['nfe']} #ideally hihger but we keep lower to ensure numerics work
nda = {'name':'nda','distribution':'fixed', 'value':pfitDict['nda']}
nce = {'name':'nce','distribution':'fixed', 'value':pfitDict['nce']}
nub = {'name':'nub','distribution':'fixed', 'value':1}
nee = {'name':'nee','distribution':'fixed', 'value':4}
neb = {'name':'neb','distribution':'fixed', 'value':4}
nfe = {'name':'nfe','distribution':'fixed', 'value':8}
n_parameters = [nub,nee,neb,nvd,nda,nce,nfe]


createParams=False
if createParams == True:
    parameterDictList = D_parameters  + V_parameters + K_parameters + mu_parameters + n_parameters
    stackedDistributions = preLhs(parameterDictList)
    lhsDist = lhs(stackedDistributions,n_samples)
    lhsDist_df = pd.DataFrame(data = lhsDist, columns=[parameter['name'] for parameter in parameterDictList])

#%%

#############
###concatenate lhs and fit#####
#############
df = df.drop(columns='Ve1')
lhsDistFit_df=pd.concat([lhsDist_df, df], axis=1)
lhsDistFit_df = lhsDistFit_df.dropna()
if test==True:
    for column in lhsDistFit_df.columns:
        plt.hist(lhsDistFit_df[column],bins=20)
        plt.xscale('log')
        plt.title(column)
        plt.show()
        
#%%
#############
###check balance of distributionst#####
#############

Km_list = ['Kda', 'Kub', 'Keb', 'Kvd', 'Kfe',  'Kce' ]
KtoV = {'Kda': 'Vd', 'Kub': 'Va', 'Keb': 'Ve', 'Kvd': 'Vb', 'Kfe': 'Vf','Kce': 'Vc' }

def checkBalance(par_dict):
    balanceDict = {}
    for Km in Km_list:
        # print(Km)
        Vx =par_dict[KtoV[Km]]
        Kxy = par_dict[Km]
        if Kxy >= 1 and Kxy <= Vx:
            balanceDict[Km] = 'Balanced'

        elif Kxy > 0.1 and Kxy < Vx*10:
            balanceDict[Km] ='Semi balanced'
        elif Kxy <= 0.1 or Kxy >= Vx*10:
            balanceDict[Km] ='Not balanced'
        else:
            print('ERROR!!!!!!!!!')

    if 'Not balanced' in balanceDict.values():
        return 'Not balanced'
    elif 'Semi balanced'  in balanceDict.values():
        # print('semibalanced')
        return 'Semi balanced'
    elif all(x == 'Balanced' for x in balanceDict.values()):
        # print('Balanced')
        return 'Balanced'
    

createBalancedParams=False
if createBalancedParams == True:
    balancedDf = pd.DataFrame()
    semiBalancedDf = pd.DataFrame()
    notBalancedDf = pd.DataFrame()

    balanceList = []    
    for parID in tqdm(lhsDistFit_df.index):
        par_dict = lhsDistFit_df.loc[parID].to_dict()
        balanceList.append(checkBalance(par_dict))
    lhsDistFit_df['balance'] = balanceList
    
    #separate 3df
    balancedDf = lhsDistFit_df[lhsDistFit_df['balance']=='Balanced']
    semiBalancedDf= lhsDistFit_df[lhsDistFit_df['balance']=='Semi balanced']
    notBalancedDf = lhsDistFit_df[lhsDistFit_df['balance']=='Not balanced']
    
    lhsDistFit_df_balancesemibalance = lhsDistFit_df[lhsDistFit_df['balance']!='Not balanced']
    print('Lenght of balancedsemibalance df: %r'%len(lhsDistFit_df_balancesemibalance))
    pkl.dump(lhsDistFit_df, open(modellingpath + '/3954/paper/input/fitted_parameterfiles/df_circuit%r_variant%s_%rparametersets.pkl'%(circuit_n,variant,n_samples), 'wb'))
    pkl.dump(lhsDistFit_df_balancesemibalance, open(modellingpath + '/3954/paper/input/fitted_parameterfiles/df_circuit%r_variant%s_%rparametersets_balancedSemiBalanced.pkl'%(circuit_n,variant,n_samples), 'wb'))


print(len(balancedDf), len(semiBalancedDf), len(notBalancedDf), len(lhsDistFit_df_balancesemibalance))
# %%
