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

test=False


#############
###plot data#####
#############
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


#############
###load dataset#####
#############
HSLtransform = 0.14*10**3



doseResponseExp1 = pkl.load(open('input/liquid_culture/curatedData/Jure_subcircuit1_circuit14_doseResponseOC14_0.5ATC.pkl','rb'))
OC14_list1= np.array(doseResponseExp1['OC14'])*HSLtransform; gfpExp_list1 = list(doseResponseExp1['mean_gfp']); rfpExp_list1 = list(doseResponseExp1['mean_rfp'])
semGreen1 = doseResponseExp1['std_gfp']; semRed1 = doseResponseExp1['std_rfp']
plotData(OC14_list1, rfpExp_list1, gfpExp_list1, semRed1, semGreen1)



doseResponseExp3 = pkl.load(open('input/liquid_culture/curatedData/Jure_subcircuit3_circuit14_doseResponseOC14_0.5ATC.pkl','rb'))
OC14_list3= np.array(doseResponseExp3['OC14'])*HSLtransform; gfpExp_list3 = list(doseResponseExp3['mean_gfp']); rfpExp_list3 = list(doseResponseExp3['mean_rfp'])
semGreen3 = doseResponseExp3['std_gfp']; semRed3 = doseResponseExp3['std_rfp']
plotData(OC14_list3, rfpExp_list3, gfpExp_list3, semRed3, semGreen3)
print(np.amax([rfpExp_list3, rfpExp_list1]))




#############
###fit dataset#####
#############

nvd = 2
nfe = 10
nda=2
nce=3

Ve = np.amax([rfpExp_list3, rfpExp_list1]) - 1


def gfp1_steadystate(OC14, Vf,Kvd): 
    muv = 0.0225 ; kv =  0.0183 ;
    F1 = 1 + Vf*(1/(1+((muv*Kvd)/(kv*OC14 + 1e-8))**nvd ))
    return F1


def rfp1_steadystate(OC14, Vf,Kvd,Kfe): 
    muv = 0.0225 ; kv =  0.0183 ;
    E1 = 1 + Ve*(1/(1+((1 + Vf*(1/(1+((muv*Kvd)/(kv*OC14 + 1e-8))**nvd )))/(Kfe+1e-8))**nfe))

    return E1

def gfp3_steadystate(OC14,  Vd,Kvd): 
    muv = 0.0225 ; kv =  0.0183 ;
    D3 = 1 + Vd*(1/(1+((muv*Kvd)/(kv*OC14 + 1e-8))**nvd ))
    return D3

def bfp3_steadystate(D,Vc,Kda):
    C3 = 1 + Vc*(1/(1+((D/(Kda+1e-8))**nda)))
    return C3


def rfp3_steadystate(D,Vc, Kda,Kce): 

    E3 = 1 + Ve*(1/(1+((bfp3_steadystate(D, Vc, Kda)/(Kce+1e-8))**nce)))
    return E3

OC14_continuous = np.logspace(-3,2, 100)*HSLtransform


def steadystate(OC14,Vc,Vd,Vf, Kvd,Kda, Kfe, Kce):
  
  if len(OC14) == 22:
      gaps = [5,5,6,6]
  else:
        gaps = [int(len(OC14)/4)]*4
  F1 = gfp1_steadystate(OC14[:np.sum(gaps[:1])],  Vf,Kvd)
  E1 = rfp1_steadystate(OC14[np.sum(gaps[:1]):np.sum(gaps[:2])],Vf,Kvd,Kfe)
  D3 = gfp3_steadystate(OC14[np.sum(gaps[:2]):np.sum(gaps[:3])], Vd,Kvd)
  E3 = rfp3_steadystate(OC14[np.sum(gaps[:3]):np.sum(gaps[:4])], Vc, Kda,Kce)
  FE = np.hstack([F1,E1, D3, E3])
  return FE


fluorescenceData = np.hstack([gfpExp_list1,rfpExp_list1, gfpExp_list3,rfpExp_list3])
OC14data_new = np.hstack([OC14_list1,OC14_list1, OC14_list3,OC14_list3])
OC14data_continuous= np.hstack([OC14_continuous,OC14_continuous, OC14_continuous,OC14_continuous])
semStacked= np.hstack([semGreen1,semRed1, semGreen3,semRed3])

popt, pcov = curve_fit(f=steadystate, xdata=OC14data_new, ydata=fluorescenceData )
# popt, pcov = curve_fit(f=steadystate, xdata=OC14data_new, ydata=fluorescenceData ,sigma =semStacked)

paramNames = ['Vc','Vd','Vf', 'Kvd','Kda', 'Kfe', 'Kce']
pfitDict = {}
for param in popt:
    pfitDict[paramNames[popt.tolist().index(param)]] = param

# if test==False:
print(pfitDict)
fluorescenceFit = steadystate(OC14data_new, *popt)
fluorescenceFit_continuous = steadystate(OC14data_continuous, *popt)
gfpFit1 = fluorescenceFit[:10]; rfpFit1 = fluorescenceFit[10:20]; gfpFit3 = fluorescenceFit[20:30]; rfpFit3 = fluorescenceFit[30:40]
gfpFit1_continuous = fluorescenceFit_continuous[:100]; rfpFit1_continuous = fluorescenceFit_continuous[100:200]; gfpFit3_continuous = fluorescenceFit_continuous[200:300]; rfpFit3_continuous = fluorescenceFit_continuous[300:400]

plotFitvsData(OC14_list1,OC14_continuous, gfpExp_list1, rfpExp_list1, semGreen1, semRed1, gfpFit1_continuous,rfpFit1_continuous)

plotFitvsData(OC14_list3,OC14_continuous, gfpExp_list3, rfpExp_list3, semGreen3, semRed3, gfpFit3_continuous,rfpFit3_continuous)
gfpFit1_continuous_copy,rfpFit1_continuous_copy, gfpFit3_continuous_copy,rfpFit3_continuous_copy = gfpFit1_continuous,rfpFit1_continuous, gfpFit3_continuous,rfpFit3_continuous 





#############
###produce distributions#####
#############

#produce multivariate gaussian
if test==True:
    sampled_parameters = np.random.multivariate_normal(popt,pcov*100, size=100, check_valid='warn')#
else:
    sampled_parameters = np.random.multivariate_normal(popt,pcov*100, size=12000000, check_valid='warn')#


def steadystateloss(OC14,Vc,Vd,Vf, Kvd,Kda, Kfe, Kce):

  F1 = gfp1_steadystate(OC14, Vf,Kvd)
  E1 = rfp1_steadystate(OC14,Vf,Kvd,Kfe)
  D3 = gfp3_steadystate(OC14, Vd,Kvd)
  E3 = rfp3_steadystate(OC14,  Vc, Kda,Kce)
  FE = [F1,E1, D3, E3]

  return FE



def func(p):
    loss_i = 0

    for count,OC14 in enumerate(OC14_list1):
      model= steadystateloss(OC14_list1[count],*p)

      loss_i+= ((gfpExp_list1[count] - model[0])**2 + (rfpExp_list1[count] - model[1])**2 +(gfpExp_list3[count] - model[2])**2 + (rfpExp_list3[count] - model[3])**2 )

    return loss_i



lossList = []
parameters_list = []
for p in tqdm(sampled_parameters):
   if np.all(p>0): #check for positive parameters
      if func(p)<100: #check for loss
        lossList.append(func(p))
        parameters_list.append(p)

df = pd.DataFrame(parameters_list, columns=paramNames)
# df['lossList'] = np.log(lossList)


print('Df with filtered parameters', len(parameters_list))

#############
###visualize distributions#####
#############
if test==True:
        
    fig,ax = plt.subplots()
    ax2=ax.twinx()
    OC14_continuous = np.logspace(-3,2,100)*HSLtransform

    for p in parameters_list:
        fluorescenceFit = steadystate(OC14data_new, *p)
        fluorescenceFit_continuous = steadystate(OC14data_continuous, *p)
        fluorescenceSingleFit_continuous = steadystate(OC14data_continuous, *popt)
        gfpFit1 = fluorescenceFit[:5]; rfpFit1 = fluorescenceFit[5:10]; gfpFit3 = fluorescenceFit[10:16]; rfpFit3 = fluorescenceFit[16:22]
        gfpFit1_continuous = fluorescenceFit_continuous[:100]; rfpFit1_continuous = fluorescenceFit_continuous[100:200]; gfpFit3_continuous = fluorescenceFit_continuous[200:300]; rfpFit3_continuous = fluorescenceFit_continuous[300:400]
        gfpSingleFit1_continuous = fluorescenceSingleFit_continuous[:100]; rfpSingleFit1_continuous = fluorescenceSingleFit_continuous[100:200]; gfpSingleFit3_continuous = fluorescenceSingleFit_continuous[200:300]; rfpSingleFit3_continuous = fluorescenceSingleFit_continuous[300:400]
    #     steadystate(OC14_continuous, *p)
        ax2.plot(OC14_continuous, gfpFit1_continuous, c='darkseagreen', alpha=0.08)
        ax.plot(OC14_continuous, rfpFit1_continuous, c='lightcoral', alpha=0.08)
        ax2.scatter(OC14_list1,gfpExp_list1 , label='data', c='green')
        ax2.errorbar(OC14_list1,gfpExp_list1,yerr=semGreen1,c='green',fmt='o')
        ax.scatter(OC14_list1,rfpExp_list1 , label='data', c='red')
        ax.errorbar(OC14_list1,rfpExp_list1,yerr=semRed1,c='red',fmt='o')
        plt.xscale('log')

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
    for p in parameters_list:
        fluorescenceFit = steadystate(OC14data_new, *p)
        fluorescenceFit_continuous = steadystate(OC14data_continuous, *p)
        gfpFit1 = fluorescenceFit[:5]; rfpFit1 = fluorescenceFit[5:10]; gfpFit3 = fluorescenceFit[10:16]; rfpFit3 = fluorescenceFit[16:22]
        gfpFit1_continuous = fluorescenceFit_continuous[:100]; rfpFit1_continuous = fluorescenceFit_continuous[100:200]; gfpFit3_continuous = fluorescenceFit_continuous[200:300]; rfpFit3_continuous = fluorescenceFit_continuous[300:400]
        
        ax2.plot(OC14_continuous, gfpFit3_continuous, c='darkseagreen', alpha=0.08)
        ax.plot(OC14_continuous, rfpFit3_continuous, c='lightcoral', alpha=0.08)
        ax2.scatter(OC14_list3,gfpExp_list3 , label='data', c='green')
        ax2.errorbar(OC14_list3,gfpExp_list3,yerr=semGreen3,c='green',fmt='o')
        ax.scatter(OC14_list3,rfpExp_list3 , label='data', c='red')
        ax.errorbar(OC14_list3,rfpExp_list3,yerr=semRed3,c='red',fmt='o')
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

    # sns.pairplot(df, hue="lossList", diag_kind='hist', diag_kws={'multiple': 'stack'},palette='rocket_r',corner=True)
    # sns.pairplot(df, hue="lossList",corner=True)
    # plt.show()

#############
###generate lhs distribution#####
#############

pfitDict['Ve'] = np.amax([rfpExp_list3, rfpExp_list1]) - 1
pfitDict['nvd'] = nvd
pfitDict['nfe'] = nfe
pfitDict['nda'] = nda
pfitDict['nce'] = nce

circuit_n=14
variant='fitted4'
#diffusion parameters


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

# [] at half activation parameters (K)
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


createParams=True
if createParams == True:
    # nsamples=len(df)
    nsamples=10000000
    if test==True:
        nsamples=100
    # nsamples=int(sys.argv[1])
    parameterDictList = D_parameters  + V_parameters + K_parameters + mu_parameters + n_parameters
    stackedDistributions = preLhs(parameterDictList)
    lhsDist = lhs(stackedDistributions,nsamples)
    lhsDist_df = pd.DataFrame(data = lhsDist, columns=[parameter['name'] for parameter in parameterDictList])



#############
###concatenate lhs and fit#####
#############

lhsDistFit_df=pd.concat([lhsDist_df, df], axis=1)
lhsDistFit_df = lhsDistFit_df.dropna()
if test==True:
    for column in lhsDistFit_df.columns:
        plt.hist(lhsDistFit_df[column],bins=20)
        plt.title(column)
        plt.show()
        

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
    

createBalancedParams=True
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
    
    # #concat to df
    # if len(balancedDf)<nsamples:
    #     balancedDf = pd.concat([balancedDf, balancedDfPre], ignore_index=True)
    # if len(semiBalancedDf)<nsamples:
    #     semiBalancedDf = pd.concat([semiBalancedDf, semiBalancedDfPre], ignore_index=True)
    # if len(notBalancedDf)<nsamples:
    #     notBalancedDf = pd.concat([notBalancedDf, notBalancedDfPre], ignore_index=True)

    lhsDistFit_df_balancesemibalance = lhsDistFit_df[lhsDistFit_df['balance']!='Not balanced']
    print('Lenght of balancedsemibalance df: %r'%len(lhsDistFit_df_balancesemibalance))
        # pkl.dump(balancedDf[:nsamples], open(modellingpath + '/3954/paper/input/balanced_parameterfiles/df_circuit%r_variant%s_%rparametersets_balanced.pkl'%(circuit_n,variant,nsamples), 'wb'))
        # pkl.dump(semiBalancedDf[:nsamples], open(modellingpath + '/3954/paper/input/balanced_parameterfiles/df_circuit%r_variant%s_%rparametersets_semiBalanced.pkl'%(circuit_n,variant,nsamples), 'wb'))
        # pkl.dump(notBalancedDf[:nsamples], open(modellingpath + '/3954/paper/input/balanced_parameterfiles/df_circuit%r_variant%s_%rparametersets_notBalanced.pkl'%(circuit_n,variant,nsamples), 'wb'))
    pkl.dump(lhsDistFit_df, open(modellingpath + '/3954/paper/input/fitted_parameterfiles/df_circuit%r_variant%s_%rparametersets.pkl'%(circuit_n,variant,nsamples), 'wb'))
    pkl.dump(lhsDistFit_df_balancesemibalance, open(modellingpath + '/3954/paper/input/fitted_parameterfiles/df_circuit%r_variant%s_%rparametersets_balancedSemiBalanced.pkl'%(circuit_n,variant,nsamples), 'wb'))


len(balancedDf), len(semiBalancedDf), len(notBalancedDf), len(lhsDistFit_df_balancesemibalance)