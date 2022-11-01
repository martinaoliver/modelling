#############
###paths#####
#############
import sys
import os

pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')
#############

circuit_n=12
variant= 0
nsamples=1000000

df= pickle.load( open(modellingpath + '/3954/paper/input/lhs_parameterfiles/df_circuit%s_variant%s_%rparametersets.pkl'%(circuit_n,variant,nsamples), "rb" ) )
df['steadystates']