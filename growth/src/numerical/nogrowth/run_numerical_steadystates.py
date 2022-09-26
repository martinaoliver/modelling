#############
###paths#####
#############
import sys
import os

from importlib_metadata import distribution
pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')
#############

from numerical.cn_nogrowth import cn_nogrowth
from numerical.cn_plot import plot1D, surfpattern
import pickle
import matplotlib.pyplot as plt

#system parameters
circuit_n = 'turinghill'
variant=0 
n_param_sets = 2000000
mechanism = 'nogrowth'

#dataframes with parameters and lsa output
# df= pickle.load( open(modellingpath + "/growth/input/parameterfiles/df_%s_variant%r_%rparametersets.pkl"%(circuit_n,variant,n_param_sets), "rb"))
multiple_lsa_df= pickle.load( open(modellingpath + "/growth/out/analytical/lsa_dataframes/lsa_df_%s_variant%r_%rparametersets.pkl"%(circuit_n,variant,n_param_sets), "rb"))
turing_df= pickle.load( open(modellingpath + '/growth/out/analytical/turing/turing_df_%s_variant%r_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), "rb"))

#solver parameters
L=50; x_gridpoints=5; J=L*x_gridpoints;I=J 
T=5000; t_gridpoints = 30; N=T*t_gridpoints #Number of timepoints

filename = lambda parID: '%s_variant%s_%s_ID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant,mechanism,parID,L,J,T,N)

for parID in turing_df.index:
    par_dict = multiple_lsa_df.loc[parID].loc[0].to_dict()
    ss_n = int(par_dict['ss_n'])

    for ss_n in range(ss_n):
        print(parID, ss_n)
        par_dict = multiple_lsa_df.loc[parID].loc[ss_n].to_dict()
        # print(par_dict)
        #run
        U_final,U_record, U0, x_grid, reduced_t_grid= cn_nogrowth(par_dict,L,J,T,N, circuit_n, tqdm_disable=True)

        pickle.dump(U_final, open(modellingpath + '/growth/out/numerical/%s/%s/data/2Dfinal_%s_ss%s.pkl'%(circuit_n,mechanism,filename(parID),ss_n), 'wb'))
        pickle.dump(U_record, open(modellingpath + '/growth/out/numerical/%s/%s/data/2Drecord_%s_ss%s.pkl'%(circuit_n,mechanism,filename(parID),ss_n), 'wb'))
