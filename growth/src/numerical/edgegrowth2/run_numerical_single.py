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

from numerical.cn_edgegrowth2 import cn_edgegrowth2
from numerical.cn_plot import plot1D, surfpattern
import pickle
import matplotlib.pyplot as plt

#system parameters
circuit_n = 'turinghill'
variant=0 
n_param_sets = 10

# par_dict = {'c1':0.1, 'c2':1,'c3':0.9,'c4':1, 'd_A': 1, 'd_B':10}
multiple_df= pickle.load( open(modellingpath + "/growth/out/analytical/lsa_dataframes/lsa_df_%s_variant%r_%rparametersets.pkl"%(circuit_n,variant,n_param_sets), "rb"))
df = multiple_df.xs(0, level=1)
#solver parameters
L=50; x_gridpoints=5; J=L*x_gridpoints;I=J 
T=5000; t_gridpoints = 30; N=T*t_gridpoints #Number of timepoints


# L=10; x_gridpoints=5; J=L*x_gridpoints;I=J 
# T=100; t_gridpoints = 30; N=T*t_gridpoints #Number of timepoints



parID= 10 #parameter set to use
par_dict = df.loc[parID].to_dict()

#run
U,U_record, U0, x_grid, reduced_t_grid, cellMatrix= cn_edgegrowth2(par_dict,L,J,T,N, circuit_n, rate=0.1, boundaryCoeff=2)
print(cellMatrix)
#plot
plot1D(U, savefig=False,filename='')
plt.show()
surfpattern(U_record, [x_grid, reduced_t_grid], 'linear',morphogen=1, rate=0, savefig=False,filename='',logResults=False,normalize=False)
surfpattern(U_record, [x_grid, reduced_t_grid], 'linear',  morphogen=0, rate=0, savefig=False,filename='',logResults=False,normalize=False)
