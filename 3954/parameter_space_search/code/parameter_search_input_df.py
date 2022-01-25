##########################
#########IMPORTS##########
##########################

#import module folder containing general functions used frequently
# modulepath = 'Users/mo2016/Documents/modelling/6eq/modules'
import os.path
import sys
modulepath = os.path.expanduser('~/Documents/modelling/6eq/modules')  # os.path.expanduser(path) : return the argument with an initial component of ~ or ~user replaced by that user’s home directory.
sys.path.append(modulepath)


#importing functions from that module folder
from parametercombination_analysis import *
import pickle
from datetime import date

date = str(date.today())
#path of folder containing code, results, parameter files
path = os.path.expanduser('~/Documents/modelling/6eq/parameter_space_search') #path of project folder: folder where code, results and parameterfiles are found.

#######################
#########CODE##########
#######################
circuit_n = int(sys.argv[1])
ID = int(sys.argv[2])

lhs_parameterfile = 'df_parameterfile_%s_ID%r.pkl' %(date,ID) #name of file containing dataframe with lhs parameter sets
lhs_df= pickle.load( open(path +  '/parameterfiles/df_parameterfile_%s/%s'%(date,lhs_parameterfile), "rb" ) )
big_turing_analysis_df(lhs_df,circuit_n) #performs turing analysis on dataframe and returns a dictionary with parameter set ID of turing I and turing II combinations
# pattern_class_dict = big_turing_analysis_df(lhs_df) #performs turing analysis on dataframe and returns a dictionary with parameter set ID of turing I and turing II combinations
# pickle.dump(pattern_class_dict, open(path + '/results/pattern_class_dict/pattern_class_dict_%s/pattern_class_dict_%s_ID%r.pkl' %(date,date,ID), 'wb')) #save dictionary into results folder
