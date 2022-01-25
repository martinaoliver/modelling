
# %config Completer.use_jedi = False
import os
import sys
pwd = os.getcwd()
root = pwd.split("Documents", 1)[0]
modulepath = os.path.expanduser(root + 'Documents/modelling/6eq/modules')  # os.path.expanduser(path) : return the argument with an initial component of ~ or ~user replaced by that user’s home directory.
sys.path.append(modulepath)
path = os.path.expanduser(root + 'Documents/modelling/6eq/parameter_space_search') #path of project folder: folder where code, results and parameterfiles are found.
from sympy import *

import pickle
from class_circuit_eq import *
# print('load df')
# df= pickle.load( open(path +  '/parameterfiles/df_circuit2_variant0_1000000parametersets.pkl', "rb" ) )
# print('df loaded')
par_ID = 0
# par_dict = df.iloc[par_ID].to_dict()
par_dict = {}
symbol_par_dict = {}
# df = {'k1':0.05, 'k2': 0.01, 'k3':2, 'n1':3, 'n2':3, 'n3':1, 'beta':2.8,'Va': alpha, 'Vb':alpha, 'mua':1, 'mub':1, 'd_A':5e-05, 'd_B':0.0025}
par_dict = {'alpha':0.95, 'beta': -0.91, 'D':0.45, 'r3':0,  'r2':0, 'gamma':4}
# for key,value in par_dict.items():
#     symbol_par_dict[key] = symbols('self.'+ key)


for key,value in par_dict.items():
    symbol_par_dict[key] = symbols('self.'+ key)
# eq = circuit1_eq(symbol_par_dict)
# eq = circuit2_eq(symbol_par_dict)
# eq = circuit3_eq(symbol_par_dict)
eq = circuit7_eq(symbol_par_dict)
eq
for key,value in par_dict.items():
    symbol_par_dict[key] = symbols(key)


# # A,B,C,D,E,F,wvn= symbols('A'), symbols('B'), symbols('C'), symbols('D'), symbols('E'), symbols('F'), symbols('wvn')
# # A,B,C,D,E,F,M1,M2,wvn= symbols('A'), symbols('B'), symbols('C'), symbols('D'), symbols('E'), symbols('F'), symbols('M1'), symbols('M2'), symbols('wvn')
# # A,B,D,F,M1,M2,wvn= symbols('A'), symbols('B'), symbols('D'), symbols('F'), symbols('M1'), symbols('M2'), symbols('wvn')
# # A,B,D,F,wvn= symbols('A'), symbols('B'), symbols('D'), symbols('F'),  symbols('wvn')
A,B,wvn= symbols('A'), symbols('B'), symbols('wvn')
# print(eq)
# print
# functions = Matrix([eq.dAdt_f(A,B,C,D,E,F),eq.dBdt_f(A,B,C,D,E,F),eq.dCdt_f(A,B,C,D,E,F), eq.dDdt_f(A,B,C,D,E,F), eq.dEdt_f(A,B,C,D,E,F),eq.dFdt_f(A,B,C,D,E,F)])
# functions = Matrix([eq.dAdt_f(A,B,C,D,E,F,M1,M2),eq.dBdt_f(A,B,C,D,E,F,M1,M2),eq.dCdt_f(A,B,C,D,E,F,M1,M2), eq.dDdt_f(A,B,C,D,E,F,M1,M2), eq.dEdt_f(A,B,C,D,E,F,M1,M2),eq.dFdt_f(A,B,C,D,E,F,M1,M2),eq.diffusing_dM1dt_f(A,B,C,D,E,F,M1,M2,wvn),eq.diffusing_dM2dt_f(A,B,C,D,E,F,M1,M2,wvn)])
# functions = Matrix([eq.dAdt_f(A,B,D,F,M1,M2),eq.dBdt_f(A,B,D,F,M1,M2), eq.dDdt_f(A,B,D,F,M1,M2),eq.dFdt_f(A,B,D,F,M1,M2),eq.diffusing_dM1dt(A,B,D,F,M1,M2,wvn),eq.diffusing_dM2dt(A,B,D,F,M1,M2,wvn)])
# functions = Matrix([eq.diffusing_dAdt(A,B,D,F,wvn),eq.diffusing_dBdt(A,B,D,F,wvn), eq.dDdt_f(A,B,D,F),eq.dFdt_f(A,B,D,F)])
functions = Matrix([eq.diffusing_dAdt(A,B,wvn),eq.diffusing_dBdt(A,B,wvn)])

# variables = Matrix([A,B,C,D,E,F])
# variables = Matrix([A,B,C,D,E,F,M1,M2])
# variables = Matrix([A,B,D,F,M1,M2])
# variables = Matrix([A,B,D,F])
variables = Matrix([A,B])

jac = functions.jacobian(variables)
jac = np.array(jac)
# J_list = ['JA', 'JB', 'JC', 'JD','JE', 'JF' ]
# J_list = ['JA', 'JB', 'JC', 'JD','JE', 'JF','JM1', 'JM2' ]
# J_list = ['JA', 'JB', 'JD', 'JF','JM1', 'JM2' ]
# J_list = ['JA', 'JB', 'JD', 'JF']
J_list = ['JA', 'JB']
count = 0
for J in (J_list):
    print(J_list[count] + ' = ' + str(list(jac[count])))
    count+=1
