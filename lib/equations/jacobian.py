
###paths#####
import sys
import os
pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')
#############
libpath = modellingpath+'\lib'
sys.path.append(libpath)
from sympy import *

import pickle
from equations.twonode_eq import * 
from equations.class_circuit_eq import * 

par_ID = 0
par_dict = {}
symbol_par_dict = {}
circuit_n=14
variant='0nd'
nsamples = 10
general_df = pickle.load(open(modellingpath +  '/3954/paper/input/lhs_parameterfiles/df_circuit%r_variant%s_%rparametersets.pkl'%(circuit_n,variant,nsamples), "rb"))
# general_df = pickle.load(open(modellingpath + '/growth/input/parameterfiles/df_schnakenberg_variant0_10parametersets.pkl', "rb"))
par_dict = general_df.iloc[1]


for key,value in par_dict.items():
    symbol_par_dict[key] = symbols('self.'+ key)
# print(symbol_par_dict)
# eq = circuit1_eq(symbol_par_dict)
# eq = circuit2(symbol_par_dict)
# eq = circuit3_eq(symbol_par_dict)
# eq = circuit7_eq(symbol_par_dict)
# eq = circuit8_eq(symbol_par_dict)
# eq = circuit9_eq(symbol_par_dict)
# eq = circuit10_eq(symbol_par_dict)
# eq = circuit12(symbol_par_dict)
eq = circuit13(symbol_par_dict)
eq = circuit14(symbol_par_dict)
# eq=turinghill(symbol_par_dict)
# eq=twonode(symbol_par_dict)
# eq=schnakenberg(symbol_par_dict)
# for key,value in par_dict.items():
#     symbol_par_dict[key] = symbols(key)


A,B,C,D,E,F,wvn= symbols('A'), symbols('B'), symbols('C'), symbols('D'), symbols('E'), symbols('F'), symbols('wvn')
# # A,B,C,D,E,F,M1,M2,wvn= symbols('A'), symbols('B'), symbols('C'), symbols('D'), symbols('E'), symbols('F'), symbols('M1'), symbols('M2'), symbols('wvn')
# # A,B,D,F,M1,M2,wvn= symbols('A'), symbols('B'), symbols('D'), symbols('F'), symbols('M1'), symbols('M2'), symbols('wvn')
# # A,B,D,F,wvn= symbols('A'), symbols('B'), symbols('D'), symbols('F'),  symbols('wvn')
# A,B,wvn= symbols('A'), symbols('B'), symbols('wvn')
# B,D,E,F,wvn= symbols('B'), symbols('D'), symbols('E'), symbols('F'), symbols('wvn')
# A,C,D,E,F,wvn= symbols('A'),symbols('C'), symbols('D'), symbols('E'), symbols('F'), symbols('wvn')
# B,H,D,E,F,wvn= symbols('B'),symbols('H'), symbols('D'), symbols('E'), symbols('F'), symbols('wvn')
# U,V,A,B,C,D,E,F,aTc,wvn = symbols('U'),symbols('V'),symbols('A'),symbols('B'),symbols('C'),symbols('D'),symbols('E'),symbols('F'),symbols('aTc'),symbols('wvn')
# U,V,A,B,C,D,E,F,wvn = symbols('U'),symbols('V'),symbols('A'),symbols('B'),symbols('C'),symbols('D'),symbols('E'),symbols('F'),symbols('wvn')
# print(eq)
# print

functions = Matrix([eq.dAdt_f([A,B,C,D,E,F], wvn),eq.dBdt_f([A,B,C,D,E,F], wvn),eq.dCdt_f([A,B,C,D,E,F]), eq.dDdt_f([A,B,C,D,E,F]), eq.dEdt_f([A,B,C,D,E,F]),eq.dFdt_f([A,B,C,D,E,F])])

# functions = Matrix([eq.dAdt_f(A,B,C,D,E,F,M1,M2),eq.dBdt_f(A,B,C,D,E,F,M1,M2),eq.dCdt_f(A,B,C,D,E,F,M1,M2), eq.dDdt_f(A,B,C,D,E,F,M1,M2), eq.dEdt_f(A,B,C,D,E,F,M1,M2),eq.dFdt_f(A,B,C,D,E,F,M1,M2),eq.diffusing_dM1dt_f(A,B,C,D,E,F,M1,M2,wvn),eq.diffusing_dM2dt_f(A,B,C,D,E,F,M1,M2,wvn)])
# functions = Matrix([eq.dAdt_f(A,B,D,F,M1,M2),eq.dBdt_f(A,B,D,F,M1,M2), eq.dDdt_f(A,B,D,F,M1,M2),eq.dFdt_f(A,B,D,F,M1,M2),eq.diffusing_dM1dt(A,B,D,F,M1,M2,wvn),eq.diffusing_dM2dt(A,B,D,F,M1,M2,wvn)])
# functions = Matrix([eq.diffusing_dAdt(A,B,D,F,wvn),eq.diffusing_dBdt(A,B,D,F,wvn), eq.dDdt_f(A,B,D,F),eq.dFdt_f(A,B,D,F)])
# interaction_matrix = np.array([[1,1],[-1,0]])
# functions = Matrix([eq.diffusing_dAdt_f([A,B],wvn, interaction_matrix),eq.diffusing_dBdt_f([A,B],wvn, interaction_matrix)])
# functions = Matrix([eq.diffusing_dAdt_f([A,B],wvn),eq.diffusing_dBdt_f([A,B],wvn)])
# functions = Matrix([eq.diffusing_dBdt_f([B,D,E,F],wvn), eq.dDdt_f([B,D,E,F]), eq.dEdt_f([B,D,E,F]),eq.dFdt_f([B,D,E,F])])
# functions = Matrix([eq.diffusing_dAdt_f([A,C,D,E,F],wvn), eq.dCdt_f([A,C,D,E,F]),eq.dDdt_f([A,C,D,E,F]), eq.dEdt_f([A,C,D,E,F]),eq.dFdt_f([A,C,D,E,F])])
# functions = Matrix([eq.function_list([B,H,D,E,F],wvn)[0],eq.function_list([B,H,D,E,F]), eq.function_list([B,H,D,E,F]), eq.function_list([B,H,D,E,F]),eq.function_list([B,H,D,E,F])])
# functions = Matrix([eq.function_list([B,H,D,E,F],wvn)[x] for x in range(5)])#,eq.function_list([B,H,D,E,F]), eq.function_list([B,H,D,E,F]), eq.function_list([B,H,D,E,F]),eq.function_list([B,H,D,E,F])])
# functions = Matrix([eq.function_list([U,V,A,B,C,D,E,F,aTc],wvn)[x] for x in range(9)])
# functions = Matrix([eq.dUdt_f([U,V,A,B,C,D,E,F],wvn), eq.dVdt_f([U,V,A,B,C,D,E,F],wvn), eq.dAdt_f([U,V,A,B,C,D,E,F]), eq.dBdt_f([U,V,A,B,C,D,E,F]), eq.dCdt_f([U,V,A,B,C,D,E,F]), eq.dDdt_f([U,V,A,B,C,D,E,F]), eq.dEdt_f([U,V,A,B,C,D,E,F]), eq.dFdt_f([U,V,A,B,C,D,E,F])])

variables = Matrix([A,B,C,D,E,F])

# variables = Matrix([A,B,C,D,E,F,M1,M2])
# variables = Matrix([A,B,D,F,M1,M2])
# variables = Matrix([A,B,D,F])
# variables = Matrix([A,B])
# variables = Matrix([B,D,E,F])
# variables = Matrix([A,C,D,E,F])
# variables = Matrix([B,H,D,E,F])
# variables = Matrix([U,V,A,B,C,D,E,F,aTc])
# variables = Matrix([U,V,A,B,C,D,E,F])

jac = functions.jacobian(variables)
jac = np.array(jac)
# print(jac)
J_list = ['JA', 'JB', 'JC', 'JD','JE', 'JF' ]
# J_list = ['JA', 'JB', 'JC', 'JD','JE', 'JF','JM1', 'JM2' ]
# J_list = ['JA', 'JB', 'JD', 'JF','JM1', 'JM2' ]
# J_list = ['JA', 'JB', 'JD', 'JF']
# J_list = ['JA', 'JB']
# J_list = [ 'JB', 'JD','JE', 'JF' ]
# J_list = [ 'JA','JC', 'JD','JE', 'JF' ]
# J_list = [ 'JB', 'JH','JD','JE', 'JF' ]
# J_list = [ 'JU', 'JV','JA','JB','JC','JD','JE', 'JF','JaTc' ]
# J_list = [ 'JU', 'JV','JA','JB','JC','JD','JE', 'JF']

count = 0
for J in (J_list):
    print(J_list[count] + ' = ' + str(list(jac[count])))
    count+=1
