import sys
import os
pwd = os.getcwd()
root = pwd.split("home", 1)[0]
modelling_home = root + 'home/Documents/modelling'
modelling_ephemeral = root + 'ephemeral/Documents/modelling'
modulepath = modelling_home + '/3954/modules'
sys.path.append(modulepath)
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
# par_dict = {'alpha':0.95, 'beta': -0.91, 'D':0.45, 'r3':0,  'r2':0, 'gamma':4}
par_dict = {'bb':0.95,'bd':0.95,'be':0.95,'bf':0.95, 'Vb': -0.91,'Vd': -0.91,'Ve': -0.91,'Vf': -0.91,'keb': -0.91, 'khd': -0.91,'kfe': -0.91,
'kee': -0.91,'muH':0.45,'muLVA':0.45,'alphaH':0.45,'d_H':0.45,'n':0.45,'kD':0.2,'R':1} #circuit10
# general_df = pickle.load(open(modelling_home + '/3954/parameter_space_search/code/interesting_parIDs/interesting_parIDs.pkl', "rb"))
# par_dict = general_df.iloc[0]



for key,value in par_dict.items():
    symbol_par_dict[key] = symbols('self.'+ key)
# print(symbol_par_dict)
# eq = circuit1_eq(symbol_par_dict)
# eq = circuit2_eq(symbol_par_dict)
# eq = circuit3_eq(symbol_par_dict)
# eq = circuit7_eq(symbol_par_dict)
# eq = circuit8_eq(symbol_par_dict)
# eq = circuit9_eq(symbol_par_dict)
eq = circuit10_eq(symbol_par_dict)
#
# for key,value in par_dict.items():
#     symbol_par_dict[key] = symbols(key)


# # A,B,C,D,E,F,wvn= symbols('A'), symbols('B'), symbols('C'), symbols('D'), symbols('E'), symbols('F'), symbols('wvn')
# # A,B,C,D,E,F,M1,M2,wvn= symbols('A'), symbols('B'), symbols('C'), symbols('D'), symbols('E'), symbols('F'), symbols('M1'), symbols('M2'), symbols('wvn')
# # A,B,D,F,M1,M2,wvn= symbols('A'), symbols('B'), symbols('D'), symbols('F'), symbols('M1'), symbols('M2'), symbols('wvn')
# # A,B,D,F,wvn= symbols('A'), symbols('B'), symbols('D'), symbols('F'),  symbols('wvn')
# A,B,wvn= symbols('A'), symbols('B'), symbols('wvn')
# B,D,E,F,wvn= symbols('B'), symbols('D'), symbols('E'), symbols('F'), symbols('wvn')
# A,C,D,E,F,wvn= symbols('A'),symbols('C'), symbols('D'), symbols('E'), symbols('F'), symbols('wvn')
B,H,D,E,F,wvn= symbols('B'),symbols('H'), symbols('D'), symbols('E'), symbols('F'), symbols('wvn')

# print(eq)
# print
# functions = Matrix([eq.dAdt_f(A,B,C,D,E,F),eq.dBdt_f(A,B,C,D,E,F),eq.dCdt_f(A,B,C,D,E,F), eq.dDdt_f(A,B,C,D,E,F), eq.dEdt_f(A,B,C,D,E,F),eq.dFdt_f(A,B,C,D,E,F)])
# functions = Matrix([eq.dAdt_f(A,B,C,D,E,F,M1,M2),eq.dBdt_f(A,B,C,D,E,F,M1,M2),eq.dCdt_f(A,B,C,D,E,F,M1,M2), eq.dDdt_f(A,B,C,D,E,F,M1,M2), eq.dEdt_f(A,B,C,D,E,F,M1,M2),eq.dFdt_f(A,B,C,D,E,F,M1,M2),eq.diffusing_dM1dt_f(A,B,C,D,E,F,M1,M2,wvn),eq.diffusing_dM2dt_f(A,B,C,D,E,F,M1,M2,wvn)])
# functions = Matrix([eq.dAdt_f(A,B,D,F,M1,M2),eq.dBdt_f(A,B,D,F,M1,M2), eq.dDdt_f(A,B,D,F,M1,M2),eq.dFdt_f(A,B,D,F,M1,M2),eq.diffusing_dM1dt(A,B,D,F,M1,M2,wvn),eq.diffusing_dM2dt(A,B,D,F,M1,M2,wvn)])
# functions = Matrix([eq.diffusing_dAdt(A,B,D,F,wvn),eq.diffusing_dBdt(A,B,D,F,wvn), eq.dDdt_f(A,B,D,F),eq.dFdt_f(A,B,D,F)])
# functions = Matrix([eq.diffusing_dAdt(A,B,wvn),eq.diffusing_dBdt(A,B,wvn)])
# functions = Matrix([eq.diffusing_dBdt_f([B,D,E,F],wvn), eq.dDdt_f([B,D,E,F]), eq.dEdt_f([B,D,E,F]),eq.dFdt_f([B,D,E,F])])
# functions = Matrix([eq.diffusing_dAdt_f([A,C,D,E,F],wvn), eq.dCdt_f([A,C,D,E,F]),eq.dDdt_f([A,C,D,E,F]), eq.dEdt_f([A,C,D,E,F]),eq.dFdt_f([A,C,D,E,F])])
# functions = Matrix([eq.function_list([B,H,D,E,F],wvn)[0],eq.function_list([B,H,D,E,F]), eq.function_list([B,H,D,E,F]), eq.function_list([B,H,D,E,F]),eq.function_list([B,H,D,E,F])])
functions = Matrix([eq.function_list([B,H,D,E,F],wvn)[x] for x in range(5)])#,eq.function_list([B,H,D,E,F]), eq.function_list([B,H,D,E,F]), eq.function_list([B,H,D,E,F]),eq.function_list([B,H,D,E,F])])

# variables = Matrix([A,B,C,D,E,F])
# variables = Matrix([A,B,C,D,E,F,M1,M2])
# variables = Matrix([A,B,D,F,M1,M2])
# variables = Matrix([A,B,D,F])
# variables = Matrix([A,B])
# variables = Matrix([B,D,E,F])
# variables = Matrix([A,C,D,E,F])
variables = Matrix([B,H,D,E,F])

jac = functions.jacobian(variables)
jac = np.array(jac)
# J_list = ['JA', 'JB', 'JC', 'JD','JE', 'JF' ]
# J_list = ['JA', 'JB', 'JC', 'JD','JE', 'JF','JM1', 'JM2' ]
# J_list = ['JA', 'JB', 'JD', 'JF','JM1', 'JM2' ]
# J_list = ['JA', 'JB', 'JD', 'JF']
# J_list = ['JA', 'JB']
# J_list = [ 'JB', 'JD','JE', 'JF' ]
# J_list = [ 'JA','JC', 'JD','JE', 'JF' ]
J_list = [ 'JB', 'JH','JD','JE', 'JF' ]

count = 0
for J in (J_list):
    print(J_list[count] + ' = ' + str(list(jac[count])))
    count+=1
