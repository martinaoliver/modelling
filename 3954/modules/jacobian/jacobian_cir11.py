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
par_dict = {'bmB': 0.1, 'bmD': 0.1, 'bmE': 0.1, 'bmF': 0.1, 'VmB': 48.6876001174938, 'VmD': 100.98190135112577, 'VmE': 318.03253736493616, 'VmF': 100.5358411892441, 'kemb': 14.067508952943552, 'khmd': 0.4060457100618679, 'kfme': 1.6136015452754962, 'keme': 64.48049019294464, 'kD': 200.2363920717794, 'aB': 1.3, 'aD': 1.3, 'aE': 1.3, 'aF': 1.3, 'alphaH': 10, 'muRNA': 12, 'muLVA': 1.08, 'muH': 1, 'R': 100, 'd_H': 0.5, 'n': 2}
# general_df = pickle.load(open(modelling_home + '/3954/parameter_space_search/code/interesting_parIDs/interesting_parIDs.pkl', "rb"))
# par_dict = general_df.iloc[0]



for key,value in par_dict.items():
    symbol_par_dict[key] = symbols('self.'+ key)

eq = circuit11_eq(symbol_par_dict)

mB,H,mD,mE,mF,B,D,E,F,wvn= symbols('mB'),symbols('H'), symbols('mD'), symbols('mE'), symbols('mF'), symbols('B'), symbols('D'), symbols('E'), symbols('F'), symbols('wvn')

functions = Matrix([eq.function_list([mB,H,mD,mE,mF,B,D,E,F],wvn)[x] for x in range(9)])#,eq.function_list([B,H,D,E,F]), eq.function_list([B,H,D,E,F]), eq.function_list([B,H,D,E,F]),eq.function_list([B,H,D,E,F])])

variables = Matrix([mB,H,mD,mE,mF,B,D,E,F])

jac = functions.jacobian(variables)
jac = np.array(jac)

J_list = [ 'JmB','JH','JmD','JmE','JmF','JB','JD','JE','JF' ]

count = 0
for J in (J_list):
    print(J_list[count] + ' = ' + str(list(jac[count])))
    count+=1
