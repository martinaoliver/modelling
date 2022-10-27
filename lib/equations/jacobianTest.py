
###paths#####
import sys
import os
pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')
#############
libpath = modellingpath+'\lib'
sys.path.append(libpath)
#############
import unittest
import jacobian as s
import numpy as np
from numpy import testing
from sympy import *
import pickle as pkl

class Jacobian(unittest.TestCase):
    def test_jacobian(self):
        general_df = pkl.load(open(modellingpath + '/growth/input/parameterfiles/df_schnakenberg_variant0_10parametersets.pkl', "rb"))
        par_dict = general_df.iloc[1]
        for key,value in par_dict.items():
            setattr(self,key,value)

        A,B,wvn= symbols('A'), symbols('B'), symbols('wvn')
        resultingJacobian = s.jac
        expectedJacobian = np.array([[2*A*B*self.c3 - self.c2 - self.d_A*wvn**2, A**2*self.c3],[-2*A*B*self.c3, -A**2*self.c3 - self.d_B*wvn**2]])
        testing.assert_array_equal(expectedJacobian,resultingJacobian)

   
   
    # def testSize_cell_matrix_record(self):
    #     self.assertEqual(np.shape(s.cell_matrix_record),(s.J,s.J,s.N), 'Expected to be JJN' )

    # def test_check_neighbours(self):
    #     cell_matrix = np.ones((5,5))
    #     y = 2
    #     x = 3
    #     result = s.check_neighbours(cell_matrix, x, y)
        
    #     expected = np.array([[1,1,1],[1,np.nan,1],[1,1,1]])
    #     testing.assert_array_equal(result, expected)
if __name__ =='__main__':
    unittest.main()

    