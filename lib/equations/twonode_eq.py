
# - This file contains the definitions for every 2node circuit used. Every circuit is composed of different
# molecular interactions.

# - The system is defined with differential equations for every specie where hill functions are used.
# - Concentrations of the morphogens and model parameters are inputs to these classes.

# - Every class contains (1) Differential equations and (2) jacobian matrix for the system.
# - The jacobian matrix is calculated in advance and defined here so it doesn't have to be
#calculated every time a parameter set is analysed (saves computational power).
#############

###paths#####
import sys
import os
pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')
#############

###imports###
import numpy as np
#############


class hill_functions():

    def __init__(self, par_dict):
        for key, value in par_dict.items():
            setattr(self, key, value)

    def noncompetitiveact(self, U, km,n=2):
        act = ((U / km) ** (n)) / (1 + (U / km) ** (n))
        return act

    def noncompetitiveinh(self, U, km,n=2):
        inh = 1 / (1 + (U / km) ** (n))
        return inh
    
    def noncompetitive_interaction(self,U,km,interaction,n=2):
        if interaction == -1:
            return 1 / (1 + (U / km) ** (n))
        if interaction == +1:
            return ((U / km) ** (n)) / (1 + (U / km) ** (n))
        if interaction == 0:
            return 1


class turinghill(hill_functions):

    def __init__(self,par_dict,stochasticity=0):
        for key,value in par_dict.items():
            setattr(self,key,value)
        setattr(self, 'stochasticity', stochasticity)


    def dAdt_f(self,species_list):
        A,B= species_list
        dadt= self.ba+self.Va*self.noncompetitive_interaction(A,self.kaa, 1)*self.noncompetitive_interaction(B,self.kba, -1)-self.mua*A
        return dadt


    def dBdt_f(self,species_list):
        A,B= species_list
        dbdt= self.bb+self.Vb*self.noncompetitive_interaction(A,self.kab,1)*self.noncompetitive_interaction(B,self.kbb,0)-self.mub*B
        return dbdt

    function_list = [dAdt_f,dBdt_f]


    def diffusing_dAdt_f(self,species_list,wvn):
        A,B= species_list
        return self.dAdt_f(species_list) - A*self.d_A*wvn**2


    def diffusing_dBdt_f(self,species_list,wvn):
        A,B= species_list
        return self.dBdt_f(species_list) - B*self.d_B*wvn**2



    def getJacobian(self,x,wvn):
        A,B=x
        JA = [-2*A**3*self.Va/(self.kaa**4*(A**2/self.kaa**2 + 1)**2*(B**2/self.kba**2 + 1)) + 2*A*self.Va/(self.kaa**2*(A**2/self.kaa**2 + 1)*(B**2/self.kba**2 + 1)) - self.d_A*wvn**2 - self.mua, -2*A**2*B*self.Va/(self.kaa**2*self.kba**2*(A**2/self.kaa**2 + 1)*(B**2/self.kba**2 + 1)**2)]
        JB = [-2*A**3*self.Vb/(self.kab**4*(A**2/self.kab**2 + 1)**2) + 2*A*self.Vb/(self.kab**2*(A**2/self.kab**2 + 1)), -self.d_B*wvn**2 - self.mub]
        return np.array([JA, JB])
        
class twonode(hill_functions):

    def __init__(self,par_dict,stochasticity=0):
        for key,value in par_dict.items():
            setattr(self,key,value)
        setattr(self, 'stochasticity', stochasticity)

    # interaction_matrix = np.array([[1,1],[-1,0]])

    def dAdt_f(self,species_list, interaction_matrix):
        A,B= species_list
        dadt= self.ba+self.Va*self.noncompetitive_interaction(A,self.kaa, interaction_matrix[0,0])*self.noncompetitive_interaction(B,self.kba, interaction_matrix[1,0])-self.mua*A
        return dadt


    def dBdt_f(self,species_list,interaction_matrix ):
        A,B= species_list
        dbdt= self.bb+self.Vb*self.noncompetitive_interaction(A,self.kab, interaction_matrix[0,1])*self.noncompetitive_interaction(B,self.kbb, interaction_matrix[1,1])-self.mub*B
        return dbdt

    def diffusing_dAdt_f(self,species_list,wvn, interaction_matrix):
        A,B= species_list
        dadt= self.ba+self.Va*self.noncompetitive_interaction(A,self.kaa, interaction_matrix[0,0])*self.noncompetitive_interaction(B,self.kba, interaction_matrix[1,0])-self.mua*A - A*self.d_A*wvn**2
        return dadt


    def diffusing_dBdt_f(self,species_list,wvn, interaction_matrix):
        A,B= species_list
        dbdt= self.bb+self.Vb*self.noncompetitive_interaction(A,self.kab, interaction_matrix[0,1])*self.noncompetitive_interaction(B,self.kbb, interaction_matrix[1,1])-self.mub*B - B*self.d_B*wvn**2
        return dbdt

    function_list = [dAdt_f,dBdt_f]

    #####this jacobian below is for turing circuit.... think how to get jacobian dependent on interaction matrix.
    # def getJacobian(self,x,wvn):
    #     A,B=x
    #     return np.array([JA, JB])

class schnakenberg():


    def __init__(self,par_dict,stochasticity=0):
        for key,value in par_dict.items():
            setattr(self,key,value)
        setattr(self, 'stochasticity', stochasticity)

# def schnakenberg(u,c=[0.1,1,0.9,1]):
#     c1,cm1,c2,c3 = c
#     f_u0 = c1 - cm1*u[0] + c3*(u[0]**2)*u[1]
#     f_u1 = c2 - c3*(u[0]**2)*u[1]
#     return f_u0,f_u1
    
    def dAdt_f(self,species_list):
        A,B= species_list
        dadt=  self.c1 - self.c2*A + self.c3*(A**2)*B
        return dadt


    def dBdt_f(self,species_list):
        A,B= species_list
        dbdt= self.c4 - self.c3*(A**2)*B
        return dbdt

    function_list = [dAdt_f,dBdt_f]


    def diffusing_dAdt_f(self,species_list,wvn):
        A,B= species_list
        return self.dAdt_f(species_list) - A*self.d_A*wvn**2
        # dadt= self.ba+self.Va*self.noncompetitive_interaction(A,self.kaa, 1)*self.noncompetitive_interaction(B,self.kba, -1)-self.mua*A - A*self.d_A*wvn**2
        # return dadt


    def diffusing_dBdt_f(self,species_list,wvn):
        A,B= species_list
        return self.dBdt_f(species_list) - B*self.d_B*wvn**2



    def getJacobian(self,x,wvn):
        A,B=x
        JA = [2*A*B*self.c3 - self.c2 - self.d_A*wvn**2, A**2*self.c3]
        JB = [-2*A*B*self.c3, -A**2*self.c3 - self.d_B*wvn**2]
        return np.array([JA, JB])
