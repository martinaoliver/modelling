import numpy as np
from numpy.random import normal
# - This file contains the definitions for every circuit used. Every circuit is composed of different
# molecular interactions.


# - The system is defined with differential equations for every specie where hill functions are used.
# - Concentrations of the morphogens and model parameters are inputs to these classes.

# - Every class contains (1) Differential equations and (2) jacobian matrix for the system.
# - The jacobian matrix is calculated in advance and defined here so it doesn't have to be
#calculated every time a parameter set is analysed (saves computational power).


class hill_functions():

    def __init__(self, par_dict):
        for key, value in par_dict.items():
            setattr(self, key, value)

    def noncompetitiveact(self, X, km,n):
        act = ((X / km) ** (n)) / (1 + (X / km) ** (n))
        act = (1 / (1 + (km / (X+1e-08)) ** (n)))
        return act

    def noncompetitiveinh(self, X, km,n):
        inh = 1 / (1 + (X / (km+1e-08)) ** (n))
        return inh


    def noncompetitivediffact(self, X, km,n, kdiff,mudiff):
        act = (1 / (1 + ((mudiff*km) / (kdiff*X)) ** (n)))
        return act
    
class subcircuit1_circuit14(hill_functions):#nodeC with inhibitions
#from circuit14

    def __init__(self,par_dict,stochasticity=0):
        for key,value in par_dict.items():
            setattr(self,key,value)
        setattr(self, 'stochasticity', stochasticity)



    def ddt(self,species_list,t, B,B1,wvn=0):
        # C,E,F = species_list
        E,F = species_list
        # dcdt= self.muLVA*(1 + self.VC - C ) 
        # dedt= self.muLVA*(1 + self.VE*self.noncompetitiveinh(C,self.Kce,self.nce)*self.noncompetitiveinh(F,self.Kfe,self.nfe)*self.noncompetitiveact(E,self.Kee,self.nee) - E ) 
        dedt= self.muLVA*(1 + self.VE*self.noncompetitiveinh(F,self.Kfe,self.nfe)*self.noncompetitiveact(E,self.Kee,self.nee) - E ) 
        dfdt= self.muLVA*(1 + self.VF*self.noncompetitiveact(B,self.Kbd,self.nbd) - F ) 
        # print(])
        # return dcdt, dedt, dfdt
        return dedt, dfdt
    

class subcircuit2_circuit14(hill_functions): #node B
#from circuit14

    def __init__(self,par_dict,stochasticity=0):
        for key,value in par_dict.items():
            setattr(self,key,value)
        setattr(self, 'stochasticity', stochasticity)



    def ddt(self,species_list,t, A,A1,wvn=0):
        B,D = species_list
        dbdt= self.muLVA*(1 + self.VB*self.noncompetitiveact(A,self.Kab, self.nab) - B ) -  B*self.Dr*wvn**2
        dddt= self.muLVA*(1 + self.VD*self.noncompetitiveact(B,self.Kbd,self.nbd) - D ) 
        
        return dbdt,dddt
    

class subcircuit3_circuit14(hill_functions): #node A and C
#from circuit14

    def __init__(self,par_dict,stochasticity=0):
        for key,value in par_dict.items():
            setattr(self,key,value)
        setattr(self, 'stochasticity', stochasticity)



    def ddt(self,species_list,t, B, B1, wvn=0):
        C,D,E = species_list
        dcdt= self.muLVA*(1 + self.VC*self.noncompetitiveinh(D,self.Kda,self.nda) - C ) 
        dddt= self.muLVA*(1 + self.VD*self.noncompetitiveact(B,self.Kbd,self.nbd) - D ) 
        dedt= self.muLVA*(1 + self.VE*self.noncompetitiveinh(C,self.Kce,self.nce)*self.noncompetitiveact(E,self.Kee,self.nee) - E ) 

        return dcdt,dddt, dedt
    

class subcircuit4_circuit14(hill_functions): #node B and C
#from circuit14

    def __init__(self,par_dict,stochasticity=0):
        for key,value in par_dict.items():
            setattr(self,key,value)
        setattr(self, 'stochasticity', stochasticity)



    def ddt(self,species_list,t, A,A1,wvn=0):
        B,C,D,E = species_list
        dbdt= self.muLVA*(1 + self.VB*self.noncompetitiveact(A,self.Kab, self.nab)*self.noncompetitiveinh(E,self.Keb, self.neb) - B ) -  B*self.Dr*wvn**2
        dcdt= self.muLVA*(1 + self.VC - C ) 
        dddt= self.muLVA*(1 + self.VD*self.noncompetitiveact(B,self.Kbd,self.nbd) - D ) 
        dedt= self.muLVA*(1 + self.VE*self.noncompetitiveinh(C,self.Kce,self.nce)*self.noncompetitiveact(E,self.Kee,self.nee) - E ) 

        return dbdt,dcdt,dddt,dedt




class subcircuitAC_1(hill_functions):

    def __init__(self,par_dict):
        for key,value in par_dict.items():
            setattr(self,key,value)


    def ddt(self,species_list,t):
        E= species_list
        C=10
        # dcdt= self.bc-self.mulva*C
        dedt= self.be+self.Ve*self.noncompetitiveinh(C,self.kce)*self.noncompetitiveact(E,self.kee)-self.mulva*E
        return dedt






