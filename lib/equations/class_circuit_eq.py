import numpy as np
from numpy.random import normal
import numba

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
        # act = ((X / km) ** (n)) / (1 + (X / km) ** (n))
        act = (1 / (1 + (km / (X + 1e-8)) ** (n)))
        return act

    def noncompetitiveinh(self, X, km,n):
        inh = 1 / (1 + (X / (km + 1e-8) ) ** (n))
        return inh


    def noncompetitivediffact(self, X, km,n, kdiff,mudiff):
        act = (1 / (1 + ((mudiff*km) / (kdiff*X)) ** (n)))
        return act
 





    #
    # def noncompetitiveact(self, U, km):
    #     act = ((U / km) ** (self.n)) / (1 + (U / km) ** (self.n))
    #     return act
    #
    # def noncompetitiveinh(self, U, km):
    #     inh = 1 / (1 + (U / km) ** (self.n))
    #     return inh
#3954 with positive feedback on nodeA
class circuit1(hill_functions):

        def __init__(self,par_dict):
            for key,value in par_dict.items():
                setattr(self,key,value)

        def dAdt_f(self,A,B,C,D,E,F):
            dadt= self.ba+self.Va*self.noncompetitiveact(A,self.kaa)*self.noncompetitiveinh(D,self.kda)-self.mua*A 
            return dadt


        def dBdt_f(self,A,B,C,D,E,F):
            dbdt= self.bb+self.Vb*self.noncompetitiveact(A,self.kaa)*self.noncompetitiveinh(E,self.keb)-self.mub*B
            return dbdt


        def dCdt_f(self,A,B,C,D,E,F):
            dcdt= self.bc+self.Vc*self.noncompetitiveact(A,self.kaa)*self.noncompetitiveinh(D,self.kda)-self.mulva*C
            return dcdt

        def dDdt_f(self,A,B,C,D,E,F):
            dddt= self.bd+self.Vd*self.noncompetitiveact(B,self.kbd)-self.mulva*D
            return dddt


        def dEdt_f(self,A,B,C,D,E,F):
            dedt= self.be+self.Ve*self.noncompetitiveinh(C,self.kce)*self.noncompetitiveinh(F,self.kfe)*self.noncompetitiveact(E,self.kee)-self.mulva*E
            return dedt

        def dFdt_f(self,A,B,C,D,E,F):
            dfdt= self.bf+self.Vf*self.noncompetitiveact(B,self.kbd)-self.mulva*F
            return dfdt

        function_list = [dAdt_f, dBdt_f, dCdt_f, dDdt_f, dEdt_f, dFdt_f]


        def getJacobian_(self,x,wvn):
            A, B, C, D, E, F = x

            JA = [-self.mua - self.Va * self.n * (A / self.kaa) ** (2 * self.n) / (
                        A * ((A / self.kaa) ** self.n + 1) ** 2 * ((D / self.kda) ** self.n + 1)) + self.Va * self.n * (
                              A / self.kaa) ** self.n / (
                              A * ((A / self.kaa) ** self.n + 1) * ((D / self.kda) ** self.n + 1)), 0, 0,
                  -self.Va * self.n * (A / self.kaa) ** self.n * (D / self.kda) ** self.n / (
                              D * ((A / self.kaa) ** self.n + 1) * ((D / self.kda) ** self.n + 1) ** 2), 0, 0]
            JB = [-self.Vb * self.n * (A / self.kaa) ** (2 * self.n) / (
                        A * ((A / self.kaa) ** self.n + 1) ** 2 * ((E / self.keb) ** self.n + 1)) + self.Vb * self.n * (
                              A / self.kaa) ** self.n / (
                              A * ((A / self.kaa) ** self.n + 1) * ((E / self.keb) ** self.n + 1)), -self.mub, 0, 0,
                  -self.Vb * self.n * (A / self.kaa) ** self.n * (E / self.keb) ** self.n / (
                              E * ((A / self.kaa) ** self.n + 1) * ((E / self.keb) ** self.n + 1) ** 2), 0]
            JC = [-self.Vc * self.n * (A / self.kaa) ** (2 * self.n) / (
                        A * ((A / self.kaa) ** self.n + 1) ** 2 * ((D / self.kda) ** self.n + 1)) + self.Vc * self.n * (
                              A / self.kaa) ** self.n / (
                              A * ((A / self.kaa) ** self.n + 1) * ((D / self.kda) ** self.n + 1)), 0, -self.mulva,
                  -self.Vc * self.n * (A / self.kaa) ** self.n * (D / self.kda) ** self.n / (
                              D * ((A / self.kaa) ** self.n + 1) * ((D / self.kda) ** self.n + 1) ** 2), 0, 0]
            JD = [0, -self.Vd * self.n * (B / self.kbd) ** (2 * self.n) / (
                        B * ((B / self.kbd) ** self.n + 1) ** 2) + self.Vd * self.n * (B / self.kbd) ** self.n / (
                              B * ((B / self.kbd) ** self.n + 1)), 0, -self.mulva, 0, 0]
            JE = [0, 0, -self.Ve * self.n * (C / self.kce) ** self.n * (E / self.kee) ** self.n / (
                        C * ((C / self.kce) ** self.n + 1) ** 2 * ((E / self.kee) ** self.n + 1) * (
                            (F / self.kfe) ** self.n + 1)), 0,
                  -self.mulva - self.Ve * self.n * (E / self.kee) ** (2 * self.n) / (
                              E * ((C / self.kce) ** self.n + 1) * ((E / self.kee) ** self.n + 1) ** 2 * (
                                  (F / self.kfe) ** self.n + 1)) + self.Ve * self.n * (E / self.kee) ** self.n / (
                              E * ((C / self.kce) ** self.n + 1) * ((E / self.kee) ** self.n + 1) * (
                                  (F / self.kfe) ** self.n + 1)),
                  -self.Ve * self.n * (E / self.kee) ** self.n * (F / self.kfe) ** self.n / (
                              F * ((C / self.kce) ** self.n + 1) * ((E / self.kee) ** self.n + 1) * (
                                  (F / self.kfe) ** self.n + 1) ** 2)]
            JF = [0, -self.Vf * self.n * (B / self.kbd) ** (2 * self.n) / (
                        B * ((B / self.kbd) ** self.n + 1) ** 2) + self.Vf * self.n * (B / self.kbd) ** self.n / (
                              B * ((B / self.kbd) ** self.n + 1)), 0, 0, 0, -self.mulva]

            return np.array([JA, JB, JC, JD, JE, JF])

#3954 without positive feedback on nodeA
class circuit2(hill_functions):

    def __init__(self,par_dict,stochasticity=0):
        for key,value in par_dict.items():
            setattr(self,key,value)
        setattr(self, 'stochasticity', stochasticity)


    def dAdt_f(self,species_list, wvn=0):
        A,B,C,D,E,F = species_list
        dadt= self.bA+self.VA*self.noncompetitiveinh(D,self.Kda,self.nda)-self.muASV*A - A*self.DA*wvn**2
        if self.stochasticity ==1:
            dadt+=dadt*normal(0,0.05,1)
        return dadt


    def dBdt_f(self,species_list, wvn=0):
        A,B,C,D,E,F = species_list
        dbdt= self.bB+self.VB*self.noncompetitiveact(A,self.Kab, self.nab)*self.noncompetitiveinh(E,self.Keb, self.neb)-self.muASV*B - B*self.DB*wvn**2
        if self.stochasticity ==1:
            dbdt+=dbdt*normal(0,0.05,1)
        return dbdt


    def dCdt_f(self,species_list):
        A,B,C,D,E,F = species_list
        dcdt= self.bC+self.VC*self.noncompetitiveinh(D,self.Kda, self.nda)-self.muLVA*C
        if self.stochasticity ==1:
            dcdt+=dcdt*normal(0,0.05,1)
        return dcdt

    def dDdt_f(self,species_list):
        A,B,C,D,E,F = species_list
        dddt= self.bD+self.VD*self.noncompetitiveact(B,self.Kbd, self.nbd)-self.muLVA*D
        if self.stochasticity ==1:
            dddt+=dddt*normal(0,0.05,1)
        return dddt

    def dEdt_f(self,species_list):
        A,B,C,D,E,F = species_list
        dedt= self.bE+self.VE*self.noncompetitiveinh(C,self.Kce, self.nce)*self.noncompetitiveinh(F,self.Kfe, self.nfe)*self.noncompetitiveact(E,self.Kee, self.nee)-self.muLVA*E
        if self.stochasticity ==1:
            dedt+=dedt*normal(0,0.05,1)
        return dedt
        
    def dFdt_f(self,species_list):
        A,B,C,D,E,F = species_list
        dfdt= self.bF+self.VF*self.noncompetitiveact(B,self.Kbd, self.nbd)-self.muLVA*F
        if self.stochasticity ==1:
            dfdt+=dfdt*normal(0,0.05,1)
        return dfdt

    function_list = [dAdt_f,dBdt_f,dCdt_f,dDdt_f,dEdt_f,dFdt_f]

    def dudt_growth(self,U, cell_matrix):
        function_list = [self.dAdt_f(U),self.dBdt_f(U),self.dCdt_f(U), self.dDdt_f(U),self.dEdt_f(U),self.dFdt_f(U)]
        dudt = [eq*cell_matrix for eq in function_list]
        return dudt
    @numba.jit(nopython=True)
    def dudt(self,U):
        dudt = [self.dAdt_f(U),self.dBdt_f(U),self.dCdt_f(U), self.dDdt_f(U),self.dEdt_f(U),self.dFdt_f(U)]
        return dudt

    def getJacobian(self,x,wvn):

        A,B,C,D,E,F = x

        # JA = [- self.d_A * wvn ** 2 - self.mua, 0, 0, -(self.Va * self.n * (D / self.kda) ** (self.n - 1)) / (self.kda * ((D / self.kda) ** self.n + 1) ** 2), 0,
        #       0]
        # JB = [(self.Vb * self.n * (A / self.kaa) ** (self.n - 1)) / (
        #             self.kaa * ((A / self.kaa) ** self.n + 1) * ((E / self.keb) ** self.n + 1)) - (
        #               self.Vb * self.n * (A / self.kaa) ** self.n * (A / self.kaa) ** (self.n - 1)) / (
        #               self.kaa * ((A / self.kaa) ** self.n + 1) ** 2 * ((E / self.keb) ** self.n + 1)),
        #       - self.d_B * wvn ** 2 - self.mua, 0, 0,
        #       -(self.Vb * self.n * (A / self.kaa) ** self.n * (E / self.keb) ** (self.n - 1)) / (
        #                   self.keb * ((A / self.kaa) ** self.n + 1) * ((E / self.keb) ** self.n + 1) ** 2), 0]
        # JC = [0, 0, -self.mulva,
        #       -(self.Vc * self.n * (D / self.kda) ** (self.n - 1)) / (self.kda * ((D / self.kda) ** self.n + 1) ** 2), 0,
        #       0]
        # JD = [0, (self.Vd * self.n * (B / self.kbd) ** (self.n - 1)) / (self.kbd * ((B / self.kbd) ** self.n + 1)) - (
        #         self.Vd * self.n * (B / self.kbd) ** self.n * (B / self.kbd) ** (self.n - 1)) / (
        #                   self.kbd * ((B / self.kbd) ** self.n + 1) ** 2), 0, -self.mulva, 0, 0]
        # JE = [0, 0, -(self.Ve * self.n * (C / self.kce) ** (self.n - 1) * (E / self.kee) ** self.n) / (
        #         self.kce * ((C / self.kce) ** self.n + 1) ** 2 * ((E / self.kee) ** self.n + 1) * (
        #             (F / self.kfe) ** self.n + 1)), 0,
        #       (self.Ve * self.n * (E / self.kee) ** (self.n - 1)) / (
        #               self.kee * ((C / self.kce) ** self.n + 1) * ((E / self.kee) ** self.n + 1) * (
        #                   (F / self.kfe) ** self.n + 1)) - self.mulva - (
        #               self.Ve * self.n * (E / self.kee) ** self.n * (E / self.kee) ** (self.n - 1)) / (
        #               self.kee * ((C / self.kce) ** self.n + 1) * ((E / self.kee) ** self.n + 1) ** 2 * (
        #                   (F / self.kfe) ** self.n + 1)),
        #       -(self.Ve * self.n * (E / self.kee) ** self.n * (F / self.kfe) ** (self.n - 1)) / (
        #               self.kfe * ((C / self.kce) ** self.n + 1) * ((E / self.kee) ** self.n + 1) * (
        #                   (F / self.kfe) ** self.n + 1) ** 2)]
        # JF = [0, (self.Vf * self.n * (B / self.kbd) ** (self.n - 1)) / (self.kbd * ((B / self.kbd) ** self.n + 1)) - (
        #         self.Vf * self.n * (B / self.kbd) ** self.n * (B / self.kbd) ** (self.n - 1)) / (
        #                   self.kbd * ((B / self.kbd) ** self.n + 1) ** 2), 0, 0, 0, -self.mulva]


        JA = [-self.DA*wvn**2 - self.muASV, 0, 0, -self.VA*self.nda*(D/self.Kda)**self.nda/(D*((D/self.Kda)**self.nda + 1)**2), 0, 0]
        JB = [-self.VB*self.nab*(A/self.Kab)**(2*self.nab)/(A*((A/self.Kab)**self.nab + 1)**2*((E/self.Keb)**self.neb + 1)) + self.VB*self.nab*(A/self.Kab)**self.nab/(A*((A/self.Kab)**self.nab + 1)*((E/self.Keb)**self.neb + 1)), -self.DB*wvn**2 - self.muASV, 0, 0, -self.VB*self.neb*(A/self.Kab)**self.nab*(E/self.Keb)**self.neb/(E*((A/self.Kab)**self.nab + 1)*((E/self.Keb)**self.neb + 1)**2), 0]
        JC = [0, 0, -self.muLVA, -self.VC*self.nda*(D/self.Kda)**self.nda/(D*((D/self.Kda)**self.nda + 1)**2), 0, 0]
        JD = [0, -self.VD*self.nbd*(B/self.Kbd)**(2*self.nbd)/(B*((B/self.Kbd)**self.nbd + 1)**2) + self.VD*self.nbd*(B/self.Kbd)**self.nbd/(B*((B/self.Kbd)**self.nbd + 1)), 0, -self.muLVA, 0, 0]
        JE = [0, 0, -self.VE*self.nce*(C/self.Kce)**self.nce*(E/self.Kee)**self.nee/(C*((C/self.Kce)**self.nce + 1)**2*((E/self.Kee)**self.nee + 1)*((F/self.Kfe)**self.nfe + 1)), 0, -self.muLVA - self.VE*self.nee*(E/self.Kee)**(2*self.nee)/(E*((C/self.Kce)**self.nce + 1)*((E/self.Kee)**self.nee + 1)**2*((F/self.Kfe)**self.nfe + 1)) + self.VE*self.nee*(E/self.Kee)**self.nee/(E*((C/self.Kce)**self.nce + 1)*((E/self.Kee)**self.nee + 1)*((F/self.Kfe)**self.nfe + 1)), -self.VE*self.nfe*(E/self.Kee)**self.nee*(F/self.Kfe)**self.nfe/(F*((C/self.Kce)**self.nce + 1)*((E/self.Kee)**self.nee + 1)*((F/self.Kfe)**self.nfe + 1)**2)]
        JF = [0, -self.VF*self.nbd*(B/self.Kbd)**(2*self.nbd)/(B*((B/self.Kbd)**self.nbd + 1)**2) + self.VF*self.nbd*(B/self.Kbd)**self.nbd/(B*((B/self.Kbd)**self.nbd + 1)), 0, 0, 0, -self.muLVA]
        return np.array([JA, JB, JC, JD, JE, JF])

#nodeC dele (C and E species removed)
class circuit3(hill_functions):

    def __init__(self,par_dict):
        for key,value in par_dict.items():
            setattr(self,key,value)


    def dAdt_f(self,A,B,D,F,M1,M2):
        dadt= self.ba+self.Va*self.noncompetitiveinh(D,self.kda)-self.muasv*A
        return dadt


    def dBdt_f(self,A,B,D,F,M1,M2):
        dbdt= self.bb+self.Vb*self.noncompetitiveact(M1,self.kM1a)-self.mulva*B
        return dbdt


    def dDdt_f(self,A,B,D,F,M1,M2):
        dddt= self.bd+self.Vd*self.noncompetitiveact(M2,self.kM2d)-self.mulva*D
        return dddt


    def dFdt_f(self,A,B,D,F,M1,M2):
        dfdt= self.bf+self.Vf*self.noncompetitiveact(M2,self.kM2d)-self.mulva*F
        return dfdt

    def dM1dt_f(self,A,B,D,F,M1,M2):
        dM1dt = self.kM1*A - self.muM1*M1
        return dM1dt

    def dM2dt_f(self,A,B,D,F,M1,M2):
        dM2dt = self.kM2*A - self.muM2*M1
        return dM2dt

    def diffusing_dM1dt(self,A,B,D,F,M1,M2,wvn):
        dM1dt = self.kM1*A - self.muM1*M1 - M1*self.d_M1*wvn**2
        return dM1dt

    def diffusing_dM2dt(self,A,B,D,F,M1,M2,wvn):
        dM2dt = self.kM2*B - self.muM2*M2 - M2*self.d_M2*wvn**2
        return dM2dt

    function_list = [dAdt_f,dBdt_f,dDdt_f,dFdt_f,dM1dt_f,dM2dt_f]

    def getJacobian(self, x, wvn):
        A,B,D,F,M1,M2 = x
        JA = [-self.muasv, 0, -self.Va * self.n * (D / self.kda) ** self.n / (D * ((D / self.kda) ** self.n + 1) ** 2),
              0, 0, 0]
        JB = [0, -self.mulva, 0, 0, -self.Vb * self.n * (M1 / self.kM1a) ** (2 * self.n) / (
                    M1 * ((M1 / self.kM1a) ** self.n + 1) ** 2) + self.Vb * self.n * (M1 / self.kM1a) ** self.n / (
                          M1 * ((M1 / self.kM1a) ** self.n + 1)), 0]
        JD = [0, 0, -self.mulva, 0, 0, -self.Vd * self.n * (M2 / self.kM2d) ** (2 * self.n) / (
                    M2 * ((M2 / self.kM2d) ** self.n + 1) ** 2) + self.Vd * self.n * (M2 / self.kM2d) ** self.n / (
                          M2 * ((M2 / self.kM2d) ** self.n + 1))]
        JF = [0, 0, 0, -self.mulva, 0, -self.Vf * self.n * (M2 / self.kM2d) ** (2 * self.n) / (
                    M2 * ((M2 / self.kM2d) ** self.n + 1) ** 2) + self.Vf * self.n * (M2 / self.kM2d) ** self.n / (
                          M2 * ((M2 / self.kM2d) ** self.n + 1))]
        JM1 = [self.kM1, 0, 0, 0, -self.d_M1 * wvn ** 2 - self.muM1, 0]
        JM2 = [0, self.kM2, 0, 0, 0, -self.d_M2 * wvn ** 2 - self.muM2]

        return np.array([JA,JB,JD,JF,JM1,JM2])


class circuit4(hill_functions):

    def __init__(self,par_dict):
        for key,value in par_dict.items():
            setattr(self,key,value)

    def dAdt_f(self, A, B, D, F):
        dadt = self.ba + self.Va * self.noncompetitiveinh(D, self.kda) - self.mua * A
        return dadt

    def dBdt_f(self, A, B, D, F):
        dbdt = self.bb + self.Vb * self.noncompetitiveact(A, self.kaa)  - self.mub * B
        return dbdt


    def dDdt_f(self, A, B, D, F):
        dddt = self.bd + self.Vd * self.noncompetitiveact(B, self.kbd) - self.mulva * D
        return dddt

    def dFdt_f(self, A, B,  D, F):
        dfdt = self.bf + self.Vf * self.noncompetitiveact(B, self.kbd) - self.mulva * F
        return dfdt

    def diffusing_dAdt(self, A, B, D, F,wvn):
        dadt = self.ba + self.Va * self.noncompetitiveinh(D, self.kda) - self.mua * A - A*self.d_A*wvn**2
        return dadt

    def diffusing_dBdt(self, A, B, D, F,wvn):
        dbdt = self.bb + self.Vb * self.noncompetitiveact(A, self.kaa)  - self.mub * B - B*self.d_B*wvn**2
        return dbdt


    function_list = [dAdt_f, dBdt_f, dDdt_f, dFdt_f]

    def getJacobian(self, x, wvn):
        A,B,D,F = x

        JA = [-self.d_A * wvn ** 2 - self.mua, 0,
              -self.Va * self.n * (D / self.kda) ** self.n / (D * ((D / self.kda) ** self.n + 1) ** 2), 0]
        JB = [-self.Vb * self.n * (A / self.kaa) ** (2 * self.n) / (
                    A * ((A / self.kaa) ** self.n + 1) ** 2) + self.Vb * self.n * (A / self.kaa) ** self.n / (
                          A * ((A / self.kaa) ** self.n + 1)), -self.d_B * wvn ** 2 - self.mub, 0, 0]
        JD = [0, -self.Vd * self.n * (B / self.kbd) ** (2 * self.n) / (
                    B * ((B / self.kbd) ** self.n + 1) ** 2) + self.Vd * self.n * (B / self.kbd) ** self.n / (
                          B * ((B / self.kbd) ** self.n + 1)), -self.mulva, 0]
        JF = [0, -self.Vf * self.n * (B / self.kbd) ** (2 * self.n) / (
                    B * ((B / self.kbd) ** self.n + 1) ** 2) + self.Vf * self.n * (B / self.kbd) ** self.n / (
                          B * ((B / self.kbd) ** self.n + 1)), 0, -self.mulva]

        return np.array([JA,JB,JD,JF])


class circuit5(hill_functions):

    def __init__(self,par_dict):
        for key,value in par_dict.items():
            setattr(self,key,value)


    def dAdt_f(self,A,B,D,F):
        dadt= self.ba+self.Va*self.noncompetitiveinh(D,self.kda)-self.mua*A
        return dadt


    def dBdt_f(self,A,B,D,F):
        dbdt= self.bb+self.Vb*self.noncompetitiveact(A,self.kaa)-self.mua*B
        return dbdt

    def dDdt_f(self,A,B,D,F):
        dddt= self.bd+self.Vd*self.noncompetitiveact(B,self.kbd)-self.mulva*D
        return dddt

    def dFdt_f(self,A,B,D,F):
        dfdt= self.bf+self.Vf*self.noncompetitiveact(B,self.kbd)-self.mulva*F
        return dfdt
    function_list = [dAdt_f,dBdt_f,dDdt_f,dFdt_f]


    def getJacobian(self,x,wvn):

        A,B,D,F = x

        JA = [- self.d_A * wvn ** 2 - self.mua, 0, 0,
              -(self.Va * self.n * (D / self.kda) ** (self.n - 1)) / (self.kda * ((D / self.kda) ** self.n + 1) ** 2), 0,
              0]
        JB = [(self.Vb * self.n * (A / self.kaa) ** (self.n - 1)) / (
                    self.kaa * ((A / self.kaa) ** self.n + 1) * ((E / self.keb) ** self.n + 1)) - (
                      self.Vb * self.n * (A / self.kaa) ** self.n * (A / self.kaa) ** (self.n - 1)) / (
                      self.kaa * ((A / self.kaa) ** self.n + 1) ** 2 * ((E / self.keb) ** self.n + 1)),
              - self.d_B * wvn ** 2 - self.mua, 0, 0,
              -(self.Vb * self.n * (A / self.kaa) ** self.n * (E / self.keb) ** (self.n - 1)) / (
                          self.keb * ((A / self.kaa) ** self.n + 1) * ((E / self.keb) ** self.n + 1) ** 2), 0]
        JC = [0, 0, -self.mulva,
              -(self.Vc * self.n * (D / self.kda) ** (self.n - 1)) / (self.kda * ((D / self.kda) ** self.n + 1) ** 2), 0,
              0]
        JD = [0, (self.Vd * self.n * (B / self.kbd) ** (self.n - 1)) / (self.kbd * ((B / self.kbd) ** self.n + 1)) - (
                self.Vd * self.n * (B / self.kbd) ** self.n * (B / self.kbd) ** (self.n - 1)) / (
                          self.kbd * ((B / self.kbd) ** self.n + 1) ** 2), 0, -self.mulva, 0, 0]
        JE = [0, 0, -(self.Ve * self.n * (C / self.kce) ** (self.n - 1) * (E / self.kee) ** self.n) / (
                self.kce * ((C / self.kce) ** self.n + 1) ** 2 * ((E / self.kee) ** self.n + 1) * (
                    (F / self.kfe) ** self.n + 1)), 0,
              (self.Ve * self.n * (E / self.kee) ** (self.n - 1)) / (
                      self.kee * ((C / self.kce) ** self.n + 1) * ((E / self.kee) ** self.n + 1) * (
                          (F / self.kfe) ** self.n + 1)) - self.mulva - (
                      self.Ve * self.n * (E / self.kee) ** self.n * (E / self.kee) ** (self.n - 1)) / (
                      self.kee * ((C / self.kce) ** self.n + 1) * ((E / self.kee) ** self.n + 1) ** 2 * (
                          (F / self.kfe) ** self.n + 1)),
              -(self.Ve * self.n * (E / self.kee) ** self.n * (F / self.kfe) ** (self.n - 1)) / (
                      self.kfe * ((C / self.kce) ** self.n + 1) * ((E / self.kee) ** self.n + 1) * (
                          (F / self.kfe) ** self.n + 1) ** 2)]
        JF = [0, (self.Vf * self.n * (B / self.kbd) ** (self.n - 1)) / (self.kbd * ((B / self.kbd) ** self.n + 1)) - (
                self.Vf * self.n * (B / self.kbd) ** self.n * (B / self.kbd) ** (self.n - 1)) / (
                          self.kbd * ((B / self.kbd) ** self.n + 1) ** 2), 0, 0, 0, -self.mulva]

        return np.array([JA, JB, JD, JF])

#PORD MODEL ALEXANDER BSC
class circuit6(hill_functions):
# df = {'k1':0.05, 'k2': 0.01, 'k3':2, 'n1':3, 'n2':3, 'n3':1, 'beta':2.8,'Va': alpha, 'Vb':alpha, 'mua':1, 'mub':1, 'd_A':5e-05, 'd_B':0.0025}

    def __init__(self,par_dict):
        # original_par_dict = {'k1':0.029, 'k2': 0.008, 'k3':1.658, 'n':2.411, 'beta':1,'Va': 1, 'Vb':1, 'mua':1, 'mub':1, 'd_A':0, 'd_B':1}
        # for key,value in original_par_dict.items():
        #     setattr(self,key,value)

        for key,value in par_dict.items():
            setattr(self,key,value)

    def dAdt_f(self, species_list):
        A, B = species_list
        dadt =  self.Va * ((self.beta + (A/self.k1)**self.n1)/(1+(A/self.k1)**self.n1 + (B/self.k2)**self.n2 )) - self.mua * A
        return dadt

    def dBdt_f(self, species_list):
        A, B = species_list
        dbdt = self.Vb * (((A/self.k3)**self.n3)/(1+(A/self.k3)**self.n3)) - self.mub * B
        return dbdt

    def diffusing_dAdt(self, A, B,wvn):
        dadt =  self.Va * ((self.beta + (A/self.k1)**self.n1)/(1+(A/self.k1)**self.n2 + (B/self.k2)**self.n2 )) - self.mua * A - A*self.d_A*wvn**2
        return dadt

    def diffusing_dBdt(self, A, B,wvn):
        dbdt = self.Vb * (((A/self.k3)**self.n3)/(1+(A/self.k3)**self.n3)) - self.mub * B - B*self.d_B*wvn**2
        return dbdt


    function_list = [dAdt_f, dBdt_f]

    def getJacobian(self, x, wvn):
        A,B = x
        JA = [-self.d_A*wvn**2 - self.mua + self.Va*self.n1*(A/self.k1)**self.n1/(A*((A/self.k1)**self.n2 + (B/self.k2)**self.n2 + 1)) - self.Va*self.n2*(A/self.k1)**self.n2*(self.beta + (A/self.k1)**self.n1)/(A*((A/self.k1)**self.n2 + (B/self.k2)**self.n2 + 1)**2),
        -self.Va*self.n2*(B/self.k2)**self.n2*(self.beta + (A/self.k1)**self.n1)/(B*((A/self.k1)**self.n2 + (B/self.k2)**self.n2 + 1)**2)]
        JB = [-self.Vb*self.n3*(A/self.k3)**(2*self.n3)/(A*((A/self.k3)**self.n3 + 1)**2) + self.Vb*self.n3*(A/self.k3)**self.n3/(A*((A/self.k3)**self.n3 + 1)),
        -self.d_B*wvn**2 - self.mub]

        return np.array([JA,JB])

#PHILIP MAINI STOS SPATIOTEMPORAL OSCILLATING SOLUTION
class circuit7(hill_functions):
    # df = {'alpha':0.95, 'beta': -0.91, 'D':0.45, 'r3':0,  'r2':0}

    def __init__(self,par_dict):
        for key,value in par_dict.items():
            setattr(self,key,value)
        self.d_A = self.D
        self.d_B=1


    def dAdt_f(self,species_list):
        A, B = species_list
        dadt= self.alpha*A + B - self.alpha*self.r3*A*B**2 - self.r2*A*B
        return dadt


    def dBdt_f(self,species_list):
        A, B = species_list
        dbdt = -self.alpha*A + self.beta*B + self.alpha*self.r3*A*B**2 + self.r2*A*B
        return dbdt


    def diffusing_dAdt(self,A,B,wvn):
        dadt= self.alpha*A + B - self.alpha*self.r3*A*B**2 - self.r2*A*B - A*self.gamma*self.D*wvn**2
        return dadt


    def diffusing_dBdt(self,A,B,wvn):
        dbdt = -self.alpha*A + self.beta*B + self.alpha*self.r3*A*B**2 + self.r2*A*B - B*self.gamma*wvn**2
        return dbdt



    function_list = [dAdt_f,dBdt_f]


    def getJacobian(self,x,wvn):

        A,B = x
        JA = [-B**2*self.alpha*self.r3 - B*self.r2 - self.D*self.gamma*wvn**2 + self.alpha, -2*A*B*self.alpha*self.r3 - A*self.r2 + 1]
        JB = [B**2*self.alpha*self.r3 + B*self.r2 - self.alpha, 2*A*B*self.alpha*self.r3 + A*self.r2 + self.beta - self.gamma*wvn**2]

        return np.array([JA, JB])

#circuit2 (3954): nodeA dele
class circuit8(hill_functions):

    def __init__(self,par_dict,stochasticity=1):
        for key,value in par_dict.items():
            setattr(self,key,value)
        setattr(self, 'stochasticity', stochasticity)

    def dBdt_f(self,species_list):
        B,E,F = species_list
        dbdt= self.bb+self.Vb**self.noncompetitiveinh(E,self.keb)-self.mua*B
        dbdt+=dbdt*normal(0,0.05,1)*self.stochasticity
        return dbdt
    def diffusing_dBdt_f(self,species_list,wvn):
        B,E,F = species_list
        dbdt= self.bb+self.Vb**self.noncompetitiveinh(E,self.keb)-self.mua*B - B*self.d_B*wvn**2
        return dbdt

    def dEdt_f(self,species_list):
        B,E,F = species_list
        dedt= self.be+self.Ve**self.noncompetitiveinh(F,self.kfe)*self.noncompetitiveact(E,self.kee)-self.mulva*E
        dedt+=dedt*normal(0,0.05,1)*self.stochasticity
        return dedt
    def dFdt_f(self,species_list):
        B,E,F = species_list
        dfdt= self.bf+self.Vf*self.noncompetitiveact(B,self.kbd)-self.mulva*F
        dfdt+=dfdt*normal(0,0.05,1)*self.stochasticity
        return dfdt
    def dudt(self,U, cell_matrix):
        function_list = [self.dBdt_f(U),self.dEdt_f(U),self.dFdt_f(U)]
        dudt = [eq*cell_matrix for eq in function_list]
        return dudt


    def getJacobian(self,x,wvn):

        B,D,E,F = x
        JB = [-self.d_B*wvn**2 - self.mua,  -self.Vb**(1/((E/self.keb)**self.n + 1))*self.n*(E/self.keb)**self.n*log(self.Vb)/(E*((E/self.keb)**self.n + 1)**2), 0]
        JE = [0,  -self.mulva - self.Ve**(1/((F/self.kfe)**self.n + 1))*self.n*(E/self.kee)**(2*self.n)/(E*((E/self.kee)**self.n + 1)**2) + self.Ve**(1/((F/self.kfe)**self.n + 1))*self.n*(E/self.kee)**self.n/(E*((E/self.kee)**self.n + 1)), -self.Ve**(1/((F/self.kfe)**self.n + 1))*self.n*(E/self.kee)**self.n*(F/self.kfe)**self.n*log(self.Ve)/(F*((E/self.kee)**self.n + 1)*((F/self.kfe)**self.n + 1)**2)]
        JF = [-self.Vf*self.n*(B/self.kbd)**(2*self.n)/(B*((B/self.kbd)**self.n + 1)**2) + self.Vf*self.n*(B/self.kbd)**self.n/(B*((B/self.kbd)**self.n + 1)), 0, -self.mulva]


        return np.array([JB, JE, JF])

#circuit2 (3954): nodeB dele (top cassete: cinI dele)
class circuit9(hill_functions):

    def __init__(self,par_dict):
        for key,value in par_dict.items():
            setattr(self,key,value)


    def dAdt_f(self,species_list):
        A,C,D,E,F = species_list
        dadt= self.ba+self.Va*self.noncompetitiveinh(D,self.kda)-self.mua*A
        return dadt

    def diffusing_dAdt_f(self,species_list,wvn):
        A,C,D,E,F = species_list
        dadt= self.ba+self.Va*self.noncompetitiveinh(D,self.kda)-self.mua*A - A*self.d_B*wvn**2
        return dadt

    def dCdt_f(self,species_list):
        A,C,D,E,F = species_list
        dcdt= self.bc+self.Vc*self.noncompetitiveinh(D,self.kda)-self.mulva*C
        return dcdt

    def dDdt_f(self,species_list):
        A,C,D,E,F = species_list
        dddt= self.bd-self.mulva*D
        return dddt

    def dEdt_f(self,species_list):
        A,C,D,E,F = species_list
        dedt= self.be+self.Ve*self.noncompetitiveinh(C,self.kce)*self.noncompetitiveinh(F,self.kfe)*self.noncompetitiveact(E,self.kee)-self.mulva*E
        return dedt

    def dFdt_f(self,species_list):
        A,C,D,E,F = species_list
        dfdt= self.bf-self.mulva*F
        return dfdt
    function_list = [dAdt_f,dCdt_f,dDdt_f,dEdt_f,dFdt_f]


    def getJacobian(self,x,wvn):

        A,C,D,E,F = x

        JA = [-self.d_B*wvn**2 - self.mua, 0, -self.Va*self.n*(D/self.kda)**self.n/(D*((D/self.kda)**self.n + 1)**2), 0, 0]
        JC = [0, -self.mulva, -self.Vc*self.n*(D/self.kda)**self.n/(D*((D/self.kda)**self.n + 1)**2), 0, 0]
        JD = [0, 0, -self.mulva, 0, 0]
        JE = [0, -self.Ve*self.n*(C/self.kce)**self.n*(E/self.kee)**self.n/(C*((C/self.kce)**self.n + 1)**2*((E/self.kee)**self.n + 1)*((F/self.kfe)**self.n + 1)), 0, -self.mulva - self.Ve*self.n*(E/self.kee)**(2*self.n)/(E*((C/self.kce)**self.n + 1)*((E/self.kee)**self.n + 1)**2*((F/self.kfe)**self.n + 1)) + self.Ve*self.n*(E/self.kee)**self.n/(E*((C/self.kce)**self.n + 1)*((E/self.kee)**self.n + 1)*((F/self.kfe)**self.n + 1)), -self.Ve*self.n*(E/self.kee)**self.n*(F/self.kfe)**self.n/(F*((C/self.kce)**self.n + 1)*((E/self.kee)**self.n + 1)*((F/self.kfe)**self.n + 1)**2)]
        JF = [0, 0, 0, 0, -self.mulva]

        return np.array([JA, JC, JD, JE, JF])


#circuit2 (3954): nodeA dele + receptor binding and diffuser
class circuit10(hill_functions):#correct **

    def __init__(self,par_dict):
        for key,value in par_dict.items():
            setattr(self,key,value)
        self.muB = self.muLVA
        self.muD = self.muLVA
        self.muE = self.muLVA
        self.muF = self.muLVA
        self.khf = self.khd

    def function_list(self,species_list,wvn):
        B,H,D,E,F = species_list


        dBdt_f= self.bb+self.Vb**self.noncompetitiveinh(E,self.keb)-self.muB*B
        dHdt_f= self.alphaH*B - self.muH*H - H*self.d_H*wvn**2
        HR=self.R/(1+(self.kD/H)) #Quasisteadystate
        dDdt_f= self.bd+self.Vd*self.noncompetitiveact(HR,self.khd)-self.muD*D
        dEdt_f= self.be+self.Ve**self.noncompetitiveinh(F,self.kfe)*self.noncompetitiveact(E,self.kee)-self.muE*E
        dFdt_f= self.bf+self.Vf*self.noncompetitiveact(HR,self.khf)-self.muF*F

        return [dBdt_f,dHdt_f,dDdt_f,dEdt_f,dFdt_f]



    def getJacobian(self,x,wvn):

        B,H,D,E,F = x
        JB = [-self.muB, 0, 0, -self.Vb**(1/((E/self.keb)**self.n + 1))*self.n*(E/self.keb)**self.n*log(self.Vb)/(E*((E/self.keb)**self.n + 1)**2), 0]
        JH = [self.alphaH, -self.d_H*wvn**2 - self.muH, 0, 0, 0]
        JD = [0, -self.Vd*self.n*(H*self.R*self.kon/(self.khd*self.koff))**(2*self.n)/(H*((H*self.R*self.kon/(self.khd*self.koff))**self.n + 1)**2) + self.Vd*self.n*(H*self.R*self.kon/(self.khd*self.koff))**self.n/(H*((H*self.R*self.kon/(self.khd*self.koff))**self.n + 1)), -self.muD, 0, 0]
        JE = [0, 0, 0, -self.muE - self.Ve**(1/((F/self.kfe)**self.n + 1))*self.n*(E/self.kee)**(2*self.n)/(E*((E/self.kee)**self.n + 1)**2) + self.Ve**(1/((F/self.kfe)**self.n + 1))*self.n*(E/self.kee)**self.n/(E*((E/self.kee)**self.n + 1)), -self.Ve**(1/((F/self.kfe)**self.n + 1))*self.n*(E/self.kee)**self.n*(F/self.kfe)**self.n*log(self.Ve)/(F*((E/self.kee)**self.n + 1)*((F/self.kfe)**self.n + 1)**2)]
        JF = [0, -self.Vf*self.n*(H*self.R*self.kon/(self.khf*self.koff))**(2*self.n)/(H*((H*self.R*self.kon/(self.khf*self.koff))**self.n + 1)**2) + self.Vf*self.n*(H*self.R*self.kon/(self.khf*self.koff))**self.n/(H*((H*self.R*self.kon/(self.khf*self.koff))**self.n + 1)), 0, 0, -self.muF]
        return np.array([JB, JH, JD, JE, JF])


#circuit2 (3954): nodeA dele + receptor binding and diffuser
class circuit10(hill_functions):

    def __init__(self,par_dict):
        for key,value in par_dict.items():
            setattr(self,key,value)
        self.muB = self.muLVA
        self.muD = self.muLVA
        self.muE = self.muLVA
        self.muF = self.muLVA
        self.khf = self.khd

    def function_list(self,species_list,wvn):
        B,H,D,E,F = species_list


        dBdt_f= self.bb+self.Vb**self.noncompetitiveinh(E,self.keb)-self.muB*B
        dHdt_f= self.alphaH*B - self.muH*H - H*self.d_H*wvn**2
        HR=self.R/(1+(self.kD/H)) #Quasisteadystate
        dDdt_f= self.bd+self.Vd*self.noncompetitiveact(HR,self.khd)-self.muD*D
        dEdt_f= self.be+self.Ve**self.noncompetitiveinh(F,self.kfe)*self.noncompetitiveact(E,self.kee)-self.muE*E
        dFdt_f= self.bf+self.Vf*self.noncompetitiveact(HR,self.khf)-self.muF*F

        return [dBdt_f,dHdt_f,dDdt_f,dEdt_f,dFdt_f]



    def getJacobian(self,x,wvn):

        B,H,D,E,F = x
        JB = [-self.muLVA, 0, 0, -self.Vb**(1/((E/self.keb)**self.n + 1))*self.n*(E/self.keb)**self.n*log(self.Vb)/(E*((E/self.keb)**self.n + 1)**2), 0]
        JH = [self.alphaH, -self.d_H*wvn**2 - self.muH, 0, 0, 0]
        JD = [0, -self.Vd*self.kD*self.n*(self.R/(self.khd*(1 + self.kD/H)))**(2*self.n)/(H**2*(1 + self.kD/H)*((self.R/(self.khd*(1 + self.kD/H)))**self.n + 1)**2) + self.Vd*self.kD*self.n*(self.R/(self.khd*(1 + self.kD/H)))**self.n/(H**2*(1 + self.kD/H)*((self.R/(self.khd*(1 + self.kD/H)))**self.n + 1)), -self.muLVA, 0, 0]
        JE = [0, 0, 0, -self.muLVA - self.Ve**(1/((F/self.kfe)**self.n + 1))*self.n*(E/self.kee)**(2*self.n)/(E*((E/self.kee)**self.n + 1)**2) + self.Ve**(1/((F/self.kfe)**self.n + 1))*self.n*(E/self.kee)**self.n/(E*((E/self.kee)**self.n + 1)), -self.Ve**(1/((F/self.kfe)**self.n + 1))*self.n*(E/self.kee)**self.n*(F/self.kfe)**self.n*log(self.Ve)/(F*((E/self.kee)**self.n + 1)*((F/self.kfe)**self.n + 1)**2)]
        JF = [0, -self.Vf*self.kD*self.n*(self.R/(self.khd*(1 + self.kD/H)))**(2*self.n)/(H**2*(1 + self.kD/H)*((self.R/(self.khd*(1 + self.kD/H)))**self.n + 1)**2) + self.Vf*self.kD*self.n*(self.R/(self.khd*(1 + self.kD/H)))**self.n/(H**2*(1 + self.kD/H)*((self.R/(self.khd*(1 + self.kD/H)))**self.n + 1)), 0, 0, -self.muLVA]
        return np.array([JB, JH, JD, JE, JF])


#circuit2 (3954): nodeA dele + receptor binding and diffuser
class circuit11(hill_functions):

    def __init__(self,par_dict):
        for key,value in par_dict.items():
            setattr(self,key,value)
        self.muB = self.muLVA
        self.muE = self.muLVA
        self.muF = self.muLVA
        self.khmf = self.khmd

    def dudt(self,species_list,cell_matrix,wvn=0):
        mB,H,mE,mF,B,E,F = species_list


        dmBdt_f= self.bmB+self.VmB*self.noncompetitiveinh(E,self.kemb,self.nE)-self.muRNA*mB
        dHdt_f= (self.alphaH*B - self.muH*H) - H*self.d_H*wvn**2
        HR=self.R/(1+(self.kD/H)) #Quasisteadystate
        dmEdt_f= self.bmE+self.VmE*self.noncompetitiveinh(F,self.kfme,self.nF)*self.noncompetitiveact(E,self.keme,self.nE)-self.muRNA*mE
        dmFdt_f= self.bmF+self.VmF*self.noncompetitiveact(HR,self.khmf,self.nH)-self.muRNA*mF
        dBdt_f = self.aB*mB - self.muB*B
        dEdt_f = self.aE*mE - self.muE*E
        dFdt_f = self.aF*mF - self.muF*F
        function_list = [dmBdt_f,dHdt_f,dmEdt_f,dmFdt_f,dBdt_f,dEdt_f,dFdt_f]

        dudt = [eq*cell_matrix for eq in function_list]
        return dudt
        # def dudt(self,U, cell_matrix):
        #     function_list = [self.dBdt_f(U),self.dEdt_f(U),self.dFdt_f(U)]
        #     dudt = [eq*cell_matrix for eq in function_list]
        #     return dudt

        # return [dmBdt_f,dHdt_f,dmDdt_f,dmEdt_f,dmFdt_f,dBdt_f,dDdt_f,dEdt_f,dFdt_f]



    def getJacobian(self,x,wvn):

        mB,H,mD,mE,mF,B,D,E,F = x


        JmB = [-self.muRNA, 0, 0, 0, 0, 0, 0, -self.VmB*self.n*(E/self.kemb)**self.n/(E*((E/self.kemb)**self.n + 1)**2), 0]
        JH = [0, -self.d_H*wvn**2 - self.muH, 0, 0, 0, self.alphaH, 0, 0, 0]
        JmD = [0, -self.VmD*self.kD*self.n*(self.R/(self.khmd*(1 + self.kD/H)))**(2*self.n)/(H**2*(1 + self.kD/H)*((self.R/(self.khmd*(1 + self.kD/H)))**self.n + 1)**2) + self.VmD*self.kD*self.n*(self.R/(self.khmd*(1 + self.kD/H)))**self.n/(H**2*(1 + self.kD/H)*((self.R/(self.khmd*(1 + self.kD/H)))**self.n + 1)), -self.muRNA, 0, 0, 0, 0, 0, 0]
        JmE = [0, 0, 0, -self.muRNA, 0, 0, 0, -self.VmE*self.n*(E/self.keme)**(2*self.n)/(E*((E/self.keme)**self.n + 1)**2*((F/self.kfme)**self.n + 1)) + self.VmE*self.n*(E/self.keme)**self.n/(E*((E/self.keme)**self.n + 1)*((F/self.kfme)**self.n + 1)), -self.VmE*self.n*(E/self.keme)**self.n*(F/self.kfme)**self.n/(F*((E/self.keme)**self.n + 1)*((F/self.kfme)**self.n + 1)**2)]
        JmF = [0, -self.VmF*self.kD*self.n*(self.R/(self.khmd*(1 + self.kD/H)))**(2*self.n)/(H**2*(1 + self.kD/H)*((self.R/(self.khmd*(1 + self.kD/H)))**self.n + 1)**2) + self.VmF*self.kD*self.n*(self.R/(self.khmd*(1 + self.kD/H)))**self.n/(H**2*(1 + self.kD/H)*((self.R/(self.khmd*(1 + self.kD/H)))**self.n + 1)), 0, 0, -self.muRNA, 0, 0, 0, 0]
        JB = [self.aB, 0, 0, 0, 0, -self.muLVA, 0, 0, 0]
        JD = [0, 0, self.aD, 0, 0, 0, -self.muLVA, 0, 0]
        JE = [0, 0, 0, self.aE, 0, 0, 0, -self.muLVA, 0]
        JF = [0, 0, 0, 0, self.aF, 0, 0, 0, -self.muLVA]

        return np.array([JmB,JH,JmD,JmE,JmF,JB,JD,JE,JF])

#diffusers, proteins and regulators modelled
class circuit12(hill_functions):

    def __init__(self,par_dict,stochasticity=0):
        for key,value in par_dict.items():
            setattr(self,key,value)
        setattr(self, 'stochasticity', stochasticity)


    def dUdt_f(self,X,wvn=0):
        U,V,A,B,C,D,E,F,aTc = X
        dUdt= self.k1*A - self.muU*U - U*self.DU*wvn**2
        if self.stochasticity ==1:
            dUdt+=dUdt*normal(0,0.05,1)
        return dUdt

    def dVdt_f(self,X,wvn=0):
        U,V,A,B,C,D,E,F,aTc = X
        dVdt= self.k2*A - self.muV*V - V*self.DV*wvn**2
        if self.stochasticity ==1:
            dVdt+=dVdt*normal(0,0.05,1)
        return dVdt


    def dAdt_f(self,X):
        U,V,A,B,C,D,E,F,aTc = X
        iptg_regulated_K = self.Kda * (1+(self.iptg/self.Kiptg)**self.niptg)
        dAdt= self.bA+self.VA*self.noncompetitiveinh(D,iptg_regulated_K,self.nda)-self.muASV*A
        if self.stochasticity ==1:
            dAdt+=dAdt*normal(0,0.05,1)
        return dAdt


    def dBdt_f(self,X):
        U,V,A,B,C,D,E,F,aTc = X
        dBdt= self.bB+self.VB*self.noncompetitiveact(U,self.Kab,self.nab)*self.noncompetitiveinh(E,self.Keb,self.nab)-self.muASV*B
        if self.stochasticity ==1:
            dBdt+=dBdt*normal(0,0.05,1)
        return dBdt


    def dCdt_f(self,X):
        U,V,A,B,C,D,E,F,aTc = X
        iptg_regulated_K = self.Kda * (1+(self.iptg/self.Kiptg)**self.niptg)
        dCdt= self.bC+self.VC*self.noncompetitiveinh(D,iptg_regulated_K, self.nda)-self.muLVA*C
        if self.stochasticity ==1:
            dCdt+=dCdt*normal(0,0.05,1)
        return dCdt

    def dDdt_f(self,X):
        U,V,A,B,C,D,E,F,aTc = X
        dDdt= self.bD+self.VD*self.noncompetitiveact(V,self.Kbd,self.nbd)-self.muLVA*D
        if self.stochasticity ==1:
            dDdt+=dDdt*normal(0,0.05,1)
        return dDdt

    def dEdt_f(self,X):
        U,V,A,B,C,D,E,F,aTc = X
        aTc_regulated_K = self.Kce * (1+(aTc/self.KaTc)**self.naTc)
        dEdt= self.bE+self.VE*self.noncompetitiveinh(C,aTc_regulated_K,self.nce)*self.noncompetitiveinh(F,self.Kfe, self.nfe)*self.noncompetitiveact(E,self.Kee,self.nee)-self.muLVA*E
        if self.stochasticity ==1:
            dEdt+=dEdt*normal(0,0.05,1)
        return dEdt
        
    def dFdt_f(self,X):
        U,V,A,B,C,D,E,F,aTc = X
        dFdt= self.bF+self.VF*self.noncompetitiveact(V,self.Kbd,self.nbd)-self.muLVA*F
        if self.stochasticity ==1:
            dFdt+=dFdt*normal(0,0.05,1)
        return dFdt
    
    def daTcdt_f(self,X):
        U,V,A,B,C,D,E,F,aTc = X
        daTcdt= -self.muaTc*aTc 
        if self.stochasticity ==1:
            daTcdt+=daTcdt*normal(0,0.05,1)
        return daTcdt

    # def function_list(self,X,wvn):
        # return [self.dUdt_f,self.dVdt_f,self.dAdt_f,self.dBdt_f,self.dCdt_f,self.dDdt_f,self.dEdt_f,self.dFdt_f,self.daTc_f]
    function_list = [dUdt_f,dVdt_f,dAdt_f,dBdt_f,dCdt_f,dDdt_f,dEdt_f,dFdt_f,daTcdt_f]

    def dudt_growth(self,X, cell_matrix):
        function_list = [self.dUdt_f(X),self.dVdt_f(X),self.dAdt_f(X),self.dBdt_f(X),self.dCdt_f(X), self.dDdt_f(X),self.dEdt_f(X),self.dFdt_f(X), self.daTcdt_f(X)]
        dudt = [eq*cell_matrix for eq in function_list]
        return dudt

    def dudt(self,X):
        dudt = [self.dUdt_f(X),self.dVdt_f(X),self.dAdt_f(X),self.dBdt_f(X),self.dCdt_f(X), self.dDdt_f(X),self.dEdt_f(X),self.dFdt_f(X), self.daTcdt_f(X)]
        return dudt


    # def dudt(self,X,growth=False,cell_matrix=None):
    #     dXdt = [self.dUdt_f(X),self.dVdt_f(X),self.dAdt_f(X),self.dBdt_f(X),self.dCdt_f(X), self.dDdt_f(X),self.dEdt_f(X),self.dFdt_f(X), self.daTcdt_f(X)]
    #     if growth==False:
    #         return dXdt
    #     if growth==True:
    #         dXdt_growth = [eq*cell_matrix for eq in dXdt]


    # def dudt_growth(self,U, cell_matrix):
    #     function_list = [self.dAdt_f(U),self.dBdt_f(U),self.dCdt_f(U), self.dDdt_f(U),self.dEdt_f(U),self.dFdt_f(U)]
    #     dudt = [eq*cell_matrix for eq in function_list]
    #     return dudt

    def getJacobian(self,X,wvn):

        U,V,A,B,C,D,E,F,aTc = X
        JU = [-self.DU*wvn**2 - self.muU, 0, self.k1, 0, 0, 0, 0, 0, 0]
        JV = [0, -self.DV*wvn**2 - self.muV, self.k2, 0, 0, 0, 0, 0, 0]
        JA = [0, 0, -self.muASV, 0, 0, -self.VA*self.nda*(D/(self.Kda*((self.iptg/self.Kiptg)**self.niptg + 1)))**self.nda/(D*((D/(self.Kda*((self.iptg/self.Kiptg)**self.niptg + 1)))**self.nda + 1)**2), 0, 0, 0]
        JB = [self.VB*self.nab*(self.Kab/U)**self.nab/(U*((E/self.Keb)**self.nab + 1)*((self.Kab/U)**self.nab + 1)**2), 0, 0, -self.muASV, 0, 0, -self.VB*self.nab*(E/self.Keb)**self.nab/(E*((E/self.Keb)**self.nab + 1)**2*((self.Kab/U)**self.nab + 1)), 0, 0]
        JC = [0, 0, 0, 0, -self.muLVA, -self.VC*self.nda*(D/(self.Kda*((self.iptg/self.Kiptg)**self.niptg + 1)))**self.nda/(D*((D/(self.Kda*((self.iptg/self.Kiptg)**self.niptg + 1)))**self.nda + 1)**2), 0, 0, 0]
        JD = [0, self.VD*self.nbd*(self.Kbd/V)**self.nbd/(V*((self.Kbd/V)**self.nbd + 1)**2), 0, 0, 0, -self.muLVA, 0, 0, 0]
        JE = [0, 0, 0, 0, -self.VE*self.nce*(C/(self.Kce*((aTc/self.KaTc)**self.naTc + 1)))**self.nce/(C*((self.Kee/E)**self.nee + 1)*((F/self.Kfe)**self.nfe + 1)*((C/(self.Kce*((aTc/self.KaTc)**self.naTc + 1)))**self.nce + 1)**2), 0, -self.muLVA + self.VE*self.nee*(self.Kee/E)**self.nee/(E*((self.Kee/E)**self.nee + 1)**2*((F/self.Kfe)**self.nfe + 1)*((C/(self.Kce*((aTc/self.KaTc)**self.naTc + 1)))**self.nce + 1)), -self.VE*self.nfe*(F/self.Kfe)**self.nfe/(F*((self.Kee/E)**self.nee + 1)*((F/self.Kfe)**self.nfe + 1)**2*((C/(self.Kce*((aTc/self.KaTc)**self.naTc + 1)))**self.nce + 1)), self.VE*self.naTc*self.nce*(aTc/self.KaTc)**self.naTc*(C/(self.Kce*((aTc/self.KaTc)**self.naTc + 1)))**self.nce/(aTc*((self.Kee/E)**self.nee + 1)*((F/self.Kfe)**self.nfe + 1)*((aTc/self.KaTc)**self.naTc + 1)*((C/(self.Kce*((aTc/self.KaTc)**self.naTc + 1)))**self.nce + 1)**2)]
        JF = [0, self.VF*self.nbd*(self.Kbd/V)**self.nbd/(V*((self.Kbd/V)**self.nbd + 1)**2), 0, 0, 0, 0, 0, -self.muLVA, 0]
        JaTc = [0, 0, 0, 0, 0, 0, 0, 0, -self.muaTc]
        return np.array([JU,JV,JA,JB,JC,JD,JE,JF,JaTc ])



#diffusers, proteins modelled
class circuit13(hill_functions):

    def __init__(self,par_dict,stochasticity=0):
        for key,value in par_dict.items():
            setattr(self,key,value)
        setattr(self, 'stochasticity', stochasticity)

    def dUdt_f(self,X,wvn=0):
        U,V,A,B,C,D,E,F = X
        dUdt= self.k1*A - self.muU*U - U*self.DU*wvn**2
        if self.stochasticity ==1:
            dUdt+=dUdt*normal(0,0.05,1)
        return dUdt

    def dVdt_f(self,X,wvn=0):
        U,V,A,B,C,D,E,F = X
        dVdt= self.k2*A - self.muV*V - V*self.DV*wvn**2
        if self.stochasticity ==1:
            dVdt+=dVdt*normal(0,0.05,1)
        return dVdt


    def dAdt_f(self,X):
        U,V,A,B,C,D,E,F = X
        dAdt= self.bA+self.VA*self.noncompetitiveinh(D,self.Kda,self.nda)-self.muASV*A
        if self.stochasticity ==1:
            dAdt+=dAdt*normal(0,0.05,1)
        return dAdt


    def dBdt_f(self,X):
        U,V,A,B,C,D,E,F = X
        dBdt= self.bB+self.VB*self.noncompetitiveact(U,self.Kab,self.nab)*self.noncompetitiveinh(E,self.Keb,self.nab)-self.muASV*B
        if self.stochasticity ==1:
            dBdt+=dBdt*normal(0,0.05,1)
        return dBdt


    def dCdt_f(self,X):
        U,V,A,B,C,D,E,F = X
        dCdt= self.bC+self.VC*self.noncompetitiveinh(D,self.Kda, self.nda)-self.muLVA*C
        if self.stochasticity ==1:
            dCdt+=dCdt*normal(0,0.05,1)
        return dCdt

    def dDdt_f(self,X):
        U,V,A,B,C,D,E,F = X
        dDdt= self.bD+self.VD*self.noncompetitiveact(V,self.Kbd,self.nbd)-self.muLVA*D
        if self.stochasticity ==1:
            dDdt+=dDdt*normal(0,0.05,1)
        return dDdt

    def dEdt_f(self,X):
        U,V,A,B,C,D,E,F = X
        dEdt= self.bE+self.VE*self.noncompetitiveinh(C,self.Kce,self.nce)*self.noncompetitiveinh(F,self.Kfe, self.nfe)*self.noncompetitiveact(E,self.Kee,self.nee)-self.muLVA*E
        if self.stochasticity ==1:
            dEdt+=dEdt*normal(0,0.05,1)
        return dEdt
        
    def dFdt_f(self,X):
        U,V,A,B,C,D,E,F = X
        dFdt= self.bF+self.VF*self.noncompetitiveact(V,self.Kbd,self.nbd)-self.muLVA*F
        if self.stochasticity ==1:
            dFdt+=dFdt*normal(0,0.05,1)
        return dFdt


    function_list = [dUdt_f,dVdt_f,dAdt_f,dBdt_f,dCdt_f,dDdt_f,dEdt_f,dFdt_f]

    def dudt_growth(self,X, cell_matrix):
        function_list = [self.dUdt_f(X),self.dVdt_f(X),self.dAdt_f(X),self.dBdt_f(X),self.dCdt_f(X), self.dDdt_f(X),self.dEdt_f(X),self.dFdt_f(X)]
        dudt = [eq*cell_matrix for eq in function_list]
        return dudt

    def dudt(self,X):
        dudt = [self.dUdt_f(X),self.dVdt_f(X),self.dAdt_f(X),self.dBdt_f(X),self.dCdt_f(X), self.dDdt_f(X),self.dEdt_f(X),self.dFdt_f(X)]
        return dudt

    def getJacobian(self,X,wvn):

        U,V,A,B,C,D,E,F = X
        JU = [-self.DU*wvn**2 - self.muU, 0, self.k1, 0, 0, 0, 0, 0]
        JV = [0, -self.DV*wvn**2 - self.muV, self.k2, 0, 0, 0, 0, 0]
        JA = [0, 0, -self.muASV, 0, 0, -self.VA*self.nda*(D/self.Kda)**self.nda/(D*((D/self.Kda)**self.nda + 1)**2), 0, 0]
        JB = [-self.VB*self.nab*(U/self.Kab)**(2*self.nab)/(U*((E/self.Keb)**self.nab + 1)*((U/self.Kab)**self.nab + 1)**2) + self.VB*self.nab*(U/self.Kab)**self.nab/(U*((E/self.Keb)**self.nab + 1)*((U/self.Kab)**self.nab + 1)), 0, 0, -self.muASV, 0, 0, -self.VB*self.nab*(E/self.Keb)**self.nab*(U/self.Kab)**self.nab/(E*((E/self.Keb)**self.nab + 1)**2*((U/self.Kab)**self.nab + 1)), 0]
        JC = [0, 0, 0, 0, -self.muLVA, -self.VC*self.nda*(D/self.Kda)**self.nda/(D*((D/self.Kda)**self.nda + 1)**2), 0, 0]
        JD = [0, -self.VD*self.nbd*(V/self.Kbd)**(2*self.nbd)/(V*((V/self.Kbd)**self.nbd + 1)**2) + self.VD*self.nbd*(V/self.Kbd)**self.nbd/(V*((V/self.Kbd)**self.nbd + 1)), 0, 0, 0, -self.muLVA, 0, 0]
        JE = [0, 0, 0, 0, -self.VE*self.nce*(C/self.Kce)**self.nce*(E/self.Kee)**self.nee/(C*((C/self.Kce)**self.nce + 1)**2*((E/self.Kee)**self.nee + 1)*((F/self.Kfe)**self.nfe + 1)), 0, -self.muLVA - self.VE*self.nee*(E/self.Kee)**(2*self.nee)/(E*((C/self.Kce)**self.nce + 1)*((E/self.Kee)**self.nee + 1)**2*((F/self.Kfe)**self.nfe + 1)) + self.VE*self.nee*(E/self.Kee)**self.nee/(E*((C/self.Kce)**self.nce + 1)*((E/self.Kee)**self.nee + 1)*((F/self.Kfe)**self.nfe + 1)), -self.VE*self.nfe*(E/self.Kee)**self.nee*(F/self.Kfe)**self.nfe/(F*((C/self.Kce)**self.nce + 1)*((E/self.Kee)**self.nee + 1)*((F/self.Kfe)**self.nfe + 1)**2)]
        JF = [0, -self.VF*self.nbd*(V/self.Kbd)**(2*self.nbd)/(V*((V/self.Kbd)**self.nbd + 1)**2) + self.VF*self.nbd*(V/self.Kbd)**self.nbd/(V*((V/self.Kbd)**self.nbd + 1)), 0, 0, 0, 0, 0, -self.muLVA]
        return np.array([JU,JV,JA,JB,JC,JD,JE,JF ])

#dimensionless equations with 8to6 transition
class circuit14(hill_functions):

    def __init__(self,par_dict,stochasticity=0):
        for key,value in par_dict.items():
            setattr(self,key,value)
        setattr(self, 'stochasticity', stochasticity)


    def dAdt_f(self,species_list, wvn=0):
        A,B,C,D,E,F = species_list
        dadt= 1 + self.VA*self.noncompetitiveinh(D,self.Kda,self.nda) -A - A*wvn**2
        if self.stochasticity ==1:
            dadt+=dadt*normal(0,0.05,1)
        return dadt


    def dBdt_f(self,species_list, wvn=0):
        A,B,C,D,E,F = species_list
        # dbdt= self.muLVA*(1 + self.VB*self.noncompetitivediffact(A,self.Kub,self.nub, self.ku, self.muu)*self.noncompetitiveinh(E,self.Keb, self.neb) - B ) -  B*self.Dr*wvn**2
        dbdt= self.muLVA*(1 + self.VB*self.noncompetitiveact(A,self.Kab, self.nab)*self.noncompetitiveinh(E,self.Keb, self.neb) - B ) -  B*self.Dr*wvn**2
        if self.stochasticity ==1:
            dbdt+=dbdt*normal(0,0.05,1)
        return dbdt


    def dCdt_f(self,species_list):
        A,B,C,D,E,F = species_list
        dcdt= self.muLVA*(1 + self.VC*self.noncompetitiveinh(D,self.Kda,self.nda) - C ) 
        if self.stochasticity ==1:
            dcdt+=dcdt*normal(0,0.05,1)
        return dcdt

    def dDdt_f(self,species_list):
        A,B,C,D,E,F = species_list
        dddt= self.muLVA*(1 + self.VD*self.noncompetitiveact(B,self.Kbd,self.nbd) - D ) 
        if self.stochasticity ==1:
            dddt+=dddt*normal(0,0.05,1)
        return dddt

    def dEdt_f(self,species_list):
        A,B,C,D,E,F = species_list
        dedt= self.muLVA*(1 + self.VE*self.noncompetitiveinh(C,self.Kce,self.nce)*self.noncompetitiveinh(F,self.Kfe,self.nfe)*self.noncompetitiveact(E,self.Kee,self.nee) - E ) 
        if self.stochasticity ==1:
            dedt+=dedt*normal(0,0.05,1)
        return dedt
        
    def dFdt_f(self,species_list):
        A,B,C,D,E,F = species_list
        dfdt= self.muLVA*(1 + self.VF*self.noncompetitiveact(B,self.Kbd,self.nbd) - F ) 
        if self.stochasticity ==1:
            dfdt+=dfdt*normal(0,0.05,1)
        return dfdt

    function_list = [dAdt_f,dBdt_f,dCdt_f,dDdt_f,dEdt_f,dFdt_f]
    
    # @numba.jit(nopython=True)
    def dudt_growth(self,U, cell_matrix):
        function_list = [self.dAdt_f(U),self.dBdt_f(U),self.dCdt_f(U), self.dDdt_f(U),self.dEdt_f(U),self.dFdt_f(U)]
        dudt = [eq*cell_matrix for eq in function_list]
        return dudt
    # @numba.jit(nopython=True)
    def dudt(self,U):
        dudt = [self.dAdt_f(U),self.dBdt_f(U),self.dCdt_f(U), self.dDdt_f(U),self.dEdt_f(U),self.dFdt_f(U)]
        return dudt

    def getJacobian(self,x,wvn):

        A,B,C,D,E,F = x

        JA = [-wvn**2 - 1, 0, 0, -self.VA*self.nda*(D/self.Kda)**self.nda/(D*((D/self.Kda)**self.nda + 1)**2), 0, 0]
        JB = [self.VB*self.muLVA*self.nab*(self.Kab/A)**self.nab/(A*((self.Kab/A)**self.nab + 1)**2*((E/self.Keb)**self.neb + 1)), -self.Dr*wvn**2 - self.muLVA, 0, 0, -self.VB*self.muLVA*self.neb*(E/self.Keb)**self.neb/(E*((self.Kab/A)**self.nab + 1)*((E/self.Keb)**self.neb + 1)**2), 0]
        JC = [0, 0, -self.muLVA, -self.VC*self.muLVA*self.nda*(D/self.Kda)**self.nda/(D*((D/self.Kda)**self.nda + 1)**2), 0, 0]
        JD = [0, self.VD*self.muLVA*self.nbd*(self.Kbd/B)**self.nbd/(B*((self.Kbd/B)**self.nbd + 1)**2), 0, -self.muLVA, 0, 0]
        JE = [0, 0, -self.VE*self.muLVA*self.nce*(C/self.Kce)**self.nce/(C*((C/self.Kce)**self.nce + 1)**2*((self.Kee/E)**self.nee + 1)*((F/self.Kfe)**self.nfe + 1)), 0, self.muLVA*(-1 + self.VE*self.nee*(self.Kee/E)**self.nee/(E*((C/self.Kce)**self.nce + 1)*((self.Kee/E)**self.nee + 1)**2*((F/self.Kfe)**self.nfe + 1))), -self.VE*self.muLVA*self.nfe*(F/self.Kfe)**self.nfe/(F*((C/self.Kce)**self.nce + 1)*((self.Kee/E)**self.nee + 1)*((F/self.Kfe)**self.nfe + 1)**2)]
        JF = [0, self.VF*self.muLVA*self.nbd*(self.Kbd/B)**self.nbd/(B*((self.Kbd/B)**self.nbd + 1)**2), 0, 0, 0, -self.muLVA]
        return np.array([JA, JB, JC, JD, JE, JF])
