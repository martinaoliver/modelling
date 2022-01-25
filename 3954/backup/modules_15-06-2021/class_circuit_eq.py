import numpy as np

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

    def noncompetitiveact(self, U, km):
        act = ((U / km) ** (self.n)) / (1 + (U / km) ** (self.n))
        return act

    def noncompetitiveinh(self, U, km):
        inh = 1 / (1 + (U / km) ** (self.n))
        return inh

class circuit1_eq(hill_functions):

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


class circuit2_eq(hill_functions):

    def __init__(self,par_dict):
        for key,value in par_dict.items():
            setattr(self,key,value)


    def dAdt_f(self,species_list):
        A,B,C,D,E,F = species_list
        dadt= self.ba+self.Va*self.noncompetitiveinh(D,self.kda)-self.mua*A
        return dadt


    def dBdt_f(self,species_list):
        A,B,C,D,E,F = species_list
        dbdt= self.bb+self.Vb*self.noncompetitiveact(A,self.kaa)*self.noncompetitiveinh(E,self.keb)-self.mua*B
        return dbdt


    def dCdt_f(self,species_list):
        A,B,C,D,E,F = species_list
        dcdt= self.bc+self.Vc*self.noncompetitiveinh(D,self.kda)-self.mulva*C
        return dcdt

    def dDdt_f(self,species_list):
        A,B,C,D,E,F = species_list
        dddt= self.bd+self.Vd*self.noncompetitiveact(B,self.kbd)-self.mulva*D
        return dddt

    def dEdt_f(self,species_list):
        A,B,C,D,E,F = species_list
        dedt= self.be+self.Ve*self.noncompetitiveinh(C,self.kce)*self.noncompetitiveinh(F,self.kfe)*self.noncompetitiveact(E,self.kee)-self.mulva*E
        return dedt
    def dFdt_f(self,species_list):
        A,B,C,D,E,F = species_list
        dfdt= self.bf+self.Vf*self.noncompetitiveact(B,self.kbd)-self.mulva*F
        return dfdt
    function_list = [dAdt_f,dBdt_f,dCdt_f,dDdt_f,dEdt_f,dFdt_f]


    def getJacobian(self,x,wvn):

        A,B,C,D,E,F = x

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

        return np.array([JA, JB, JC, JD, JE, JF])


class circuit3_eq(hill_functions):

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


class circuit4_eq(hill_functions):

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


class circuit5_eq(hill_functions):

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
class circuit6_eq(hill_functions):
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
class circuit7_eq(hill_functions):
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
