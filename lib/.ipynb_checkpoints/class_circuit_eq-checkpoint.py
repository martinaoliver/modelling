import numpy as np


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


        def getJacobian_nodiffusion(self,x):
            J1 = [(self.Va*self.n*(x[0]/self.kaa)**(self.n - 1))/(self.kaa*((x[0]/self.kaa)**self.n + 1)*((x[3]/self.kda)**self.n + 1)) - self.mua - (self.Va*self.n*(x[0]/self.kaa)**self.n*(x[0]/self.kaa)**(self.n - 1))/(self.kaa*((x[0]/self.kaa)**self.n + 1)**2*((x[3]/self.kda)**self.n + 1)),                                                                                                       0,                                                                                         0, -(self.Va*self.n*(x[0]/self.kaa)**self.n*(x[3]/self.kda)**(self.n - 1))/(self.kda*((x[0]/self.kaa)**self.n + 1)*((x[3]/self.kda)**self.n + 1)**2),                                                                                                                                                                               0,                                                                                         0]
            J2 = [      (self.Vb*self.n*(x[0]/self.kaa)**(self.n - 1))/(self.kaa*((x[0]/self.kaa)**self.n + 1)*((x[4]/self.keb)**self.n + 1)) - (self.Vb*self.n*(x[0]/self.kaa)**self.n*(x[0]/self.kaa)**(self.n - 1))/(self.kaa*((x[0]/self.kaa)**self.n + 1)**2*((x[4]/self.keb)**self.n + 1)),                                                                                                    -self.mub,                                                                                         0,                                                                         0,                                                                                                       -(self.Vb*self.n*(x[0]/self.kaa)**self.n*(x[4]/self.keb)**(self.n - 1))/(self.keb*((x[0]/self.kaa)**self.n + 1)*((x[4]/self.keb)**self.n + 1)**2),                                                                                         0]
            J3 = [      (self.Vc*self.n*(x[0]/self.kaa)**(self.n - 1))/(self.kaa*((x[0]/self.kaa)**self.n + 1)*((x[3]/self.kda)**self.n + 1)) - (self.Vc*self.n*(x[0]/self.kaa)**self.n*(x[0]/self.kaa)**(self.n - 1))/(self.kaa*((x[0]/self.kaa)**self.n + 1)**2*((x[3]/self.kda)**self.n + 1)),                                                                                                       0,                                                                                    -self.mulva, -(self.Vc*self.n*(x[0]/self.kaa)**self.n*(x[3]/self.kda)**(self.n - 1))/(self.kda*((x[0]/self.kaa)**self.n + 1)*((x[3]/self.kda)**self.n + 1)**2),                                                                                                                                                                               0,                                                                                         0]
            J4 = [                                                                                                                                            0, (self.Vd*self.n*(x[1]/self.kbd)**(self.n - 1))/(self.kbd*((x[1]/self.kbd)**self.n + 1)) - (self.Vd*self.n*(x[1]/self.kbd)**self.n*(x[1]/self.kbd)**(self.n - 1))/(self.kbd*((x[1]/self.kbd)**self.n + 1)**2),                                                                                         0,                                                                    -self.mulva,                                                                                                                                                                               0,                                                                                         0]
            J5 = [                                                                                                                                            0,                                                                                                       0, -(self.Ve*self.n*(x[2]/self.kce)**(self.n - 1)*(x[4]/self.kee)**self.n)/(self.kce*((x[2]/self.kce)**self.n + 1)**2*((x[4]/self.kee)**self.n + 1)*((x[5]/self.kfe)**self.n + 1)),                                                                         0, (self.Ve*self.n*(x[4]/self.kee)**(self.n - 1))/(self.kee*((x[2]/self.kce)**self.n + 1)*((x[4]/self.kee)**self.n + 1)*((x[5]/self.kfe)**self.n + 1)) - self.mulva - (self.Ve*self.n*(x[4]/self.kee)**self.n*(x[4]/self.kee)**(self.n - 1))/(self.kee*((x[2]/self.kce)**self.n + 1)*((x[4]/self.kee)**self.n + 1)**2*((x[5]/self.kfe)**self.n + 1)), -(self.Ve*self.n*(x[4]/self.kee)**self.n*(x[5]/self.kfe)**(self.n - 1))/(self.kfe*((x[2]/self.kce)**self.n + 1)*((x[4]/self.kee)**self.n + 1)*((x[5]/self.kfe)**self.n + 1)**2)]
            J6 = [                                                                                                                                            0, (self.Vf*self.n*(x[1]/self.kbd)**(self.n - 1))/(self.kbd*((x[1]/self.kbd)**self.n + 1)) - (self.Vf*self.n*(x[1]/self.kbd)**self.n*(x[1]/self.kbd)**(self.n - 1))/(self.kbd*((x[1]/self.kbd)**self.n + 1)**2),                                                                                         0,                                                                         0,                                                                                                                                                                               0,                                                                                    -self.mulva]
            return  np.array([J1,J2,J3,J4,J5,J6])

        def getJacobian_diffusion(self, A, B, C, D, E, F, wvn):  # circuit1_eq
            JA = [
                (self.Va * self.n * (A / self.kaa) ** (self.n - 1)) / (self.kaa * ((A / self.kaa) ** self.n + 1) * (
                        (D / self.kda) ** self.n + 1)) - self.d_A * wvn ** 2 - self.mua - (
                        self.Va * self.n * (A / self.kaa) ** self.n * (A / self.kaa) ** (self.n - 1)) / (
                        self.kaa * ((A / self.kaa) ** self.n + 1) ** 2 * ((D / self.kda) ** self.n + 1)), 0, 0,
                -(self.Va * self.n * (A / self.kaa) ** self.n * (D / self.kda) ** (self.n - 1)) / (
                        self.kda * ((A / self.kaa) ** self.n + 1) * ((D / self.kda) ** self.n + 1) ** 2), 0, 0]
            JB = [(self.Vb * self.n * (A / self.kaa) ** (self.n - 1)) / (
                    self.kaa * ((A / self.kaa) ** self.n + 1) * ((E / self.keb) ** self.n + 1)) - (
                          self.Vb * self.n * (A / self.kaa) ** self.n * (A / self.kaa) ** (self.n - 1)) / (
                          self.kaa * ((A / self.kaa) ** self.n + 1) ** 2 * ((E / self.keb) ** self.n + 1)),
                  - self.d_B * wvn ** 2 - self.mub, 0, 0,
                  -(self.Vb * self.n * (A / self.kaa) ** self.n * (E / self.keb) ** (self.n - 1)) / (
                          self.keb * ((A / self.kaa) ** self.n + 1) * ((E / self.keb) ** self.n + 1) ** 2), 0]
            JC = [(self.Vc * self.n * (A / self.kaa) ** (self.n - 1)) / (
                    self.kaa * ((A / self.kaa) ** self.n + 1) * ((D / self.kda) ** self.n + 1)) - (
                          self.Vc * self.n * (A / self.kaa) ** self.n * (A / self.kaa) ** (self.n - 1)) / (
                          self.kaa * ((A / self.kaa) ** self.n + 1) ** 2 * ((D / self.kda) ** self.n + 1)), 0,
                  -self.mulva, -(self.Vc * self.n * (A / self.kaa) ** self.n * (D / self.kda) ** (self.n - 1)) / (
                          self.kda * ((A / self.kaa) ** self.n + 1) * ((D / self.kda) ** self.n + 1) ** 2), 0, 0]
            JD = [0,
                  (self.Vd * self.n * (B / self.kbd) ** (self.n - 1)) / (
                              self.kbd * ((B / self.kbd) ** self.n + 1)) - (
                          self.Vd * self.n * (B / self.kbd) ** self.n * (B / self.kbd) ** (self.n - 1)) / (
                          self.kbd * ((B / self.kbd) ** self.n + 1) ** 2), 0, -self.mulva, 0, 0]
            JE = [0, 0, -(self.Ve * self.n * (C / self.kce) ** (self.n - 1) * (E / self.kee) ** self.n) / (
                    self.kce * ((C / self.kce) ** self.n + 1) ** 2 * ((E / self.kee) ** self.n + 1) * (
                    (F / self.kfe) ** self.n + 1)), 0, (self.Ve * self.n * (E / self.kee) ** (self.n - 1)) / (
                          self.kee * ((C / self.kce) ** self.n + 1) * ((E / self.kee) ** self.n + 1) * (
                          (F / self.kfe) ** self.n + 1)) - self.mulva - (
                          self.Ve * self.n * (E / self.kee) ** self.n * (E / self.kee) ** (self.n - 1)) / (
                          self.kee * ((C / self.kce) ** self.n + 1) * ((E / self.kee) ** self.n + 1) ** 2 * (
                          (F / self.kfe) ** self.n + 1)),
                  -(self.Ve * self.n * (E / self.kee) ** self.n * (F / self.kfe) ** (self.n - 1)) / (
                          self.kfe * ((C / self.kce) ** self.n + 1) * ((E / self.kee) ** self.n + 1) * (
                          (F / self.kfe) ** self.n + 1) ** 2)]
            JF = [0,
                  (self.Vf * self.n * (B / self.kbd) ** (self.n - 1)) / (
                              self.kbd * ((B / self.kbd) ** self.n + 1)) - (
                          self.Vf * self.n * (B / self.kbd) ** self.n * (B / self.kbd) ** (self.n - 1)) / (
                          self.kbd * ((B / self.kbd) ** self.n + 1) ** 2), 0, 0, 0, -self.mulva]
            return np.array([JA, JB, JC, JD, JE, JF])
class circuit2_eq(hill_functions):

    def __init__(self,par_dict):
        for key,value in par_dict.items():
            setattr(self,key,value)


    def dAdt_f(self,A,B,C,D,E,F):
        dadt= self.ba+self.Va*self.noncompetitiveinh(D,self.kda)-self.mua*A
        return dadt


    def dBdt_f(self,A,B,C,D,E,F):
        dbdt= self.bb+self.Vb*self.noncompetitiveact(A,self.kaa)*self.noncompetitiveinh(E,self.keb)-self.mua*B
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


    def getJacobian_nodiffusion(self,x):
        J1 = [-self.mua, 0, 0, -(self.Va * self.n * (x[3] / self.kda) ** (self.n - 1)) / (self.kda * ((x[3] / self.kda) ** self.n + 1) ** 2), 0, 0]

        J2 = [(self.Vb * self.n * (x[0] / self.kaa) ** (self.n - 1)) / ( self.kaa * ((x[0] / self.kaa) ** self.n + 1) * ((x[4] / self.keb) ** self.n + 1)) - (
                      self.Vb * self.n * (x[0] / self.kaa) ** self.n * (x[0] / self.kaa) ** (self.n - 1)) / (
                      self.kaa * ((x[0] / self.kaa) ** self.n + 1) ** 2 * ((x[4] / self.keb) ** self.n + 1)), -self.mua, 0,  0,
              -(self.Vb * self.n * (x[0] / self.kaa) ** self.n * (x[4] / self.keb) ** (self.n - 1)) / ( self.keb * ((x[0] / self.kaa) ** self.n + 1) * ((x[4] / self.keb) ** self.n + 1) ** 2), 0]

        J3 = [(self.Vc * self.n * (x[0] / self.kaa) ** (self.n - 1)) / ( self.kaa * ((x[0] / self.kaa) ** self.n + 1) * ((x[3] / self.kda) ** self.n + 1)) - (
                      self.Vc * self.n * (x[0] / self.kaa) ** self.n * (x[0] / self.kaa) ** (self.n - 1)) / (
                      self.kaa * ((x[0] / self.kaa) ** self.n + 1) ** 2 * ((x[3] / self.kda) ** self.n + 1)), 0,
              -self.mulva, -(self.Vc * self.n * (x[0] / self.kaa) ** self.n * (x[3] / self.kda) ** (self.n - 1)) / (
                          self.kda * ((x[0] / self.kaa) ** self.n + 1) * ((x[3] / self.kda) ** self.n + 1) ** 2), 0, 0]

        J4 = [0, (self.Vd * self.n * (x[1] / self.kbd) ** (self.n - 1)) / (self.kbd * ((x[1] / self.kbd) ** self.n + 1)) - (
                      self.Vd * self.n * (x[1] / self.kbd) ** self.n * (x[1] / self.kbd) ** (self.n - 1)) / (
                          self.kbd * ((x[1] / self.kbd) ** self.n + 1) ** 2), 0, -self.mulva, 0, 0]

        J5 = [0, 0, -(self.Ve * self.n * (x[2] / self.kce) ** (self.n - 1) * (x[4] / self.kee) ** self.n) / (
                self.kce * ((x[2] / self.kce) ** self.n + 1) ** 2 * ((x[4] / self.kee) ** self.n + 1) * (
                    (x[5] / self.kfe) ** self.n + 1)), 0,  (self.Ve * self.n * (x[4] / self.kee) ** (self.n - 1)) / (
                      self.kee * ((x[2] / self.kce) ** self.n + 1) * ((x[4] / self.kee) ** self.n + 1) * (
                          (x[5] / self.kfe) ** self.n + 1)) - self.mulva - ( self.Ve * self.n * (x[4] / self.kee) ** self.n * (x[4] / self.kee) ** (self.n - 1)) / (
                      self.kee * ((x[2] / self.kce) ** self.n + 1) * ((x[4] / self.kee) ** self.n + 1) ** 2 * (
                          (x[5] / self.kfe) ** self.n + 1)),   -(self.Ve * self.n * (x[4] / self.kee) ** self.n * (x[5] / self.kfe) ** (self.n - 1)) / (
                      self.kfe * ((x[2] / self.kce) ** self.n + 1) * ((x[4] / self.kee) ** self.n + 1) * (   (x[5] / self.kfe) ** self.n + 1) ** 2)]
        J6 = [0,   (self.Vf * self.n * (x[1] / self.kbd) ** (self.n - 1)) / (self.kbd * ((x[1] / self.kbd) ** self.n + 1)) - (
                      self.Vf * self.n * (x[1] / self.kbd) ** self.n * (x[1] / self.kbd) ** (self.n - 1)) / (
                          self.kbd * ((x[1] / self.kbd) ** self.n + 1) ** 2), 0, 0, 0, -self.mulva]
        return  np.array([J1,J2,J3,J4,J5,J6])

    def getJacobian_diffusion(self, A, B, C, D, E, F, wvn):  # circuit1_eq
        JA = [- self.d_A * wvn ** 2 - self.mua, 0, 0, -(self.Va * self.n * (D / self.kda) ** (self.n - 1)) / (
                self.kda * ((D / self.kda) ** self.n + 1) ** 2), 0, 0]
        JB = [(self.Vb * self.n * (A / self.kaa) ** (self.n - 1)) / (
                self.kaa * ((A / self.kaa) ** self.n + 1) * ((E / self.keb) ** self.n + 1)) - (
                      self.Vb * self.n * (A / self.kaa) ** self.n * (A / self.kaa) ** (self.n - 1)) / (
                      self.kaa * ((A / self.kaa) ** self.n + 1) ** 2 * ((E / self.keb) ** self.n + 1)),
              - self.d_B * wvn ** 2 - self.mua, 0, 0,
              -(self.Vb * self.n * (A / self.kaa) ** self.n * (E / self.keb) ** (self.n - 1)) / (
                      self.keb * ((A / self.kaa) ** self.n + 1) * ((E / self.keb) ** self.n + 1) ** 2), 0]
        JC = [(self.Vc * self.n * (A / self.kaa) ** (self.n - 1)) / (
                self.kaa * ((A / self.kaa) ** self.n + 1) * ((D / self.kda) ** self.n + 1)) - (
                      self.Vc * self.n * (A / self.kaa) ** self.n * (A / self.kaa) ** (self.n - 1)) / (
                      self.kaa * ((A / self.kaa) ** self.n + 1) ** 2 * ((D / self.kda) ** self.n + 1)), 0,
              -self.mulva, -(self.Vc * self.n * (A / self.kaa) ** self.n * (D / self.kda) ** (self.n - 1)) / (
                      self.kda * ((A / self.kaa) ** self.n + 1) * ((D / self.kda) ** self.n + 1) ** 2), 0, 0]
        JD = [0,
              (self.Vd * self.n * (B / self.kbd) ** (self.n - 1)) / (self.kbd * ((B / self.kbd) ** self.n + 1)) - (
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
        JF = [0,
              (self.Vf * self.n * (B / self.kbd) ** (self.n - 1)) / (self.kbd * ((B / self.kbd) ** self.n + 1)) - (
                      self.Vf * self.n * (B / self.kbd) ** self.n * (B / self.kbd) ** (self.n - 1)) / (
                      self.kbd * ((B / self.kbd) ** self.n + 1) ** 2), 0, 0, 0, -self.mulva]
        return np.array([JA, JB, JC, JD, JE, JF])

