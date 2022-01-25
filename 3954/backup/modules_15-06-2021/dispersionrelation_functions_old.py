import numpy as np
from numpy import linalg as LA
from class_circuit_eq import *


#This class can return the jacobian which is stored in circuit*_eq class.
class jacobian():
    def __init__(self,par_dict, circuit_n):
        for key, value in par_dict.items():
            setattr(self, key, value)
        setattr(self, 'circuit_n', circuit_n)

        self.parent_list = [circuit1_eq, circuit2_eq, circuit3_eq, circuit4_eq, circuit5_eq, circuit6_eq, circuit7_eq]

    def getJacobian(self,x,wvn):  # circuit1_eq
        return self.parent_list[self.circuit_n-1].getJacobian(self,x,wvn)


def stability_no_diffusion(par_dict,circuit_n, x):

    # - Obtain the jacobian (pre-calculated) of the corresponding circuit. Must input the model parameters (par_dict),
    # the steady state and the wavenumber (in this case 0 as there is no diffusion).
    jac = jacobian(par_dict, circuit_n).getJacobian(x, 0)

    eigenvalsteadystate, eigenvec = LA.eig(jac)#calculate the eigenvalues of the jacobian without diffusion
    # print(eigenvalsteadystate)
    if np.any(eigenvalsteadystate.real>0):
       # print('Unstable without diffusion')
       return 'Unstable without diffusion', eigenvalsteadystate
    else:
       return 'Stable without diffusion', eigenvalsteadystate

def stability_diffusion(par_dict,circuit_n, x,top_dispersion):
    n_species = len(x)

    # - Define which wavenumbers will be analysed. L = 100 in this case (100mm) (The experimental system is a 10cm plate)
    # so wavelengths bigger than that are of no interest.
    # - In this case we will sample 5000+1 different wavenumbers. If you sample less you might not find your turing instability.
    # If you sample more, is more computationally expensive.

    # wvn_list = np.array(list(range(0,5000+1)))*np.pi/100
    wvn_list = np.array(list(range(0,top_dispersion+1)))*np.pi/100
    count = 0
    eigenvalues = np.zeros((len(wvn_list),n_species) ,dtype=np.complex_)

    for wvn in wvn_list:
        jac = jacobian(par_dict, circuit_n).getJacobian(x, wvn) #obtain jacobian for corresponding system. This time with a determined wvn.
        eigenval, eigenvec = LA.eig(jac) #calculate the eigenvalues of the jacobian with diffusion
        # sort eigenvalues so the one with the instability is at position -1.
        idx = np.argsort(eigenval)
        eigenval= eigenval[idx]
        eigenvalues[count]  = eigenval
        oscillations = np.any(np.iscomplex(eigenvalues))
        count +=1


    maxeig = np.amax(eigenvalues) #highest eigenvalue of all wvn's and all the 6 eigenvalues.
    maxeig_real = np.amax(eigenvalues.real) #real part of maxeig

    #Applying Turing instability criteria
    if maxeig_real <= 0:
        # print('Stable with diffusion')
        pattern_class = 'Stable with diffusion'

    elif maxeig_real > 0:#Unstability with diffusion
        if np.any(eigenvalues[-1,:] == maxeig_real): #The maximum eigenvalue is also biggest wvn.
            # print('Turing II (instability increases with wavenumber (infinetely))')
            pattern_class = 'Turing II'

        elif np.all(eigenvalues[-1,:] != maxeig_real): #highest instability does not appear with highest wavenumber)
            if np.any(eigenvalues.imag[:,-1]>0): #if the last eigenvalue contains imaginary numbers (the last eigenvalue is that with the instability)
                # print('Oscillatory pattern')
                pattern_class = 'Turing Oscillatory'

            elif np.all(eigenvalues.imag[:,-1]<=0): #if the last eigenvalue does not contain imaginary numbers
                # print ('Turing I pattern')
                pattern_class = 'Turing I'

    return maxeig, pattern_class, eigenvalues,oscillations

def dispersionrelation(par_dict,steadystate_values_ss_n,circuit_n,calculate_unstable,top_dispersion):

    x = steadystate_values_ss_n

    #Stability without diffusion
    steadystate_stability , eigenvalsteadystate = stability_no_diffusion(par_dict,circuit_n, x)

    # - The first requirement for a Turing instability is to be stable without diffusion. If this requirement is not met,
    # we discard this steady state.
    if calculate_unstable==True:
        maxeig, pattern_class, eigenvalues, oscillations = stability_diffusion(par_dict,circuit_n, x,top_dispersion)

    else:
        if steadystate_stability == 'Unstable without diffusion':
           maxeig = float('NaN')
           pattern_class = 'Unstable without diffusion'
           oscillations = float('NaN')
           eigenvalues = []

        # - If the first requirement is met, we can now compute the dispersion relation to see how diffusion affects stability.
        #Stability with diffusion (dispersion relation)
        else:
            maxeig, pattern_class, eigenvalues, oscillations = stability_diffusion(par_dict,circuit_n, x,top_dispersion)

    return maxeig, pattern_class, eigenvalues, oscillations, eigenvalsteadystate
