from re import S
import numpy as np
from scipy.sparse import spdiags, diags
from tqdm import tqdm
import copy
from scipy.linalg import solve_banded
import matplotlib.pyplot as plt
import numba
from numba import cuda, float32

#############
###paths#####
#############
import sys
import os

from importlib_metadata import distribution
pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')
#############

from equations.class_circuit_eq import *
from equations.twonode_eq import *

def cn_nogrowth(par_dict,L,J,T,N, circuit_n, n_species=2, tqdm_disable=False):
    #spatial variables
    dx = float(L)/float(J-1)
    x_grid = np.array([j*dx for j in range(J)])
    x_gridpoints = J/L
    #temporal variables
    dt = float(T)/float(N-1)
    t_grid = np.array([n*dt for n in range(N)])
    
    dt = float(T)/float(N-1)
    t_grid = np.array([n*dt for n in range(N)])
    t_gridpoints =  N/T 

    parent_list = {'circuit1':circuit1, 'circuit2':circuit2,'circuit3':circuit3,'circuit4':circuit4,'circuit5':circuit5, 'circuit6':circuit6, 'circuit7':circuit7, 'schnakenberg':schnakenberg, 'turinghill':turinghill}
    f = parent_list[circuit_n](par_dict, stochasticity=0)

    n_species=2 #number of chemical species/variables/equations
    D = par_dict['d_A'], par_dict['d_B'] 

    U0 = []
    np.random.seed(1)

    steadystate = par_dict['ss_list']
    for index in range(n_species):
        U0.append(steadystate[index]*(1+np.random.normal(loc=0,scale=0.001,size=J)))

    #Look back at mathematical derivation above for understanding of the A and B matrices.

    def alpha(D,dt,dx,n_species):
        return [D[n]*dt/(2.*dx*dx) for n in range(n_species)]

    def A(alphan,J):
        bottomdiag = [-alphan for j in range(J-1)]
        centraldiag = [1.+alphan]+[1.+2.*alphan for j in range(J-2)]+[1.+alphan]
        topdiag = [-alphan for j in range(J-1)]
        diagonals = [bottomdiag,centraldiag,topdiag]
        A = diags(diagonals, [ -1, 0,1]).toarray()
        return A

    def B(alphan,J):
        bottomdiag = [alphan for j in range(J-1)]
        centraldiag = [1.-alphan]+[1.-2.*alphan for j in range(J-2)]+[1.-alphan]
        topdiag = [alphan for j in range(J-1)]
        diagonals = [bottomdiag,centraldiag,topdiag]
        B = diags(diagonals, [ -1, 0,1]).toarray()
        return B

    record_every_x_hours = 10

    #storage variables
    reduced_t_grid = np.arange(0,T,record_every_x_hours) 
    U = copy.deepcopy(U0) 
        #copydeepcopy is useful to make sure the original U0 concentration is not modified and we can retrieve it later on if needed. 
        #we will work with U and U_new from here onwards (U_new is the updated U after calculation).
    U_record=[]
    record_every_x_hours = 10
    for species_index in range(n_species):
        U_record.append(np.zeros([ int(T/record_every_x_hours), J])) #DO NOT SIMPLIFY TO U_record = [np.zeros([J, I, T])]*n_species


    #These two lists contain the A and B matrices for every chemical specie. They are adapted to the size of the field, 
    #meaning that if the field is J=3, the matrix will be 3x3.
    A_list = [A(alphan,J) for alphan in alpha(D,dt,dx,n_species)]  
    B_list = [B(alphan,J) for alphan in alpha(D,dt,dx,n_species)]   
    A_inv = [np.linalg.inv(a) for a in A_list] # Find inverse matrix of A. speeds calculation as only need to get the inverse once. suitable for non-growing systems only
    
    numba.jit(nopython=True)
    def cn_forloop(U,N,T,A_inv,B_list, record_every_x_hours):
        #for loop iterates over time recalculating the chemical concentrations at each timepoint (ti). 
        print('entering numba for loop')
        for ti in range(N): 
            U_new = U.copy()
            f0 = f.dudt(U_new)

            #iterate over every chemical specie when calculating concentrations. 
            for n in range(n_species):
                U_new[n] = A_inv[n].dot(B_list[n].dot(U[n]) +  f0[n]*(dt/2)) # Dot product with inverse rather than solve system of equations


            #storing results
            hour = ti / (N / T)


            if hour % record_every_x_hours == 0 :  #only grow and record every 10 hours unit time (hour)
                for n in range(n_species):
                    U_record[n][int(hour/record_every_x_hours), :] = U_new[n] #Solution added into array which records the solution over time (JxT dimensional array)
            U = U_new.copy()
        return U,U_record
    
    U,U_record = cn_forloop(U,N,T,A_inv,B_list,record_every_x_hours)
    
    return U,U_record, U0, x_grid, reduced_t_grid



