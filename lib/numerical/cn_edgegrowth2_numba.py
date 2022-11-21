import numpy as np
from scipy.sparse import spdiags, diags
from tqdm import tqdm
import copy
from scipy.linalg import solve_banded


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
import numba
from numba import cuda, float32


#in edegrowth 2 we have a fixed field and we grow the tissue with a mask. there is diffusion outside the tissue
# new cells start with the concentration similar to the concentration of the neighbouring cells due to diffusion
def cn_edgegrowth2(par_dict,L,J,T,N, circuit_n, steadystate='',growth='linear',rate=1, n_species=2, perturbation=0.01, boundaryCoeff=1, tqdm_disable=False):
    #spatial variables
    dx = float(L)/float(J)
    x_grid = np.array([j*dx for j in range(J)])
    x_gridpoints = J/L
    #temporal variables
    dt = float(T)/float(N)
    t_grid = np.array([n*dt for n in range(N)])
    t_gridpoints =  N/T 

    parent_list = {'circuit1':circuit1, 'circuit2':circuit2,'circuit3':circuit3,'circuit4':circuit4,'circuit5':circuit5, 'circuit6':circuit6, 'circuit7':circuit7, 'schnakenberg':schnakenberg, 'turinghill':turinghill}
    f = parent_list[circuit_n](par_dict, stochasticity=0)

    n_species=2 #number of chemical species/variables/equations
    D = par_dict['d_A'], par_dict['d_B'] 

    
    np.random.seed(1)

    


    steadystate = par_dict['ss_list']
    #create cell matrix
    @numba.jit(nopython=True) #numba here adds a bit of improvement
    def cellMatrixFunction(currentJ,J=J):
        len_half_pad = int((J - currentJ) / 2)
        if (len_half_pad*2+currentJ)==J:
            cellMatrix = np.concatenate((np.zeros(len_half_pad),np.ones(currentJ),np.zeros(len_half_pad)))
        elif (len_half_pad*2+currentJ)<J:
            cellMatrix = np.concatenate((np.zeros(len_half_pad),np.ones(currentJ),np.zeros(len_half_pad+1)))
        else:
            print('ERROR: cellMatrix too large')
        return cellMatrix
    initialJ=1
    cellMatrix=cellMatrixFunction(initialJ)

    U0 = []
    for index in range(n_species):
        U0X=np.random.uniform(low=steadystate[index] - perturbation, high=steadystate[index] + perturbation, size=J)
        U0.append(U0X)   
    U0 = [U0x * cellMatrix for U0x in U0]

    #Look back at mathematical derivation above for understanding of the A and B matrices.

    def alpha(D,dt,dx,n_species):
        return [D[n]*dt/(2.*dx*dx) for n in range(n_species)]

    def A(alphan,J):
        bottomdiag = [-alphan for j in range(J-1)]
        centraldiag = [1.+boundaryCoeff*alphan]+[1.+2.*alphan for j in range(J-2)]+[1.+boundaryCoeff*alphan]
        topdiag = [-alphan for j in range(J-1)]
        diagonals = [bottomdiag,centraldiag,topdiag]
        A = diags(diagonals, [ -1, 0,1]).toarray()
        return A

    def B(alphan,J):
        bottomdiag = [alphan for j in range(J-1)]
        centraldiag = [1.-boundaryCoeff*alphan]+[1.-2.*alphan for j in range(J-2)]+[1.-boundaryCoeff*alphan]
        topdiag = [alphan for j in range(J-1)]
        diagonals = [bottomdiag,centraldiag,topdiag]
        B = diags(diagonals, [ -1, 0,1]).toarray()
        return B



    #storage variables
    record_every_x_hours = 10
    reduced_t_grid = np.arange(0,T,record_every_x_hours)     

    U = copy.deepcopy(U0) 
        #copydeepcopy is useful to make sure the original U0 concentration is not modified and we can retrieve it later on if needed. 
        #we will work with U and U_new from here onwards (U_new is the updated U after calculation).
    U_record=[]
    record_every_x_hours = 10
    for species_index in range(n_species):
        U_record.append(np.zeros([ int(T/record_every_x_hours), J])) #DO NOT SIMPLIFY TO U_record = [np.zeros([J, I, T])]*n_species



    def exponential_growth(t, s=0.0001, initialL=dx):
        return (initialL*np.exp(s*t))


    def linear_growth(t,s=0.00005, initialL=dx):
        return initialL + t*s


    #These two lists contain the A and B matrices for every chemical specie. They are adapted to the size of the field, 
    #meaning that if the field is J=3, the matrix will be 3x3.

    A_list = [A(alphan,J) for alphan in alpha(D,dt,dx,n_species)]  
    B_list = [B(alphan,J) for alphan in alpha(D,dt,dx,n_species)]   
    A_inv = [np.linalg.inv(a) for a in A_list] # Find inverse matrix of A. speeds calculation as only need to get the inverse once. suitable for non-growing systems only

    newL = initialJ*dx

    #for loop iterates over time recalculating the chemical concentrations at each timepoint (ti). 
    for ti in tqdm(range(N), disable = tqdm_disable): 
        previousL = copy.deepcopy(newL)
        hour = ti / (N / T)
        if growth == 'exponential':
            newL = exponential_growth(hour,s=rate)

        if growth == 'linear':
            newL = linear_growth(hour,s=rate)

        if newL != previousL:
            currentJ = int(newL*x_gridpoints)
            cellMatrix = cellMatrixFunction(currentJ)



        U_new = copy.deepcopy(U)
        f0 = f.dudt_growth(U_new, cellMatrix)
        
        #iterate over every chemical specie when calculating concentrations. 
        for n in range(n_species):
            U_new[n] = A_inv[n].dot(B_list[n].dot(U[n]) +  f0[n]*(dt/2)) # Dot product with inverse rather than solve system of equations
            if any(x<0 for x in U_new[n]):
                print('negative')
            
        
        hour = ti / (N / T)
        if hour % record_every_x_hours == 0 :  #only grow and record every 10 hours unit time (hour)
            for n in range(n_species):
                U_record[n][int(hour/record_every_x_hours), :] = U_new[n] #Solution added into array which records the solution over time (JxT dimensional array)
        
        U = copy.deepcopy(U_new)
    


    return U,U_record, U0, x_grid, reduced_t_grid, cellMatrix




