import numpy
# import matplotlib as mpl
# mpl.use('tkagg')

import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import spdiags, diags
import matplotlib.pyplot as plt
from tqdm import tqdm
import copy
from scipy.sparse import linalg
from scipy.linalg import solve_banded
from class_circuit_eq import *

def adi(par_dict,L_x,L_y,J,I,T,N, circuit_n, n_species,D,tqdm_disable=False,stochasticity=0):

    parent_list = [circuit1_eq, circuit2_eq,circuit3_eq,circuit4_eq,circuit5_eq,circuit6_eq,circuit7_eq,circuit8_eq,circuit9_eq, circuit10_eq, circuit11_eq]
    f = parent_list[circuit_n-1](par_dict, stochasticity=stochasticity)

    #spatial variables
    dx = float(L_x)/float(J-1); dy = float(L_y)/float(I-1)
    x_grid = numpy.array([j*dx for j in range(J)]); y_grid = numpy.array([i*dy for i in range(I)])
    diffusing_species =np.nonzero(D)[0]
    nondiffusing_species = np.nonzero(D==0)[0]

    #time variables
    
    dt = float(T)/float(N-1)

    t_grid = numpy.array([n*dt for n in range(N)])

    alpha = [D[n]*dt/(2.*dx*dx) for n in range(n_species)]
    print(dt/(dx)**2)


    #Define initial conditions and cell matrix
    U0 = []
    perturbation=0.001
    steadystates=[0.1]*n_species
    np.random.seed(1)

    cell_matrix = np.zeros(shape=(I,J))
    cell_matrix[int(I/2), int(J/2)] = 1
    for index in range(n_species):
        U0.append(np.random.uniform(low=steadystates[index] - perturbation, high=steadystates[index] + perturbation, size=(I, J)))
    # U0 = U0*cell_matrix

    #A matrix (right-hand side of Ax=b)
    def A(alphan):
        bottomdiag = [-alphan for j in range(J-1)]
        centraldiag = [1.+alphan]+[1.+2.*alphan for j in range(J-2)]+[1.+alphan]
        topdiag = [-alphan for j in range(J-1)]
        diagonals = [bottomdiag,centraldiag,topdiag]
        A = diags(diagonals, [ -1, 0,1]).toarray()
        return A
    def diagonal_form(a, upper = 1, lower= 1):
        """
        a is a numpy square matrix
        this function converts a square matrix to diagonal ordered form
        returned matrix in ab shape which can be used directly for scipy.linalg.solve_banded
        """
        n = a.shape[1]
        assert(np.all(a.shape ==(n,n)))

        ab = np.zeros((2*n-1, n))

        for i in range(n):
            ab[i,(n-1)-i:] = np.diagonal(a,(n-1)-i)

        for i in range(n-1):
            ab[(2*n-2)-i,:i+1] = np.diagonal(a,i-(n-1))

        mid_row_inx = int(ab.shape[0]/2)
        upper_rows = [mid_row_inx - i for i in range(1, upper+1)]
        upper_rows.reverse()
        upper_rows.append(mid_row_inx)
        lower_rows = [mid_row_inx + i for i in range(1, lower+1)]
        keep_rows = upper_rows+lower_rows
        ab = ab[keep_rows,:]


        return ab
    ab_list = [diagonal_form(A(alphan)) for alphan in alpha]
    A_list = [A(alphan) for alphan in alpha]


    #b vector (left-hand side of Ax=b)
    def b(axis,ij,alphan,Un):
        b_t_stencil = np.array( [0] + [(1-alphan)] + [alphan])
        b_c_stencil = np.array( [alphan] + [(1-2*alphan)] + [alphan])
        b_b_stencil = np.array( [alphan] + [(1-alphan)] + [0])

        b = np.zeros(J)
        if axis == 'y':
            i = ij
            if i > 0 and i < I-1:
                for j in range(0,J):
                    ux_three = [Un[j,i-1], Un[j,i], Un[j,i+1]]
                    sub_b = np.sum(ux_three*b_c_stencil)
                    b[j] = sub_b


            if i == 0:
                for j in range(0,J):
                    ux_three = [0 , Un[j,i], Un[j,i+1]]
                    sub_b = np.sum(ux_three*b_t_stencil)
                    b[j] = sub_b

            if i == I-1:
                for j in range(0,J):
                    ux_three = [Un[j,i-1], Un[j,i] , 0]
                    sub_b = np.sum(ux_three*b_b_stencil)
                    b[j] = sub_b

        if axis == 'x':
            j = ij
            if j > 0 and  j < J-1:
                for i in range(0,I):
                    uy_three = [Un[j-1,i], Un[j,i], Un[j+1,i]]
                    sub_b = np.sum(uy_three*b_c_stencil)
                    b[i] = sub_b

            if j == 0:
                for i in range(0,I):
                    uy_three = [0, Un[j,i], Un[j+1,i]]
                    sub_b = np.sum(uy_three*b_t_stencil)
                    b[i] = sub_b

            if j == J-1:
                for i in range(0,I):
                    uy_three = [Un[j-1,i], Un[j,i], 0]
                    sub_b = np.sum(uy_three*b_b_stencil)
                    b[i] = sub_b

        return b



    U = copy.deepcopy(U0)
    U_record = []
    for species_index in range(n_species):
        U_record.append(np.zeros([J, I, T])) #DO NOT SIMPLIFY TO U_record = [np.zeros([J, I, T])]*n_species

    A_inv = [np.linalg.inv(a) for a in A_list]
    unittime=0
    for ti in tqdm(range(N), disable = tqdm_disable):
        #First step: solve in y direction from n -> n+1/2
        U_half = copy.deepcopy(U)
        f0 = f.dudt(U)
        for i in range(I):
            for n in diffusing_species:
                # U_half[n][:,i] = solve_banded((1, 1), ab_list[n], b('y',i,alpha[n],U[n]) +  f0[n][:,i]*(dt/2)) #CN step in one dimension to get banded(tridiagonal) A matrix
                U_half[n][:,i] = A_inv[n].dot(b('y',i,alpha[n],U[n])+  f0[n][:,i]*(dt/2)) # Dot product with inverse rather than solve system of equations
            for n in nondiffusing_species: #Species with D=0, no CN included, basic euler
                U_half[n][:,i] =  U[n][:,i] + f0[n][:,i]*(dt/2)

        #Second step: solve in x direction from n+1/2 -> n+1
        U_new = copy.deepcopy(U_half)
        f1 = f.dudt(U_half)
        for j in range(J):
            for n in diffusing_species:
                # U_new[n][j,:] = solve_banded((1, 1), ab_list[n], b('x',j,alpha[n],U_half[n]) + f1[n][j,:]*(dt/2))
                U_new[n][j,:] = A_inv[n].dot(b('x',j,alpha[n],U_half[n])+  f1[n][j,:]*(dt/2)) # Dot product with inverse rather than solve system of equations
            for n in nondiffusing_species:
                U_new[n][j,:] =  U_half[n][j,:] + f1[n][j,:]*(dt/2)

        hour = ti / (N / T)

        if hour % 1 == 0:  #only consider division at unit time (hour)
            #append results into top_array for records
            for species_index in range(n_species):
                U_record[species_index][:, :, int(hour)] = U_new[species_index] #issue in this line



        U = copy.deepcopy(U_new)

    return U_record, U
