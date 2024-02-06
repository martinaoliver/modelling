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
# from equations.class_circuit_eq import *
from scipy.ndimage import laplace
import numba
from numba import cuda, float32
# @numba.jit(nopython=True)

# def adi(par_dict,L_x,L_y,J,I,T,N, circuit_n, n_species,D,tqdm_disable=False,stochasticity=0, steadystates=0):
def adi(L_x,L_y,J,I,T,N, n_species,D,tqdm_disable=False,stochasticity=0, steadystates=0):
    #for dt/dx^2 <1 (stability criterion): t_gridpoints approx < xgridpoints^2
    # parent_list = {'circuit1':circuit1, 'circuit2':circuit2,'circuit3':circuit3,'circuit4':circuit4,'circuit5':circuit5, 'circuit6':circuit6, 'circuit7':circuit7, 'turinghill':turinghill}
    # f = parent_list[circuit_n](par_dict, stochasticity=stochasticity)
    def f_Turing(U,n_species=2):
        dudt = [0]*n_species
        dudt[0]= 5*U[0] - 6*U[1] + 1
        dudt[1] = 6*U[0]- 7*U[1] + 1
        return dudt

    @numba.jit(nopython=True)
    def f(c, t, f_args):
        A, B = f_args
        u = c[0, :, :]
        v = c[1, :, :]
        u2 = u**2
        u2v = u2 * v
        fu = A - (B + 1) * u + u2v
        fv = B * u - u2v
        return np.stack((fu, fv))


    # @numba.jit(nopython=True)
    # def f(U,a=2,b=3):
    #     dudt = [0]*2
    #     dudt[0]= a-(b+1)*U[0] + (U[0]**2)*U[1]
    #     dudt[1] = b*U[0] - (U[0]**2)*U[1]
    #     return dudt
    D=[0.02,0.4]
    #spatial variables
    dx = float(L_x)/float(J-1); dy = float(L_y)/float(I-1)
    x_grid = numpy.array([j*dx for j in range(J)]); y_grid = numpy.array([i*dy for i in range(I)])
    diffusing_species =np.nonzero(D)[0]
    nondiffusing_species = np.nonzero(D==0)[0]

    #time variables
    
    dt = float(T)/float(N-1)

    t_grid = numpy.array([n*dt for n in range(N)])

    alpha = [D[n]*dt/(2.*dx*dx) for n in range(n_species)]

    
    #Define initial conditions and cell matrix
    U0 = []
    perturbation=0.001
    np.random.seed(1)

    # if np.all(steadystates)==0:
    #     steadystates=[0.1]*n_species
    # cell_matrix = np.zeros(shape=(I,J))
    # cell_matrix[int(I/2), int(J/2)] = 1
    steadystates = [2,3/2]
    # for index in range(n_species):
        # U0.append(np.random.uniform(low=steadystates[index] - perturbation, high=steadystates[index] + perturbation, size=(I, J)))
    U0 =  np.random.normal(scale=.1, size=(2, I, J)) 
    U0[0, :, :] += 2
    U0[1, :, :] += 3/2
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
    @numba.jit(nopython=True)
    def b(axis,ij,alphan,Un):
        b_t_stencil = np.array( [0] + [(1-alphan)] + [alphan])
        b_c_stencil = np.array( [alphan] + [(1-2*alphan)] + [alphan])
        b_b_stencil = np.array( [alphan] + [(1-alphan)] + [0])

        b = np.zeros(J)
        if axis == 'y':
            i = ij
            if i > 0 and i < I-1:
                for j in range(0,J):
                    # ux_three = [Un[j,i-1], Un[j,i], Un[j,i+1]]
                    ux_three = np.array( [Un[j,i-1]] +  [Un[j,i]]+  [Un[j,i+1]])
                    sub_b = sum(ux_three*b_c_stencil)
                    b[j] = sub_b

            if i == 0:
                for j in range(0,J):
                    ux_three = np.array([0] +[Un[j,i]]+[Un[j,i+1]])
                    sub_b = sum(ux_three*b_t_stencil)
                    b[j] = sub_b

            if i == I-1:
                for j in range(0,J):
                    ux_three = np.array([Un[j,i-1]] + [Un[j,i]] + [0])
                    sub_b = sum(ux_three*b_b_stencil)
                    b[j] = sub_b

        if axis == 'x':
            j = ij
            if j > 0 and  j < J-1:
                for i in range(0,I):
                    uy_three = np.array([Un[j-1,i]]+ [Un[j,i]] + [Un[j+1,i]])
                    sub_b = sum(uy_three*b_c_stencil)
                    b[i] = sub_b

            if j == 0:
                for i in range(0,I):
                    uy_three = np.array([0] +[Un[j,i]] +  [Un[j+1,i]])
                    sub_b = sum(uy_three*b_t_stencil)
                    b[i] = sub_b

            if j == J-1:
                for i in range(0,I):
                    uy_three = np.array([Un[j-1,i]] + [Un[j,i]] + [0])
                    sub_b = sum(uy_three*b_b_stencil)
                    b[i] = sub_b

        return b



    U = U0.copy()
    U_record = []
    for species_index in range(n_species):
        U_record.append(np.zeros([J, I, T])) #DO NOT SIMPLIFY TO U_record = [np.zeros([J, I, T])]*n_species

    A_inv = [np.linalg.inv(a) for a in A_list]
    # unittime=0

    # @numba.jit(nopython=True)
    def adi_forloop(U,N,A_inv):

        # print('inside forloop')
        for ti in range(N):
            #First step: solve in y direction from n -> n+1/2
            U_half = U.copy()
            # f0 = f.dudt(U)
            f0 = f(U, ti, (2,3))
            for i in range(I):
                for n in range(n_species):
                    # U_half[n][:,i] = solve_banded((1, 1), ab_list[n], b('y',i,alpha[n],U[n]) +  f0[n][:,i]*(dt/2)) #CN step in one dimension to get banded(tridiagonal) A matrix
                    U_half[n][:,i] = A_inv[n].dot(b('y',i,alpha[n],U[n])+  f0[n][:,i]*(dt/2)) # Dot product with inverse rather than solve system of equations
                    # U_half[n][:,i] = A_inv[n].dot(1+  f0[n][:,i]*(dt/2)) # Dot product with inverse rather than solve system of equations

            #Second step: solve in x direction from n+1/2 -> n+1
            U_new = U_half.copy()
            # f1 = f.dudt(U_half)
            f1 = f(U_half, ti, (2,3))
            for j in range(J):
                for n in range(n_species):
                    # U_new[n][j,:] = solve_banded((1, 1), ab_list[n], b('x',j,alpha[n],U_half[n]) + f1[n][j,:]*(dt/2))
                    U_new[n][j,:] = A_inv[n].dot(b('x',j,alpha[n],U_half[n])+  f1[n][j,:]*(dt/2)) # Dot product with inverse rather than solve system of equations
                    # U_new[n][j,:] = A_inv[n].dot(1+  f1[n][j,:]*(dt/2)) # Dot product with inverse rather than solve system of equations
                    # check whether dot or * is faster in numba

            # hour = ti / (N / T)
            # # print(np.amax(np.abs(D[0]*laplace(U[0]) + f.dudt(U)[0]))) 
            # if hour % 10 == 0:  #only consider division at unit time (hour)
            #     #append results into top_array for records
            #     for species_index in range(n_species):
            #         U_record[species_index][:, :, int(hour)] = U_new[species_index] #issue in this line
            # if hour % 50 == 0:   
            #     uprime_list=[np.amax(np.abs(U_new[n]-U[n])) for n in range(n_species)]
            #     # print(np.amax(uprime_list))

            #     if all(i <= 10e-4 for i in uprime_list):
            #         if hour > 200:
            #             # print('Error1')
            #             # print(np.amax(np.abs(U_new[-2]-U[-2])))
            #             # print('Error2')
            #             # print([np.amax(np.abs(D[i]*laplace(U[i]) + f.dudt(U)[i])) for i in range(n_species)])
            #             print('converged')
            #             return U_record, U
            #             break


            U = U_new.copy()
        
        return U_record, U
    U_record, U = adi_forloop(U,N,A_inv)
        #check if solution is correct
        # print(np.amax(D[-1]*laplace(U) + f.dudt(U)[-1]))
    return U_record, U






n_species=2


#solver parameters
L_x = 5; L_y=L_x; dx =0.1; dy=dx;J = int(L_x/dx); I=J
T =25; dt = 0.005; N = int(T/dt)
T =25; dt = 0.01; N = int(T/dt)


suggesteddt = float(dx*dx*2)
print(dt, suggesteddt)
x_grid = numpy.array([j*dx for j in range(J)])
y_grid = numpy.array([i*dy for i in range(I)])
t_grid = numpy.array([n*dt for n in range(N)])

D=[0.02,0.4]





U_record, U = adi(L_x,L_y,J,I,T,N, n_species, D)


plt.imshow(-U[0],cmap='PiYG')
plt.savefig('Turing.pdf')
plt.show()
