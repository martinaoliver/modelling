#hello new
import time
start_time = time.time()
import numpy as np
import matplotlib.pyplot as plt
np.set_printoptions(precision=3)
from class_circuit_eq import circuit1_eq
from class_circuit_eq import circuit2_eq
import pickle
from tqdm import tqdm
from matplotlib.colors import LinearSegmentedColormap
import datetime
from PIL import Image
from PIL import Image, ImageFilter, ImageChops
from scipy.optimize import curve_fit
from scipy.sparse.linalg import spsolve
from scipy.linalg import solve_banded
from scipy.sparse import csr_matrix
from cProfile import Profile
from numba import jit, prange
import os
#############
###SOLVERS###
#############

def crank_nicolson(par_dict,steadystates,L,J,T,N,circuit_n, boundary_coef = 1, perturbation = 0.001, sinusoidal = False):

    parent_list = [circuit1_eq, circuit2_eq]
    eq = parent_list[circuit_n-1](par_dict)

    d_A = eq.d_A
    d_B = eq.d_B


    A = steadystates[0]
    B = steadystates[1]
    C = steadystates[2]
    D = steadystates[3]
    E = steadystates[4]
    F = steadystates[5]

    date = datetime.date.today()

    dx = float(L)/float(J-1) # our one unit domain is divided into intervals of lenght dx
    x_grid = np.array([j*dx for j in range(J)]) #array with positions in our grid

    dt = float(T)/float(N-1) # lenght of the N time intervals found time T.
    t_grid = np.array([n*dt for n in range(N)])

    #sigma present in crank nicolson expression
    sigma_a = float(d_A*dt)/float((2.*dx*dx))
    sigma_b = float(d_B*dt)/float((2.*dx*dx))
    sigma_c = float(0*dt)/float((2.*dx*dx))
    sigma_d = float(0*dt)/float((2.*dx*dx))
    sigma_e = float(0*dt)/float((2.*dx*dx))
    sigma_f = float(0*dt)/float((2.*dx*dx))



    ssA = A
    ssB = B
    ssC = C
    ssD = D
    ssE = E
    ssF = F

    # if sinusoidal == False:
    A = np.random.uniform(low=ssA-ssA*perturbation,high=ssA+ssA*perturbation,size = J)
    B = np.random.uniform(low=ssB-ssB*perturbation,high=ssB+ssB*perturbation,size = J)
    C = np.random.uniform(low=ssC-ssC*perturbation,high=ssC+ssC*perturbation,size = J)
    D = np.random.uniform(low=ssD-ssD*perturbation,high=ssD+ssD*perturbation,size = J)
    E = np.random.uniform(low=ssE-ssE*perturbation,high=ssE+ssE*perturbation,size = J)
    F = np.random.uniform(low=ssF-ssF*perturbation,high=ssF+ssF*perturbation,size = J)

    # def sinusoidal_concentrations(J, n_waves=4):
    #     x = np.linspace(0, 3.14 * 2 * n_waves, J)  # start,stop,n_steps
    #     y = (np.sin(x)) * 0.0001 + 0.001
    #     return y
    #
    # if sinusoidal == True:
    #     A = sinusoidal_concentrations(J)
    #     B = sinusoidal_concentrations(J)
    #     C = sinusoidal_concentrations(J)
    #     D = sinusoidal_concentrations(J)
    #     E = sinusoidal_concentrations(J)
    #     F = sinusoidal_concentrations(J)

    A_a = np.diagflat([-sigma_a for i in range(J-1)], -1) +\
          np.diagflat([1.+(1+boundary_coef)*sigma_a]+[1.+2.*sigma_a for i in range(J-2)]+[1.+(1+boundary_coef)*sigma_a]) +\
          np.diagflat([-sigma_a for i in range(J-1)], 1)

    B_a = np.diagflat([sigma_a for i in range(J-1)], -1) +\
          np.diagflat([1.-(1+boundary_coef)*sigma_a]+[1.-2.*sigma_a for i in range(J-2)]+[1.-(1+boundary_coef)*sigma_a]) +\
          np.diagflat([sigma_a for i in range(J-1)], 1)


    A_b = np.diagflat([-sigma_b for i in range(J-1)], -1) +\
          np.diagflat([1.+(1+boundary_coef)*sigma_b]+[1.+2.*sigma_b for i in range(J-2)]+[1.+(1+boundary_coef)*sigma_b]) +\
          np.diagflat([-sigma_b for i in range(J-1)], 1)

    B_b = np.diagflat([sigma_b for i in range(J-1)], -1) +\
          np.diagflat([1.-(1+boundary_coef)*sigma_b]+[1.-2.*sigma_b for i in range(J-2)]+[1.-(1+boundary_coef)*sigma_b]) +\
          np.diagflat([sigma_b for i in range(J-1)], 1)


    A_c = np.diagflat([-sigma_c for i in range(J-1)], -1) +\
          np.diagflat([1.+(1+boundary_coef)*sigma_c]+[1.+2.*sigma_c for i in range(J-2)]+[1.+(1+boundary_coef)*sigma_c]) +\
          np.diagflat([-sigma_c for i in range(J-1)], 1)

    B_c = np.diagflat([sigma_c for i in range(J-1)], -1) +\
          np.diagflat([1.-(1+boundary_coef)*sigma_c]+[1.-2.*sigma_c for i in range(J-2)]+[1.-(1+boundary_coef)*sigma_c]) +\
          np.diagflat([sigma_c for i in range(J-1)], 1)


    A_d = np.diagflat([-sigma_d for i in range(J-1)], -1) +\
          np.diagflat([1.+(1+boundary_coef)*sigma_d]+[1.+2.*sigma_d for i in range(J-2)]+[1.+(1+boundary_coef)*sigma_d]) +\
          np.diagflat([-sigma_d for i in range(J-1)], 1)

    B_d = np.diagflat([sigma_d for i in range(J-1)], -1) +\
          np.diagflat([1.-(1+boundary_coef)*sigma_d]+[1.-2.*sigma_d for i in range(J-2)]+[1.-(1+boundary_coef)*sigma_d]) +\
          np.diagflat([sigma_d for i in range(J-1)], 1)


    A_e = np.diagflat([-sigma_e for i in range(J-1)], -1) +\
          np.diagflat([1.+(1+boundary_coef)*sigma_e]+[1.+2.*sigma_e for i in range(J-2)]+[1.+(1+boundary_coef)*sigma_e]) +\
          np.diagflat([-sigma_e for i in range(J-1)], 1)

    B_e = np.diagflat([sigma_e for i in range(J-1)], -1) +\
          np.diagflat([1.-(1+boundary_coef)*sigma_e]+[1.-2.*sigma_e for i in range(J-2)]+[1.-(1+boundary_coef)*sigma_e]) +\
          np.diagflat([sigma_e for i in range(J-1)], 1)


    A_f = np.diagflat([-sigma_f for i in range(J-1)], -1) +\
          np.diagflat([1.+(1+boundary_coef)*sigma_f]+[1.+2.*sigma_f for i in range(J-2)]+[1.+(1+boundary_coef)*sigma_f]) +\
          np.diagflat([-sigma_f for i in range(J-1)], 1)

    B_f = np.diagflat([sigma_f for i in range(J-1)], -1) +\
          np.diagflat([1.-(1+boundary_coef)*sigma_f]+[1.-2.*sigma_f for i in range(J-2)]+[1.-(1+boundary_coef)*sigma_f]) +\
          np.diagflat([sigma_f for i in range(J-1)], 1)


    A_record = np.zeros([J,N])
    B_record = np.zeros([J,N])
    C_record = np.zeros([J,N])
    D_record = np.zeros([J,N])
    E_record = np.zeros([J,N])
    F_record = np.zeros([J,N])



    for ti in tqdm(range(1,N)):


        A_record[:,ti-1] = A
        B_record[:,ti-1] = B
        C_record[:,ti-1] = C
        D_record[:,ti-1] = D
        E_record[:,ti-1] = E
        F_record[:,ti-1] = F

        A_new = np.linalg.solve(A_a, B_a.dot(A) + eq.dAdt_f(A,B,C,D,E,F)*dt)
        B_new = np.linalg.solve(A_b, B_b.dot(B) + eq.dBdt_f(A,B,C,D,E,F)*dt)
        C_new = np.linalg.solve(A_c, B_c.dot(C) + eq.dCdt_f(A,B,C,D,E,F)*dt)
        D_new = np.linalg.solve(A_d, B_d.dot(D) + eq.dDdt_f(A,B,C,D,E,F)*dt)
        E_new = np.linalg.solve(A_e, B_e.dot(E) + eq.dEdt_f(A,B,C,D,E,F)*dt)
        F_new = np.linalg.solve(A_f, B_f.dot(F) + eq.dFdt_f(A,B,C,D,E,F)*dt)

        A = A_new
        B = B_new
        C = C_new
        D = D_new
        E = E_new
        F = F_new

    grids = (x_grid,t_grid)
    records = (A_record,B_record,C_record,D_record,E_record,F_record)
    final_concentration = (A,B,C,D,E,F)

    return records, final_concentration, grids
def adi(par_dict,steadystates,L_x,L_y,J,I,T,N,circuit_n,boundary_coef=0, perturbation = 0.001):

    parent_list = [circuit1_eq, circuit2_eq]
    eq = parent_list[circuit_n - 1](par_dict)


    d_A = eq.d_A
    d_B = eq.d_B

    ssA = steadystates[0]
    ssB = steadystates[1]
    ssC = steadystates[2]
    ssD = steadystates[3]
    ssE = steadystates[4]
    ssF = steadystates[5]

    date = datetime.date.today()

    #dx needs to be large and dt small
    dx = float(L_x)/float(J-1)
    dy = float(L_y)/float(I-1)
    dt = float(T)/float(N-1)

    x_grid = np.array([j*dx for j in range(J)])
    y_grid = np.array([i*dy for i in range(I)])
    t_grid = np.array([n*dt for n in range(N)])

    # def uniformly(I, J, ss, perturbation):
    #     A = np.random.uniform(low=ssA - perturbation, high=ssA + perturbation, size=(I, J))
    #     return A
    #
    # sss = []
    # A_F = []
    # for i in range(len(sss)):
    #     A_F.append(uniformly(I,J, sss[i], perturbation)

    A = np.random.uniform(low=ssA-perturbation,high=ssA+perturbation,size = (I,J))
    B = np.random.uniform(low=ssB-perturbation,high=ssB+perturbation,size = (I,J))
    C = np.random.uniform(low=ssC-perturbation,high=ssC+perturbation,size = (I,J))
    D = np.random.uniform(low=ssD-perturbation,high=ssD+perturbation,size = (I,J))
    E = np.random.uniform(low=ssE-perturbation,high=ssE+perturbation,size = (I,J))
    F = np.random.uniform(low=ssF-perturbation,high=ssF+perturbation,size = (I,J))

    alpha_x_A = d_A*dt/(2.*dx*dx)
    alpha_y_A = d_A*dt/(2.*dy*dy)

    alpha_x_B = d_B*dt/(2.*dx*dx)
    alpha_y_B = d_B*dt/(2.*dy*dy)

    alpha_x_C = 0*dt/(2.*dx*dx)
    alpha_y_C = 0*dt/(2.*dy*dy)

    alpha_x_D = 0*dt/(2.*dx*dx)
    alpha_y_D = 0*dt/(2.*dy*dy)

    alpha_x_E = 0*dt/(2.*dx*dx)
    alpha_y_E = 0*dt/(2.*dy*dy)

    alpha_x_F = 0*dt/(2.*dx*dx)
    alpha_y_F = 0*dt/(2.*dy*dy)


    A_A = np.diagflat([-alpha_x_A for j in range(J-1)], -1)+\
          np.diagflat([1.+(1+boundary_coef)*alpha_x_A]+[1.+2.*alpha_x_A for j in range(J-2)]+[1.+(1+boundary_coef)*alpha_x_A], 0)+\
          np.diagflat([-alpha_x_A for j in range(J-1)], 1)

    C_A = np.diagflat([-alpha_y_A for i in range(I-1)], -1)+\
          np.diagflat([1.+(1+boundary_coef)*alpha_y_A]+[1.+2.*alpha_y_A for i in range(I-2)]+[1.+(1+boundary_coef)*alpha_y_A], 0)+\
          np.diagflat([-alpha_y_A for i in range(I-1)], 1)

    A_B = np.diagflat([-alpha_x_B for j in range(J-1)], -1)+\
          np.diagflat([1.+(1+boundary_coef)*alpha_x_B]+[1.+2.*alpha_x_B for j in range(J-2)]+[1.+(1+boundary_coef)*alpha_x_B], 0)+\
          np.diagflat([-alpha_x_B for j in range(J-1)], 1)

    C_B = np.diagflat([-alpha_y_B for i in range(I-1)], -1)+\
          np.diagflat([1.+(1+boundary_coef)*alpha_y_B]+[1.+2.*alpha_y_B for i in range(I-2)]+[1.+(1+boundary_coef)*alpha_y_B], 0)+\
          np.diagflat([-alpha_y_B for i in range(I-1)], 1)

    A_C = np.diagflat([-alpha_x_C for j in range(J-1)], -1)+\
          np.diagflat([1.+(1+boundary_coef)*alpha_x_C]+[1.+2.*alpha_x_C for j in range(J-2)]+[1.+(1+boundary_coef)*alpha_x_C], 0)+\
          np.diagflat([-alpha_x_C for j in range(J-1)], 1)

    C_C = np.diagflat([-alpha_y_C for i in range(I-1)], -1)+\
          np.diagflat([1.+(1+boundary_coef)*alpha_y_C]+[1.+2.*alpha_y_C for i in range(I-2)]+[1.+(1+boundary_coef)*alpha_y_C], 0)+\
          np.diagflat([-alpha_y_C for i in range(I-1)], 1)

    A_D = np.diagflat([-alpha_x_D for j in range(J-1)], -1)+\
          np.diagflat([1.+(1+boundary_coef)*alpha_x_D]+[1.+2.*alpha_x_D for j in range(J-2)]+[1.+(1+boundary_coef)*alpha_x_D], 0)+\
          np.diagflat([-alpha_x_D for j in range(J-1)], 1)

    C_D = np.diagflat([-alpha_y_D for i in range(I-1)], -1)+\
          np.diagflat([1.+(1+boundary_coef)*alpha_y_D]+[1.+2.*alpha_y_D for i in range(I-2)]+[1.+(1+boundary_coef)*alpha_y_D], 0)+\
          np.diagflat([-alpha_y_D for i in range(I-1)], 1)

    A_E = np.diagflat([-alpha_x_E for j in range(J-1)], -1)+\
          np.diagflat([1.+(1+boundary_coef)*alpha_x_E]+[1.+2.*alpha_x_E for j in range(J-2)]+[1.+(1+boundary_coef)*alpha_x_E], 0)+\
          np.diagflat([-alpha_x_E for j in range(J-1)], 1)

    C_E = np.diagflat([-alpha_y_E for i in range(I-1)], -1)+\
          np.diagflat([1.+(1+boundary_coef)*alpha_y_E]+[1.+2.*alpha_y_E for i in range(I-2)]+[1.+(1+boundary_coef)*alpha_y_E], 0)+\
          np.diagflat([-alpha_y_E for i in range(I-1)], 1)

    A_F = np.diagflat([-alpha_x_F for j in range(J-1)], -1)+\
          np.diagflat([1.+(1+boundary_coef)*alpha_x_F]+[1.+2.*alpha_x_F for j in range(J-2)]+[1.+(1+boundary_coef)*alpha_x_F], 0)+\
          np.diagflat([-alpha_x_F for j in range(J-1)], 1)

    C_F = np.diagflat([-alpha_y_F for i in range(I-1)], -1)+\
          np.diagflat([1.+(1+boundary_coef)*alpha_y_F]+[1.+2.*alpha_y_F for i in range(I-2)]+[1.+(1+boundary_coef)*alpha_y_F], 0)+\
          np.diagflat([-alpha_y_F for i in range(I-1)], 1)


    b_t_stencil_A = np.array([[(1.-(1+boundary_coef)*alpha_y_A) for j in range(J)],
                                 [alpha_y_A for j in range(J)]])
    b_c_stencil_A = np.array([[alpha_y_A for j in range(J)],
                                 [1.-2.*alpha_y_A for j in range(J)],
                                 [alpha_y_A for j in range(J)]])
    b_b_stencil_A = np.array([[alpha_y_A for j in range(J)],
                                 [(1.-(1+boundary_coef)*alpha_y_A) for j in range(J)]])


    b_t_stencil_B = np.array([[(1.-(1+boundary_coef)*alpha_y_B) for j in range(J)],
                                 [alpha_y_B for j in range(J)]])
    b_c_stencil_B = np.array([[alpha_y_B for j in range(J)],
                                 [1.-2.*alpha_y_B for j in range(J)],
                                 [alpha_y_B for j in range(J)]])
    b_b_stencil_B = np.array([[alpha_y_B for j in range(J)],
                                 [(1.-(1+boundary_coef)*alpha_y_B) for j in range(J)]])


    b_t_stencil_C = np.array([[(1.-(1+boundary_coef)*alpha_y_C) for j in range(J)],
                                 [alpha_y_C for j in range(J)]])
    b_c_stencil_C = np.array([[alpha_y_C for j in range(J)],
                                 [1.-2.*alpha_y_C for j in range(J)],
                                 [alpha_y_C for j in range(J)]])
    b_b_stencil_C = np.array([[alpha_y_C for j in range(J)],
                                 [(1.-(1+boundary_coef)*alpha_y_C) for j in range(J)]])


    b_t_stencil_D = np.array([[(1.-(1+boundary_coef)*alpha_y_D) for j in range(J)],
                                 [alpha_y_D for j in range(J)]])
    b_c_stencil_D = np.array([[alpha_y_D for j in range(J)],
                                 [1.-2.*alpha_y_D for j in range(J)],
                                 [alpha_y_D for j in range(J)]])
    b_b_stencil_D = np.array([[alpha_y_D for j in range(J)],
                                 [(1.-(1+boundary_coef)*alpha_y_D) for j in range(J)]])


    b_t_stencil_E = np.array([[(1.-(1+boundary_coef)*alpha_y_E) for j in range(J)],
                                 [alpha_y_E for j in range(J)]])
    b_c_stencil_E = np.array([[alpha_y_E for j in range(J)],
                                 [1.-2.*alpha_y_E for j in range(J)],
                                 [alpha_y_E for j in range(J)]])
    b_b_stencil_E = np.array([[alpha_y_E for j in range(J)],
                                 [(1.-(1+boundary_coef)*alpha_y_E) for j in range(J)]])


    b_t_stencil_F = np.array([[(1.-(1+boundary_coef)*alpha_y_F) for j in range(J)],
                                 [alpha_y_F for j in range(J)]])
    b_c_stencil_F = np.array([[alpha_y_F for j in range(J)],
                                 [1.-2.*alpha_y_F for j in range(J)],
                                 [alpha_y_F for j in range(J)]])
    b_b_stencil_F = np.array([[alpha_y_F for j in range(J)],
                                 [(1.-(1+boundary_coef)*alpha_y_F) for j in range(J)]])


    d_l_stencil_A = np.array([[1.-(1+boundary_coef)*alpha_x_A, alpha_x_A] for i in range(I)])
    d_c_stencil_A = np.array([[alpha_x_A, 1.-2.*alpha_x_A, alpha_x_A] for i in range(I)])
    d_r_stencil_A = np.array([[alpha_x_A,1.-(1+boundary_coef)*alpha_x_A] for i in range(I)])

    d_l_stencil_B = np.array([[1.-(1+boundary_coef)*alpha_x_B, alpha_x_B] for i in range(I)])
    d_c_stencil_B = np.array([[alpha_x_B, 1.-2.*alpha_x_B, alpha_x_B] for i in range(I)])
    d_r_stencil_B = np.array([[alpha_x_B, 1.-(1+boundary_coef)*alpha_x_B] for i in range(I)])

    d_l_stencil_C = np.array([[1.-(1+boundary_coef)*alpha_x_C, alpha_x_C] for i in range(I)])
    d_c_stencil_C = np.array([[alpha_x_C, 1.-2.*alpha_x_C, alpha_x_C] for i in range(I)])
    d_r_stencil_C = np.array([[alpha_x_C, 1.-(1+boundary_coef)*alpha_x_C] for i in range(I)])

    d_l_stencil_D = np.array([[1.-(1+boundary_coef)*alpha_x_D, alpha_x_D] for i in range(I)])
    d_c_stencil_D = np.array([[alpha_x_D, 1.-2.*alpha_x_D, alpha_x_D] for i in range(I)])
    d_r_stencil_D = np.array([[alpha_x_D,1.-(1+boundary_coef)*alpha_x_D] for i in range(I)])

    d_l_stencil_E = np.array([[1.-(1+boundary_coef)*alpha_x_E, alpha_x_E] for i in range(I)])
    d_c_stencil_E = np.array([[alpha_x_E, 1.-2.*alpha_x_E, alpha_x_E] for i in range(I)])
    d_r_stencil_E = np.array([[alpha_x_E, 1.-(1+boundary_coef)*alpha_x_E] for i in range(I)])

    d_l_stencil_F = np.array([[1.-(1+boundary_coef)*alpha_x_F, alpha_x_F] for i in range(I)])
    d_c_stencil_F = np.array([[alpha_x_F, 1.-2.*alpha_x_F, alpha_x_F] for i in range(I)])
    d_r_stencil_F = np.array([[alpha_x_F, 1.-(1+boundary_coef)*alpha_x_F] for i in range(I)])



    f_curr_A = eq.dAdt_f(A,B,C,D,E,F)
    f_curr_B = eq.dBdt_f(A,B,C,D,E,F)
    f_curr_C = eq.dCdt_f(A,B,C,D,E,F)
    f_curr_D = eq.dDdt_f(A,B,C,D,E,F)
    f_curr_E = eq.dEdt_f(A,B,C,D,E,F)
    f_curr_F = eq.dFdt_f(A,B,C,D,E,F)


    def b_A(i):
        if i <= I-2 and i >= 1:
            A_y = A[[i+1, i, i-1], :]
            return np.sum(A_y*b_c_stencil_A, axis=0)
        elif i == I-1:
            A_y = A[[I-1, I-2], :]
            return np.sum(A_y*b_t_stencil_A,axis=0)
        elif i == 0:
            A_y = A[[1, 0], :]
            return np.sum(A_y*b_b_stencil_A,axis=0)

    def b_B(i):
        if i <= I-2 and i >= 1:
            B_y = B[[i+1, i, i-1], :]
            return np.sum(B_y*b_c_stencil_B,axis=0)
        elif i == I-1:
            B_y = B[[I-1, I-2], :]
            return np.sum(B_y*b_t_stencil_B,axis=0)
        elif i == 0:
            B_y = B[[1, 0], :]
            return np.sum(B_y*b_b_stencil_B, axis=0)

    def b_C(i):
        if i <= I-2 and i >= 1:
            C_y = C[[i+1, i, i-1], :]
            return np.sum(C_y*b_c_stencil_C,axis=0)
        elif i == I-1:
            C_y = C[[I-1, I-2], :]
            return np.sum(C_y*b_t_stencil_C,axis=0)
        elif i == 0:
            C_y = C[[1, 0], :]
            return np.sum(C_y*b_b_stencil_C, axis=0)


    def b_D(i):
        if i <= I-2 and i >= 1:
            D_y = D[[i+1, i, i-1], :]
            return np.sum(D_y*b_c_stencil_D, axis=0)
        elif i == I-1:
            D_y = D[[I-1, I-2], :]
            return np.sum(D_y*b_t_stencil_D,axis=0)
        elif i == 0:
            D_y = D[[1, 0], :]
            return np.sum(D_y*b_b_stencil_D,axis=0)

    def b_E(i):
        if i <= I-2 and i >= 1:
            E_y = E[[i+1, i, i-1], :]
            return np.sum(E_y*b_c_stencil_E,axis=0)
        elif i == I-1:
            E_y = E[[I-1, I-2], :]
            return np.sum(E_y*b_t_stencil_E,axis=0)
        elif i == 0:
            E_y = E[[1, 0], :]
            return np.sum(E_y*b_b_stencil_E, axis=0)

    def b_F(i):
        if i <= I-2 and i >= 1:
            F_y = F[[i+1, i, i-1], :]
            return np.sum(F_y*b_c_stencil_F,axis=0)
        elif i == I-1:
            F_y = F[[I-1, I-2], :]
            return np.sum(F_y*b_t_stencil_F,axis=0)
        elif i == 0:
            F_y = F[[1, 0], :]
            return np.sum(F_y*b_b_stencil_F, axis=0)





    def d_A(j):
        if j <= J-2 and j >= 1:
            A_x = A[:, [j-1, j, j+1]]
            return np.sum(A_x*d_c_stencil_A,axis=1)
        if j == 0:
            A_x = A[:, [0, 1]]
            return np.sum(A_x*d_l_stencil_A, axis=1)
        if j == J-1:
            A_x = A[:, [J-2, J-1]]
            return np.sum(A_x*d_r_stencil_A, axis=1)

    def d_B(j):
        if j <= J-2 and j >= 1:
            B_x = B[:, [j-1, j, j+1]]
            return np.sum(B_x*d_c_stencil_B, axis=1)
        if j == 0:
            B_x = B[:, [0, 1]]
            return np.sum(B_x*d_l_stencil_B,axis=1)
        if j == J-1:
            B_x = B[:, [J-2, J-1]]
            return np.sum(B_x*d_r_stencil_B, axis=1)

    def d_C(j):
        if j <= J-2 and j >= 1:
            C_x = C[:, [j-1, j, j+1]]
            return np.sum(C_x*d_c_stencil_C, axis=1)
        if j == 0:
            C_x = C[:, [0, 1]]
            return np.sum(C_x*d_l_stencil_C,axis=1)
        if j == J-1:
            C_x = C[:, [J-2, J-1]]
            return np.sum(C_x*d_r_stencil_C, axis=1)

    def d_D(j):
        if j <= J-2 and j >= 1:
            D_x = D[:, [j-1, j, j+1]]
            return np.sum(D_x*d_c_stencil_D,axis=1)
        if j == 0:
            D_x = D[:, [0, 1]]
            return np.sum(D_x*d_l_stencil_D, axis=1)
        if j == J-1:
            D_x = D[:, [J-2, J-1]]
            return np.sum(D_x*d_r_stencil_D, axis=1)

    def d_E(j):
        if j <= J-2 and j >= 1:
            E_x = E[:, [j-1, j, j+1]]
            return np.sum(E_x*d_c_stencil_E, axis=1)
        if j == 0:
            E_x = E[:, [0, 1]]
            return np.sum(E_x*d_l_stencil_E,axis=1)
        if j == J-1:
            E_x = E[:, [J-2, J-1]]
            return np.sum(E_x*d_r_stencil_E, axis=1)

    def d_F(j):
        if j <= J-2 and j >= 1:
            F_x = F[:, [j-1, j, j+1]]
            return np.sum(F_x*d_c_stencil_F, axis=1)
        if j == 0:
            F_x = F[:, [0, 1]]
            return np.sum(F_x*d_l_stencil_F,axis=1)
        if j == J-1:
            F_x = F[:, [J-2, J-1]]
            return np.sum(F_x*d_r_stencil_F, axis=1)





    A_new = np.zeros([J,I])
    B_new = np.zeros([J,I])
    C_new = np.zeros([J,I])
    D_new = np.zeros([J,I])
    E_new = np.zeros([J,I])
    F_new = np.zeros([J,I])

    A_record = np.zeros([J,I,T])
    B_record = np.zeros([J,I,T])
    C_record = np.zeros([J,I,T])
    D_record = np.zeros([J,I,T])
    E_record = np.zeros([J,I,T])
    F_record = np.zeros([J,I,T])
    unittime = 0
    for n in tqdm(range(N)):

        if n % (N/T) == 0:
            A_record[:,:,unittime] = A
            B_record[:,:,unittime] = B
            C_record[:,:,unittime] = C
            D_record[:,:,unittime] = D
            E_record[:,:,unittime] = E
            F_record[:,:,unittime] = F
            unittime +=1

        f_curr_A = eq.dAdt_f(A,B,C,D,E,F)
        f_curr_B = eq.dBdt_f(A,B,C,D,E,F)
        f_curr_C = eq.dCdt_f(A,B,C,D,E,F)
        f_curr_D = eq.dDdt_f(A,B,C,D,E,F)
        f_curr_E = eq.dEdt_f(A,B,C,D,E,F)
        f_curr_F = eq.dFdt_f(A,B,C,D,E,F)

        for i in range(I):
            A_new[i, :] = np.linalg.solve(A_A, b_A(i)+ f_curr_A[i,:]*dt/2)
            B_new[i, :] = np.linalg.solve(A_B, b_B(i)+ f_curr_B[i,:]*dt/2)
            C_new[i, :] = np.linalg.solve(A_C, b_C(i)+ f_curr_C[i,:]*dt/2)
            D_new[i, :] = np.linalg.solve(A_D, b_D(i)+ f_curr_D[i,:]*dt/2)
            E_new[i, :] = np.linalg.solve(A_E, b_E(i)+ f_curr_E[i,:]*dt/2)
            F_new[i, :] = np.linalg.solve(A_F, b_F(i)+ f_curr_F[i,:]*dt/2)

        for j in range(J):
            A_new[:, j] = np.linalg.solve(C_A, d_A(j)+ f_curr_A[:,j]*dt/2)
            B_new[:, j] = np.linalg.solve(C_B, d_B(j)+ f_curr_B[:,j]*dt/2)
            C_new[:, j] = np.linalg.solve(C_C, d_C(j)+ f_curr_C[:,j]*dt/2)
            D_new[:, j] = np.linalg.solve(C_D, d_D(j)+ f_curr_D[:,j]*dt/2)
            E_new[:, j] = np.linalg.solve(C_E, d_E(j)+ f_curr_E[:,j]*dt/2)
            F_new[:, j] = np.linalg.solve(C_F, d_F(j)+ f_curr_F[:,j]*dt/2)

        A = A_new
        B = B_new
        C = C_new
        D = D_new
        E = E_new
        F = F_new

    grids = (x_grid,y_grid,t_grid)
    records = (A_record,B_record,C_record,D_record,E_record,F_record)
    final_concentration = (A,B,C,D,E,F)


    return records, final_concentration, grids

def matrix_circle(lenght):
    radius = lenght/2
    y,x = np.ogrid[-radius: radius+1, -radius: radius+1]
    circle = x**2+y**2 <= radius**2
    return circle
def openimage(shape_name, size, filter_size = 3, lower_threshold = 20, upper_threshold = 100):
    image = Image.open('%s.png'%shape_name)

    #remove black or white pixels and turn them blue
    pixels = image.load() # create the pixel map
    for i in range(image.size[0]): # for every pixel:
        for j in range(image.size[1]):
            if pixels[i,j] == (0,0,0,255): # if not black:
                pixels[i,j] = (0,0,70,255) # change to white
            elif pixels[i,j] == (255,255,255,255):
                pixels[i,j] = (0,0,70,255)
    #convert blue pixels into white and the rest into black (detect cells)
    multibands = image.split() #split into RGB (cells almost have no blue).
    blue = multibands[2]
    blue = blue.filter(ImageFilter.GaussianBlur(filter_size)) #smooth: gaussian filter to difuminate (quitar puntitos y agujeritos). higher filter --> smoother
    blue_th = blue.point(lambda i: i > lower_threshold and i < upper_threshold and 255)
    blue_th = ImageChops.invert(blue_th)

    #To retain the shape
    wpercent = (size/float(image.size[0]))
    hsize = int((float(image.size[1])*float(wpercent)))
    img = blue_th.resize((size, hsize))

    # To make it square
    image_array = np.array(img)
    image_array = image_array > 1

    side_len = max(image_array.shape)
    square_array = np.zeros((side_len, side_len))
    square_array[:image_array.shape[0], :image_array.shape[1]] = image_array
    return square_array
class circle:
    def __init__(self, kernel_size):
        self._kernel_size = kernel_size
        self._kernel_radius = (self._kernel_size - 1) // 2

        x, y = np.ogrid[-self._kernel_radius:self._kernel_radius+1, -self._kernel_radius:self._kernel_radius+1]
        self._dist = np.sqrt(x**2 + y**2)

    def circle_matrix(self, radius):
        mask = self._dist - radius
        mask = np.clip(mask, 0, 1, out=mask)
        mask *= -1
        mask += 1
        padded_array = np.zeros(shape = (self._kernel_size,self._kernel_size))
        padded_array[:mask.shape[0],:mask.shape[1]] = mask
        return padded_array
def colony_diameter_evolution(x_gridpoints,t_gridpoints,temperature_factor):
    # ppt_size = [0, 0.1, 0.2, 1, 3.58, 4.12, 4.37, 4.93, 5.28, 5.56, 5.72, 5.77]
    ppt_size = [0.5,1,1.5,2,3.58,4.12,4.37,4.93,5.28,5.56,5.72,5.77]
    ydata = [i*x_gridpoints / 0.75 for i in ppt_size]
    xdata = [0, 10, 20, 40, 65, 77, 102, 110, 124, 135, 149, 158]

    def sigmoid(x, L, x0, k, b):
        y = L / (1 + np.exp(-k * (x - x0))) + b
        return (y)

    p0 = [max(ydata), np.median(xdata), 1, min(ydata)]  # this is an mandatory initial guess
    popt, pcov = curve_fit(sigmoid, xdata, ydata, p0, method='lm')
    popt[0] = popt[0]*temperature_factor

    # x = np.linspace(0,150,150*t_gridpoints)
    x = np.linspace(0,1000,1000*t_gridpoints)
    colony_diameter_evolution = sigmoid(x, *popt)
    fig = plt.figure()
    plt.plot(x[:int(160*t_gridpoints)],colony_diameter_evolution[:int(160*t_gridpoints)]/x_gridpoints)
    return colony_diameter_evolution



def adi_shape(par_dict, steadystates, L_x, L_y, J, I, T, N, circuit_n, shape_name, boundary_coef=1, perturbation=0.001,timeresolution='fast',temperature_factor=1,n_species = 6):

    parent_list = [circuit1_eq, circuit2_eq]
    eq = parent_list[circuit_n - 1](par_dict)

    d_A = eq.d_A
    d_B = eq.d_B
    diffusion_rate_list = [d_A,d_B,0,0,0,0]

    date = datetime.date.today()

    # dx needs to be large and dt small
    dx = float(L_x) / float(J - 1)
    dy = float(L_y) / float(I - 1)
    dt = float(T) / float(N - 1)

    x_grid = np.array([j * dx for j in range(J)])
    y_grid = np.array([i * dy for i in range(I)])
    t_grid = np.array([n * dt for n in range(N)])

    x_gridpoints = int(J/L_x)
    t_gridpoints = int(N/T)

    species_list = []
    for index in range(n_species):
        species_list.append(np.random.uniform(low=steadystates[index] - perturbation, high=steadystates[index] + perturbation, size=(I, J)))

    alpha_x = []
    alpha_y = []
    for index in range(n_species):
        alpha_x.append(diffusion_rate_list[index] * dt / (2. * dx * dx))
        alpha_y.append(diffusion_rate_list[index] * dt / (2. * dy * dy))

    def build_AC_matrix(alpha, space, boundary_coef):

        diag_1 = [0] + [-alpha for j in range(space - 1)]
        diag_2 = [1. + (1 + boundary_coef) * alpha] + [1. + 2. * alpha for j in range(space - 2)] + [
            1. + (1 + boundary_coef) * alpha]
        diag_3 = [-alpha for j in range(space - 1)] + [0]

        AC_matrix = np.array([diag_1, diag_2, diag_3])

        return AC_matrix

    matrix_A_list = []
    matrix_C_list = []

    for m in range(n_species):
        matrix_A_list.append(build_AC_matrix(alpha_x[m], J, boundary_coef))
        matrix_C_list.append(build_AC_matrix(alpha_y[m], I, boundary_coef))

    def build_stencil(n_species):

        b_t_stencil_list = []
        b_c_stencil_list = []
        b_b_stencil_list = []
        d_l_stencil_list = []
        d_c_stencil_list = []
        d_r_stencil_list = []

        for specie_index in range(n_species):
            #top center bottom stencils for species X
            b_t_stencil_list.append( np.array([[(1. - (1 + boundary_coef) * alpha_y[specie_index]) for j in range(J)],
                                      [alpha_y[specie_index] for j in range(J)]]))
            b_c_stencil_list.append(np.array([[alpha_y[specie_index] for j in range(J)],
                                      [1. - 2. * alpha_y[specie_index] for j in range(J)],
                                      [alpha_y[specie_index] for j in range(J)]]))
            b_b_stencil_list.append(np.array([[alpha_y[specie_index] for j in range(J)],
                                      [(1. - (1 + boundary_coef) * alpha_y[specie_index]) for j in range(J)]]))


            #left center right stencils for species X
            d_l_stencil_list.append(np.array([[1. - (1 + boundary_coef) * alpha_x[specie_index], alpha_x[specie_index]] for i in range(I)]))
            d_c_stencil_list.append(np.array([[alpha_x[specie_index], 1. - 2. * alpha_x[specie_index], alpha_x[specie_index]] for i in range(I)]))
            d_r_stencil_list.append(np.array([[alpha_x[specie_index], 1. - (1 + boundary_coef) * alpha_x[specie_index]] for i in range(I)]))

        return b_t_stencil_list,b_c_stencil_list,b_b_stencil_list,d_l_stencil_list,d_c_stencil_list,d_r_stencil_list

    b_t_stencil_list,b_c_stencil_list,b_b_stencil_list,d_l_stencil_list,d_c_stencil_list,d_r_stencil_list = build_stencil(n_species)


    def b_X(species_index,i):
        if i <= I - 2 and i >= 1:
            X_y = species_list[species_index][[i + 1, i, i - 1], :]
            return np.sum(X_y * b_c_stencil_list[species_index], axis=0)
        elif i == I - 1:
            X_y = species_list[species_index][[I - 1, I - 2], :]
            return np.sum(X_y * b_t_stencil_list[species_index], axis=0)
        elif i == 0:
            X_y = species_list[species_index][[1, 0], :]
            return np.sum(X_y * b_b_stencil_list[species_index], axis=0)

    def d_X(species_index,j):
        if j <= J - 2 and j >= 1:
            X_x = species_list[species_index][:, [j - 1, j, j + 1]]
            return np.sum(X_x * d_c_stencil_list[species_index], axis=1)
        if j == 0:
            X_x = species_list[species_index][:, [0, 1]]
            return np.sum(X_x * d_l_stencil_list[species_index], axis=1)
        if j == J - 1:
            X_x = species_list[species_index][:, [J - 2, J - 1]]
            return np.sum(X_x * d_r_stencil_list[species_index], axis=1)

    X_new_list = []
    X_record_list = []
    for species_index in range(n_species):
        X_new_list.append(np.zeros([J, I]))
        if timeresolution != 'slow':
            X_record_list.append(np.zeros([J, I, T]))
        elif timeresolution == 'slow':
            X_record_list.append(np.zeros([J, I, N]))

    A_new,B_new,C_new,D_new,E_new,F_new = X_new_list


    unittime = 0

    if shape_name == 'circle':
        # shape = matrix_circle(len(x_grid) - 1)  # if len(x_grid)-1 as input, the output will be a matrix of len(x_grid)
        shape = circle(len(x_grid)).circle_matrix(int(len(x_grid))/2-1)  # if len(x_grid)-1 as input, the output will be a matrix of len(x_grid)
        print('shape = ' + str(np.shape(shape)))

    elif shape_name == 'growing_colony':
        colony_diameter_list = colony_diameter_evolution(x_gridpoints,t_gridpoints,temperature_factor)
        # shape = matrix_circle_growth(0, J - 1, colony_diameter_list)
        # shape = circle(len(x_grid)).circle_matrix(0)
        shape = circle(len(x_grid)).circle_matrix(1)

    else:
        shape = openimage('shapes/%s' % shape_name, size=J)



    for n in tqdm(range(N)):

        # save results over time. these results can be saved a high time resolutions (save every time point) or low time resolution (save at every time unit).
        if timeresolution != 'slow':
            if n % (N / T) == 0:
                for species_index in range(n_species):
                    X_record_list[species_index][:, :, unittime] = species_list[species_index]
                unittime += 1

        elif timeresolution == 'slow':
            for species_index in range(n_species):
                X_record_list[species_index][:, :, unittime] = species_list[species_index]
            unittime += 1
        # redefine equations with new concentrations
        f_curr_A = eq.dAdt_f(species_list[0], species_list[1], species_list[2], species_list[3], species_list[4], species_list[5])
        f_curr_B = eq.dBdt_f(species_list[0], species_list[1], species_list[2], species_list[3], species_list[4], species_list[5])
        f_curr_C = eq.dCdt_f(species_list[0], species_list[1], species_list[2], species_list[3], species_list[4], species_list[5])
        f_curr_D = eq.dDdt_f(species_list[0], species_list[1], species_list[2], species_list[3], species_list[4], species_list[5])
        f_curr_E = eq.dEdt_f(species_list[0], species_list[1], species_list[2], species_list[3], species_list[4], species_list[5])
        f_curr_F = eq.dFdt_f(species_list[0], species_list[1], species_list[2], species_list[3], species_list[4], species_list[5])
        # recalculate new concentrations

        for i in range(I):
            cells_i = np.where(shape[i] != 0)
            A_new[i, cells_i] = (solve_banded((1, 1), matrix_A_list[0], b_X(0, i) + f_curr_A[i, :] * dt / 2))[cells_i]
            B_new[i, cells_i] = (solve_banded((1, 1), matrix_A_list[1], b_X(1, i) + f_curr_B[i, :] * dt / 2))[cells_i]
            C_new[i, cells_i] = (solve_banded((1, 1), matrix_A_list[2], b_X(2, i) + f_curr_C[i, :] * dt / 2))[cells_i]
            D_new[i, cells_i] = (solve_banded((1, 1), matrix_A_list[3], b_X(3, i) + f_curr_D[i, :] * dt / 2))[cells_i]
            E_new[i, cells_i] = (solve_banded((1, 1), matrix_A_list[4], b_X(4, i) + f_curr_E[i, :] * dt / 2))[cells_i]
            F_new[i, cells_i] = (solve_banded((1, 1), matrix_A_list[5], b_X(5, i) + f_curr_F[i, :] * dt / 2))[cells_i]

        for j in range(J):
            cells_j = np.where(shape[j] != 0)
            A_new[cells_j, j] = (solve_banded((1, 1), matrix_C_list[0], d_X(0, j) + f_curr_A[:, j] * dt / 2))[cells_j]
            B_new[cells_j, j] = (solve_banded((1, 1), matrix_C_list[1], d_X(1, j) + f_curr_B[:, j] * dt / 2))[cells_j]
            C_new[cells_j, j] = (solve_banded((1, 1), matrix_C_list[2], d_X(2, j) + f_curr_C[:, j] * dt / 2))[cells_j]
            D_new[cells_j, j] = (solve_banded((1, 1), matrix_C_list[3], d_X(3, j) + f_curr_D[:, j] * dt / 2))[cells_j]
            E_new[cells_j, j] = (solve_banded((1, 1), matrix_C_list[4], d_X(4, j) + f_curr_E[:, j] * dt / 2))[cells_j]
            F_new[cells_j, j] = (solve_banded((1, 1), matrix_C_list[5], d_X(5, j) + f_curr_F[:, j] * dt / 2))[cells_j]


        X_new_list = [A_new,B_new,C_new,D_new,E_new,F_new]
        for species_index in range(n_species):
            X_new_list[species_index] = np.multiply(X_new_list[species_index],shape)
        A_new, B_new, C_new, D_new, E_new, F_new = X_new_list
        species_list = X_new_list

        if shape_name == 'growing_colony':
        # redefine shape at every timepoint to define non-cell pixels as zero concentration.
        #     if n % (N / T) == 0:
        #
        #         current_T = int(n / t_gridpoints)
        #         radius = int(colony_diameter_list[current_T]/2)
        #         # shape = matrix_circle_growth(current_T, J - 1, colony_diameter_list)
        #         shape = circle(len(x_grid)).circle_matrix(radius)




        #
            radius = colony_diameter_list[n]/2 #radius as float
        # shape = matrix_circle_growth(current_T, J - 1, colony_diameter_list)
            shape = circle(len(x_grid)).circle_matrix(radius)



    grids = (x_grid, y_grid, t_grid)
    records = X_record_list
    final_concentration = X_new_list
    # plt.imshow(A)
    # plt.show()

    return records, final_concentration, grids


#GENERAL FUNCTIONS FOR 1D AND 2D
def stability_test(L,J,T,N):
    dx = float(L)/float(J-1)
    dt = float(T)/float(N-1)
    stability_ratio = dt/(dx**2)

    if stability_ratio<1/2:
        print('Acceptable stability ratio: ' + str(stability_ratio))
    else:
        print('Not acceptable stability ratio: ' + str(stability_ratio))


def t_gridpoints_stability(L,J,T):
    dx = float(L)/float(J-1)
    N = T/(0.49*(dx**2)) + 1
    return int(N/T)

def save_numerical_results(records,grids,filename,path):

    pickle_out = open(path + '/results/simulation/records_%s.pkl' % filename,"wb")
    pickle.dump(records, pickle_out)
    pickle_out.close()

    time_execution = (time.time() - start_time)

    pickle_out = open(path + '/results/simulation/grids_%s.pkl' % filename,"wb")
    pickle.dump(grids, pickle_out)
    pickle_out.close()

def numerical_convergence_1D(results,T,tolerance = 0.0001):
    frame_variation_list = []
    J = len(results[0])
    N = len(results[0][0])

    for n in range(N):
        frame_variation = results[3][int(J/2)][N-1] - results[3][int(J/2)][n]
        frame_variation_list.append(frame_variation)
        if abs(frame_variation) < tolerance:
            break
    convergence_time = int(n/T)
    return int(convergence_time)




#############
###plotting##
#############

#1D PLOTTING
def plot_1D_final_concentration(final_concentration,grids,par_ID):
    A = final_concentration[0]
    B = final_concentration[1]
    C = final_concentration[2]
    D = final_concentration[3]
    E = final_concentration[4]
    F = final_concentration[5]

    x_grid = grids[0]


    fig, (ax1,ax2,ax3,ax4,ax5,ax6) = plt.subplots(6, figsize=(10,10))
    fig.suptitle('Crank Nicolson 1D Numerical solution  par_ID%r'%par_ID)
    ax1.plot(x_grid, A, 'r', label = 'A')
    ax1.set(title = 'A', xlabel = 'x (space units)', ylabel = '[A]', ylim=(0, round(np.amax(A), 3) + 0.15*np.amax(A)))
    ax2.plot(x_grid, B,'g', label = 'B')
    ax2.set(title = 'B', xlabel = 'x (space units)', ylabel = '[B]', ylim=(0, round(np.amax(B), 3) + 0.15*np.amax(B)))
    ax3.plot(x_grid, C, 'b',label = 'C')
    ax3.set(title = 'C', xlabel = 'x (space units)', ylabel = '[C]', ylim=(0, round(np.amax(C), 3) + 0.15*np.amax(C)))
    ax4.plot(x_grid, D, 'c', label = 'D')
    ax4.set(title = 'D', xlabel = 'x (space units)', ylabel = '[D]', ylim=(0, round(np.amax(D), 3) + 0.15*np.amax(D)))
    ax5.plot(x_grid, E,'r', label = 'E')
    ax5.set(title = 'E', xlabel = 'x (space units)', ylabel = '[E]', ylim=(0, round(np.amax(E), 3) + 0.15*np.amax(E)))
    ax6.plot(x_grid, F, 'y',label = 'F')
    ax6.set(title = 'F', xlabel = 'x (space units)', ylabel = '[F]', ylim=(0, round(np.amax(F), 3) + 0.15*np.amax(F)))

    fig.tight_layout()
    fig.subplots_adjust(top = 0.85)



def surfpattern(results,grids,morphogen = 0, circuit_n=1):
    # pickle_in_results = open(working_path + "/results/simulation/CN_cir%r_par%s_lenght%rx%r_time%rx%r.mat" % (circuit_n, par_ID,L,J,T,N),"rb")
    # results = pickle.load(pickle_in_results)
    # pickle_in_grids = open(working_path + "/results/simulation/CN_grids_cir%r_par%s_lenght%rx%r_time%rx%r.mat" % (circuit_n, par_ID,L,J,T,N),"rb")
    # grids = pickle.load(pickle_in_grids)
    results = np.transpose(results[morphogen])
    x_grid = grids[0]
    t_grid = grids[1]
    values = results.reshape(len(t_grid),len(x_grid))
    x, t = np.meshgrid(x_grid, t_grid)
    return x , t , values
def reduce_1D_overtime_concentration(par_ID, L, J, T, N, trim, working_path, circuit_n=1,steadystatenumber = 1):

    pickle_in_results = open(working_path + "/results/simulation/CN_cir%r_par%s_lenght%rx%r_time%rx%r.mat" % (circuit_n, par_ID,L,J,T,N),"rb")
    results = pickle.load(pickle_in_results)
    pickle_in_grids = open(working_path + "/results/simulation/CN_grids_cir%r_par%s_lenght%rx%r_time%rx%r.mat" % (circuit_n, par_ID,L,J,T,N),"rb")
    grids = pickle.load(pickle_in_grids)

    x_grid=grids[0]
    t_grid = grids[1]

    sizesample = int(T)

    minimising = np.linspace(0,N-1-trim,sizesample)


    a_results = results[0][:,trim:]
    b_results = results[1][:,trim:]
    c_results = results[2][:,trim:]
    d_results = results[3][:,trim:]
    e_results = results[4][:,trim:]
    f_results = results[5][:,trim:]

    a_reduced = np.zeros([J,sizesample])
    b_reduced = np.zeros([J,sizesample])
    c_reduced = np.zeros([J,sizesample])
    d_reduced = np.zeros([J,sizesample])
    e_reduced = np.zeros([J,sizesample])
    f_reduced = np.zeros([J,sizesample])

    count = 0
    for i in minimising:
        a_reduced[:,count]=a_results[:,int(i)]
        b_reduced[:,count]=b_results[:,int(i)]
        c_reduced[:,count]=c_results[:,int(i)]
        d_reduced[:,count]=d_results[:,int(i)]
        e_reduced[:,count]=e_results[:,int(i)]
        f_reduced[:,count]=f_results[:,int(i)]
        count+=1


    return a_reduced, b_reduced, c_reduced, d_reduced, e_reduced, f_reduced, x_grid, t_grid

#2D PLOTTING
def plot_2D_final_concentration(final_concentration,grids):
    A = final_concentration[0]
    B = final_concentration[1]
    C = final_concentration[2]
    D = final_concentration[3]
    E = final_concentration[4]
    F = final_concentration[5]

    x_grid = grids[0]
    y_grid = grids[1]


    fig, axs = plt.subplots(2,3,figsize=(7.5,4))

    black_yellow= [(45/255,45/255,45/255),(255/255,255/255,0/255)]
    black_red = [(45/255,45/255,45/255),(255/255,0/255,0/255)]
    black_green = [(45/255,45/255,45/255),(50/255,205/255,50/255)]
    yellow_cmap = LinearSegmentedColormap.from_list('black_yellow',black_yellow)
    red_cmap = LinearSegmentedColormap.from_list('black_red',black_red)
    green_cmap = LinearSegmentedColormap.from_list('black_green',black_green)


    i=0
    im1 = axs[0,0].pcolormesh(x_grid, y_grid, A,vmin = 0, vmax = round(np.amax(A), 3), label = 'A', shading = 'flat',cmap=yellow_cmap)
    im2 = axs[0,1].pcolormesh(x_grid, y_grid, B,vmin = 0, vmax = round(np.amax(B), 3), label = 'B', shading = 'flat',cmap=yellow_cmap)
    im3 = axs[0,2].pcolormesh(x_grid, y_grid, C,vmin = 0, vmax = round(np.amax(C), 3), label = 'C', shading = 'flat',cmap=yellow_cmap)
    im4 = axs[1,0].pcolormesh(x_grid, y_grid, D,vmin = 0, vmax = round(np.amax(D), 3), label = 'D', shading = 'flat',cmap=green_cmap)
    im5 = axs[1,1].pcolormesh(x_grid, y_grid, E,vmin = 0, vmax = round(np.amax(E), 3), label = 'E', shading = 'flat',cmap=red_cmap)
    im6 = axs[1,2].pcolormesh(x_grid, y_grid, F,vmin = 0, vmax = round(np.amax(F), 3), label = 'F', shading = 'flat',cmap=green_cmap)





    for ax in axs.flat:
        ax.label_outer()

    count1=0
    morphogens = ('A','B','C','D','E','F')
    ims = (im1,im2,im3,im4,im5,im6)
    for ax in axs.flat:
        ax.set(title=morphogens[count1])
        fig.colorbar(ims[count1], ax=ax)

        count1+=1

    fig.tight_layout()


    plt.show()
    plt.close()



def reduce_2D_overtime_concentration(par_ID, L, I, J, T, N, sizesample, working_path, trim=0, circuit_n=1,steadystatenumber = 1):

    pickle_in_results = open(working_path + "/results/simulation/ADI_cir%r_par%s_lenght%rx%r_time%rx%r.mat" % (circuit_n, par_ID,L,J,T,N),"rb")
    results = pickle.load(pickle_in_results)
    pickle_in_grids = open(working_path + "/results/simulation/ADI_grids_cir%r_par%s_lenght%rx%r_time%rx%r.mat" % (circuit_n, par_ID,L,J,T,N),"rb")
    grids = pickle.load(pickle_in_grids)

    x_grid=grids[0]
    y_grid=grids[1]

    sizesample = int(T)


    minimising = np.linspace(0,T-1-trim,sizesample)

    a_results = results[0][:,:,trim:]
    b_results = results[1][:,:,trim:]
    c_results = results[2][:,:,trim:]
    d_results = results[3][:,:,trim:]
    e_results = results[4][:,:,trim:]
    f_results = results[5][:,:,trim:]

    a_reduced = np.zeros([I,J,sizesample])
    b_reduced = np.zeros([I,J,sizesample])
    c_reduced = np.zeros([I,J,sizesample])
    d_reduced = np.zeros([I,J,sizesample])
    e_reduced = np.zeros([I,J,sizesample])
    f_reduced = np.zeros([I,J,sizesample])


    count = 0
    for i in minimising:
        a_reduced[:,:,count]=a_results[:,:,int(i)]
        b_reduced[:,:,count]=b_results[:,:,int(i)]
        c_reduced[:,:,count]=c_results[:,:,int(i)]
        d_reduced[:,:,count]=d_results[:,:,int(i)]
        e_reduced[:,:,count]=e_results[:,:,int(i)]
        f_reduced[:,:,count]=f_results[:,:,int(i)]
        count+=1

    return a_reduced, b_reduced, c_reduced, d_reduced, e_reduced, f_reduced, x_grid, y_grid

#Normalising 2D intensity matrices to 0-255 values to plot as rgb, ploting reg/green contrast
# def value_rgb_normalisation(OldValue,matrix):
#     OldMin = np.min(matrix)
#     OldMax = np.amax(matrix)
#     NewMin = 0
#     NewMax = 255
#     OldRange = (OldMax - OldMin)
#     NewRange = (NewMax - NewMin)
#     NewValue = int((((OldValue - OldMin) * NewRange) / OldRange) + NewMin)
#     return NewValue
#
# def matrix_rgb_normalisation(matrix):
#     row_n = 0
#     NewMatrix = np.zeros(matrix.shape)
#     for row in matrix:
#         column_n = 0
#         for value in row:
#             NewMatrix[column_n,row_n] = value_rgb_normalisation(value,matrix)
#
#             column_n +=1
#         row_n +=1
#     return NewMatrix

def matrix_rgb_normalisation(matrix):
    row_n = 0
    NewMatrix = np.zeros(matrix.shape)

    OldMin = np.min(matrix)
    OldMax = np.amax(matrix)
    NewMin = 0
    NewMax = 255
    OldRange = (OldMax - OldMin)
    NewRange = (NewMax - NewMin)

    for row in matrix:
        column_n = 0
        for value in row:
            NewMatrix[column_n,row_n] = int((((value - OldMin) * NewRange) / OldRange) + NewMin)
            column_n +=1
        row_n +=1
    return NewMatrix

def plot_redgreen_contrast(final_concentration, scale_factor=10):
    green = final_concentration[5]
    red = final_concentration[4]
    normalised_green = matrix_rgb_normalisation(green)
    normalised_red = matrix_rgb_normalisation(red)
    zeros = np.zeros(normalised_green.shape)
    rgb = np.dstack((normalised_red,normalised_green,zeros))
    rgb = np.rot90(rgb)
    plt.imshow(rgb.astype('uint8'), origin= 'lower')
    tick_positions = np.arange(0, len(normalised_green), len(normalised_green)/4)
    tick_labels = np.arange(0,len(normalised_green)/scale_factor,len(normalised_green)/scale_factor/4 ).round(decimals=2)
    plt.xticks(tick_positions,tick_labels)
    plt.yticks(tick_positions,tick_labels)


def redgreen_contrast_timeseries(records):
    rgb_timeseries = []
    simulation_time = len(records[0][0][0])
    for time in range (simulation_time):
        red_timeseries,green_timeseries = records[4],records[5]
        red = red_timeseries[:,:,time]
        green = green_timeseries[:,:,time]
        normalised_red = matrix_rgb_normalisation(red)
        normalised_green = matrix_rgb_normalisation(green)
        zeros = np.zeros(red.shape)
        rgb = np.dstack((normalised_red,normalised_green,zeros))
        rgb = np.rot90(rgb)
        rgb_timeseries.append(rgb)
    return rgb_timeseries

def save_numerical_results(records,grids,filename,path):

    pickle_out = open(path + '/results/simulation/records_%s.pkl' % filename,"wb")
    pickle.dump(records, pickle_out)
    pickle_out.close()