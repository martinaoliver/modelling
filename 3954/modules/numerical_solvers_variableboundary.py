#hello new
import time
start_time = time.time()
import numpy as np
np.set_printoptions(precision=3)
from class_circuit_eq import *
import pickle
from tqdm import tqdm
import matplotlib
# matplotlib.use('TKAgg')
import matplotlib.pyplot as plt

from matplotlib.colors import LinearSegmentedColormap
import datetime
from PIL import Image
from PIL import Image, ImageFilter, ImageChops
from scipy.optimize import curve_fit
from scipy.sparse.linalg import spsolve
from scipy.linalg import solve_banded
import matplotlib.animation as animation
from scipy.sparse import csr_matrix
from cProfile import Profile
from numba import jit, prange
import os
import copy
from scipy.misc import derivative

import scipy
#############
###SOLVERS###
#############

def crank_nicolson(par_dict,initial_condition,L,J,T,N,circuit_n, boundary_coef = 1,n_species=6, perturbation = 0.001, growth = True, timeresolution='fast',tqdm_disable=True):

    parent_list = [circuit1_eq, circuit2_eq,circuit3_eq,circuit4_eq,circuit5_eq,circuit6_eq,circuit7_eq,circuit8_eq,circuit9_eq]
    eq = parent_list[circuit_n-1](par_dict)

    d_A = eq.d_A
    d_B = eq.d_B
    diffusion_rate_list = [d_A,d_B,0,0,0,0]

    species_list = []
    for index in range(n_species):
        species_list.append(np.random.uniform(low=initial_condition[index] - perturbation, high=initial_condition[index] + perturbation, size=(J)))

    date = datetime.date.today()

    dx = float(L)/float(J-1) # our one unit domain is divided into intervals of lenght dx
    x_grid = np.array([j*dx for j in range(J)]) #array with positions in our grid

    dt = float(T)/float(N-1) # lenght of the N time intervals found time T.
    t_grid = np.array([n*dt for n in range(N)])

    x_gridpoints = int(J/L)
    t_gridpoints = int(N/T)

    #sigma present in crank nicolson expression
    sigma_list = []
    for index in range(n_species):
        sigma_list.append(float(diffusion_rate_list[index] * dt )/ float(2. * dx * dx))



    def build_matrix_A(sigma,J,boundary_coef):
        A = np.diagflat([-sigma for i in range(J-1)], -1) +\
              np.diagflat([1.+(1+boundary_coef)*sigma]+[1.+2.*sigma for i in range(J-2)]+[1.+(1+boundary_coef)*sigma]) +\
              np.diagflat([-sigma for i in range(J-1)], 1)
        return A

    def build_matrix_B(sigma,J,boundary_coef):
        B = np.diagflat([sigma for i in range(J-1)], -1) +\
              np.diagflat([1.-(1+boundary_coef)*sigma]+[1.-2.*sigma for i in range(J-2)]+[1.-(1+boundary_coef)*sigma]) +\
              np.diagflat([sigma for i in range(J-1)], 1)
        return B

    matrix_A_list = []
    matrix_B_list = []

    for m in range(n_species):
        matrix_A_list.append(build_matrix_A(sigma_list[m], J, boundary_coef))
        matrix_B_list.append(build_matrix_B(sigma_list[m], J, boundary_coef))


    X_record_list = []
    for species_index in range(n_species):
        if timeresolution != 'slow':
            X_record_list.append(np.zeros([J, T]))
        elif timeresolution == 'slow':
            X_record_list.append(np.zeros([J, N]))    # A_record = np.zeros([J,N])


    colony_diameter_list = colony_diameter_evolution(x_gridpoints,t_gridpoints)
    t_reduced_grid = [0]
    unittime = 0
    for ti in tqdm(range(1,N),disable=tqdm_disable):

        if timeresolution != 'slow':
            if ti % (N / T) == 0:
                for species_index in range(n_species):
                    X_record_list[species_index][:,  unittime] = species_list[species_index]
                unittime += 1
                t_reduced_grid.append(unittime)



        elif timeresolution == 'slow':
            for species_index in range(n_species):
                X_record_list[species_index][:, unittime] = species_list[species_index]
            unittime += 1
        if n_species==6:
            A,B,C,D,E,F = species_list

            A_new = np.linalg.solve(matrix_A_list[0], matrix_B_list[0].dot(A) + eq.dAdt_f(A,B,C,D,E,F)*dt)
            B_new = np.linalg.solve(matrix_A_list[1], matrix_B_list[1].dot(B) + eq.dBdt_f(A,B,C,D,E,F)*dt)
            C_new = np.linalg.solve(matrix_A_list[2], matrix_B_list[2].dot(C) + eq.dCdt_f(A,B,C,D,E,F)*dt)
            D_new = np.linalg.solve(matrix_A_list[3], matrix_B_list[3].dot(D) + eq.dDdt_f(A,B,C,D,E,F)*dt)
            E_new = np.linalg.solve(matrix_A_list[4], matrix_B_list[4].dot(E) + eq.dEdt_f(A,B,C,D,E,F)*dt)
            F_new = np.linalg.solve(matrix_A_list[5], matrix_B_list[5].dot(F) + eq.dFdt_f(A,B,C,D,E,F)*dt)

            X_new_list = [A_new, B_new, C_new, D_new, E_new, F_new]
        elif n_species==2:
            A,B = species_list
            A_new = np.linalg.solve(matrix_A_list[0], matrix_B_list[0].dot(A) + eq.dAdt_f(species_list)*dt)
            B_new = np.linalg.solve(matrix_A_list[1], matrix_B_list[1].dot(B) + eq.dBdt_f(species_list)*dt)
            X_new_list = [A_new, B_new]

        if growth == True:
            diameter = int(colony_diameter_list[ti])
            len_pad = int((J - diameter) / 2)
            shape = np.concatenate([np.zeros(len_pad),np.ones(diameter),np.zeros(len_pad)])
            if len(shape)!=len(X_new_list[0]):
                shape = np.concatenate([np.zeros(len_pad), np.ones(diameter), np.zeros(len_pad+1)])
            for species_index in range(n_species):
                X_new_list[species_index] = np.multiply(X_new_list[species_index], shape)


        # A_new, B_new, C_new, D_new, E_new, F_new = X_new_list
        species_list = X_new_list

    grids = (x_grid,t_grid,t_reduced_grid)

    return X_record_list, X_new_list, grids

def crank_nicolson_ca(par_dict,initial_condition,L,J,T,N,circuit_n, boundary_coef = 1,n_species=6, perturbation = 0.001, growth = True, timeresolution='fast',tqdm_disable=True,p_division=1):
    parent_list = [circuit1_eq, circuit2_eq,circuit3_eq,circuit4_eq,circuit5_eq,circuit6_eq,circuit7_eq,circuit8_eq,circuit9_eq]
    eq = parent_list[circuit_n-1](par_dict)

    d_A = eq.d_A
    d_B = eq.d_B
    diffusion_rate_list = [d_A,d_B,0,0,0,0]

    species_list = []
    if growth==True:
        for index in range(n_species):
            grid = np.zeros(shape=(J))
            grid[int(len(grid)/2)]=np.random.uniform(low=initial_condition[index] - perturbation, high=initial_condition[index] + perturbation, size=(1))
            species_list.append(grid)
    elif growth ==False:
        for index in range(n_species):
            species_list.append(np.random.uniform(low=initial_condition[index] - perturbation, high=initial_condition[index] + perturbation, size=(J)))



    date = datetime.date.today()

    dx = float(L)/float(J-1) # our one unit domain is divided into intervals of lenght dx
    x_grid = np.array([j*dx for j in range(J)]) #array with positions in our grid

    dt = float(T)/float(N-1) # lenght of the N time intervals found time T.
    t_grid = np.array([n*dt for n in range(N)])

    x_gridpoints = int(J/L)
    t_gridpoints = int(N/T)

    #sigma present in crank nicolson expression
    sigma_list = []
    for index in range(n_species):
        sigma_list.append(float(diffusion_rate_list[index] * dt )/ float(2. * dx * dx))



    def build_matrix_A(sigma,J,boundary_coef):
        A = np.diagflat([-sigma for i in range(J-1)], -1) +\
              np.diagflat([1.+(1+boundary_coef)*sigma]+[1.+2.*sigma for i in range(J-2)]+[1.+(1+boundary_coef)*sigma]) +\
              np.diagflat([-sigma for i in range(J-1)], 1)
        return A

    def build_matrix_B(sigma,J,boundary_coef):
        B = np.diagflat([sigma for i in range(J-1)], -1) +\
              np.diagflat([1.-(1+boundary_coef)*sigma]+[1.-2.*sigma for i in range(J-2)]+[1.-(1+boundary_coef)*sigma]) +\
              np.diagflat([sigma for i in range(J-1)], 1)
        return B

    matrix_A_list = []
    matrix_B_list = []

    for m in range(n_species):
        matrix_A_list.append(build_matrix_A(sigma_list[m], J, boundary_coef))
        matrix_B_list.append(build_matrix_B(sigma_list[m], J, boundary_coef))


    X_record_list = []
    for species_index in range(n_species):
        if timeresolution != 'slow':
            X_record_list.append(np.zeros([J, T]))
        elif timeresolution == 'slow':
            X_record_list.append(np.zeros([J, N]))    # A_record = np.zeros([J,N])


    colony_diameter_list = colony_diameter_evolution(x_gridpoints,t_gridpoints)
    t_reduced_grid = [0]
    unittime = 0


    def check_neighbours(grid,x_pos): #returns grid with the neighbouring points
        neighbours_grid = [grid[x_pos-1], np.nan, grid[x_pos+1]]
        return neighbours_grid

    def cell_automata_colony(species_list,p_division):
        new_species_list = copy.deepcopy(species_list)
        original_species_list = copy.deepcopy(species_list)
        grid=original_species_list[0]

        for x_pos in np.linspace(1,len(grid)-2,len(grid)-2):
            x_pos = int(x_pos)
            if grid[x_pos]!=0:
                neighbours_grid = check_neighbours(grid,x_pos)
                if 0 in neighbours_grid:
                    cell_division=np.random.choice([1,0],p=[p_division,1-p_division])
                    if cell_division==1:
                        index_nocells=np.where(np.array(neighbours_grid )== 0)[0]
                        index_newcell = np.random.choice(index_nocells)
                        count=0
                        for species,new_species in zip(original_species_list,new_species_list):
                            new_species_list[count][index_newcell+x_pos-1] = species[x_pos]/2
                            new_species_list[count][x_pos] = species[x_pos]/2
                            count+=1


        return new_species_list


    for ti in tqdm(range(1,N),disable=tqdm_disable):
        cell_grid = np.where(species_list[0]!=0, 1, species_list[0])
        # print('cell_Grid')
        # print(cell_grid)
        if timeresolution != 'slow':
            if ti % (N / T) == 0:
                for species_index in range(n_species):
                    X_record_list[species_index][:,  unittime] = species_list[species_index]
                unittime += 1
                t_reduced_grid.append(unittime)
        elif timeresolution == 'slow':
            for species_index in range(n_species):
                X_record_list[species_index][:, unittime] = species_list[species_index]
            unittime += 1
        if n_species==6:
            A,B,C,D,E,F = species_list

            A_new = np.linalg.solve(matrix_A_list[0], matrix_B_list[0].dot(A) + eq.dAdt_f(species_list)*dt)
            B_new = np.linalg.solve(matrix_A_list[1], matrix_B_list[1].dot(B) + eq.dBdt_f(species_list)*dt)
            C_new = np.linalg.solve(matrix_A_list[2], matrix_B_list[2].dot(C) + eq.dCdt_f(species_list)*dt)
            D_new = np.linalg.solve(matrix_A_list[3], matrix_B_list[3].dot(D) + eq.dDdt_f(species_list)*dt)
            E_new = np.linalg.solve(matrix_A_list[4], matrix_B_list[4].dot(E) + eq.dEdt_f(species_list)*dt)
            F_new = np.linalg.solve(matrix_A_list[5], matrix_B_list[5].dot(F) + eq.dFdt_f(species_list)*dt)

            X_new_list = [A_new, B_new, C_new, D_new, E_new, F_new]
        elif n_species==2:
            A,B = species_list
            A_new = np.linalg.solve(matrix_A_list[0], matrix_B_list[0].dot(A) + eq.dAdt_f(species_list)*dt)
            B_new = np.linalg.solve(matrix_A_list[1], matrix_B_list[1].dot(B) + eq.dBdt_f(species_list)*dt)
            X_new_list = [A_new, B_new]
            # print('pre-division')
            # print(X_new_list)
        if growth == True:
            # print('cell definition')
            for species_index in range(n_species):
                X_new_list[species_index] = np.multiply(X_new_list[species_index], cell_grid)
            # print(X_new_list)

            # print('post-division')
            X_new_list = cell_automata_colony(X_new_list,p_division)
            X_new_list
            # print(X_new_list)

        species_list = X_new_list

    grids = (x_grid,t_grid,t_reduced_grid)

    return X_record_list, X_new_list, grids

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
def colony_diameter_evolution(x_gridpoints,t_gridpoints):
    # ppt_size = [0, 0.1, 0.2, 1, 3.58, 4.12, 4.37, 4.93, 5.28, 5.56, 5.72, 5.77]
    ppt_size = [0.5,1,1.5,2,3.58,4.12,4.37,4.93,5.28,5.56,5.72,5.77]
    ydata = [i*x_gridpoints / 0.75 for i in ppt_size]
    xdata = [0, 10, 20, 40, 65, 77, 102, 110, 124, 135, 149, 158]

    def sigmoid(x, L, x0, k, b):
        y = L / (1 + np.exp(-k * (x - x0))) + b
        return (y)

    p0 = [max(ydata), np.median(xdata), 1, min(ydata)]  # this is an mandatory initial guess
    popt, pcov = curve_fit(sigmoid, xdata, ydata, p0, method='lm')

    # x = np.linspace(0,150,150*t_gridpoints)
    x = np.linspace(0,1000,1000*t_gridpoints)
    colony_diameter_evolution = sigmoid(x, *popt)
    # fig = plt.figure()
    # plt.plot(x[:int(160*t_gridpoints)],colony_diameter_evolution[:int(160*t_gridpoints)]/x_gridpoints)
    return colony_diameter_evolution



def adi(par_dict, steadystates, L_x, L_y, J, I, T, N, circuit_n, shape_name, boundary_coef=1, perturbation=0.001,timeresolution='fast',n_species = 6,tqdm_disable=True):

    parent_list = [circuit1_eq, circuit2_eq,circuit3_eq,circuit4_eq,circuit5_eq,circuit6_eq,circuit7_eq,circuit8_eq,circuit9_eq]
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

    # A_new,B_new,C_new,D_new,E_new,F_new = X_new_list


    unittime = 0

    if shape_name == 'circle':
        # shape = matrix_circle(len(x_grid) - 1)  # if len(x_grid)-1 as input, the output will be a matrix of len(x_grid)
        shape = circle(len(x_grid)).circle_matrix(int(len(x_grid))/2-1)  # if len(x_grid)-1 as input, the output will be a matrix of len(x_grid)
        print('shape = ' + str(np.shape(shape)))

    elif shape_name == 'growing_colony':
        colony_diameter_list = colony_diameter_evolution(x_gridpoints,t_gridpoints)

        shape = circle(len(x_grid)).circle_matrix(1)
    elif shape_name == 'no_growth_square':
        shape = np.ones([J, I])
    else:
        shape = openimage('shapes/%s' % shape_name, size=J)

    if n_species==6:
        A_new,B_new,C_new,D_new,E_new,F_new = X_new_list

    if n_species==2:
        A_new,B_new = X_new_list



    for n in tqdm(range(N),disable = tqdm_disable):

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
        if n_species==6:
            f_curr_A = eq.dAdt_f(species_list)
            f_curr_B = eq.dBdt_f(species_list)
            f_curr_C = eq.dCdt_f(species_list)
            f_curr_D = eq.dDdt_f(species_list)
            f_curr_E = eq.dEdt_f(species_list)
            f_curr_F = eq.dFdt_f(species_list)
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
        if n_species==2:
            f_curr_A = eq.dAdt_f(species_list)
            f_curr_B = eq.dBdt_f(species_list)

            for i in range(I):
                cells_i = np.where(shape[i] != 0)
                A_new[i, cells_i] = (solve_banded((1, 1), matrix_A_list[0], b_X(0, i) + f_curr_A[i, :] * dt / 2))[cells_i]
                B_new[i, cells_i] = (solve_banded((1, 1), matrix_A_list[1], b_X(1, i) + f_curr_B[i, :] * dt / 2))[cells_i]

            for j in range(J):
                cells_j = np.where(shape[j] != 0)
                A_new[cells_j, j] = (solve_banded((1, 1), matrix_C_list[0], d_X(0, j) + f_curr_A[:, j] * dt / 2))[cells_j]
                B_new[cells_j, j] = (solve_banded((1, 1), matrix_C_list[1], d_X(1, j) + f_curr_B[:, j] * dt / 2))[cells_j]

            X_new_list = [A_new,B_new]

        for species_index in range(n_species):
            X_new_list[species_index] = np.multiply(X_new_list[species_index],shape)
        if n_species==6:
            A_new,B_new,C_new,D_new,E_new,F_new = X_new_list

        if n_species==2:
            A_new,B_new = X_new_list
        species_list = X_new_list

        if shape_name == 'growing_colony':
            radius = colony_diameter_list[n]/2 #radius as float
            shape = circle(len(x_grid)).circle_matrix(radius)



    grids = (x_grid, y_grid, t_grid)
    records = X_record_list
    final_concentration = X_new_list


    return records, final_concentration, grids

def adi_ca(par_dict, initial_condition, L_x, L_y, J, I, T, N, circuit_n, boundary_coef=1, perturbation=0.001,timeresolution='fast',n_species = 6,tqdm_disable=True,p_division=0.41,seed=1):
    np.random.seed(seed)

    parent_list = [circuit1_eq, circuit2_eq,circuit3_eq,circuit4_eq,circuit5_eq,circuit6_eq,circuit7_eq,circuit8_eq,circuit9_eq]
    eq = parent_list[circuit_n - 1](par_dict)
    d_A = eq.d_A
    d_B = eq.d_B


    diffusion_rate_list = [d_A,d_B,0,0,0,0]

    if circuit_n==8:
        diffusion_rate_list = [d_A,0,0,0]
    if circuit_n==9:
        diffusion_rate_list = [d_A,0,0,0,0]

    date = datetime.date.today()

    # dx needs to be large and dt small
    dx = float(L_x) / float(J - 1) #size of cell
    dy = float(L_y) / float(I - 1)
    dt = float(T) / float(N - 1) #size of timestep

    x_grid = np.array([j * dx for j in range(J)])
    y_grid = np.array([i * dy for i in range(I)])
    t_grid = np.array([n * dt for n in range(N)])

    x_gridpoints = int(J/L_x) #x_gridpoints per space unit
    t_gridpoints = int(N/T) #t_gridpoints per timeunit

    species_list = []
    # if growth ==False:
    #     for index in range(n_species):
    #         species_list.append(np.random.uniform(low=steadystates[index] - perturbation, high=steadystates[index] + perturbation, size=(I, J)))
    for index in range(n_species):
        grid = np.zeros(shape=(I,J))
        grid[int(len(grid)/2),int(len(grid)/2)]=np.random.uniform(low=initial_condition[index] - perturbation, high=initial_condition[index] + perturbation, size=(1))
        species_list.append(grid)
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

    unittime = 0

    if n_species==6:
        A_new,B_new,C_new,D_new,E_new,F_new = X_new_list

    elif n_species==2:
        A_new,B_new = X_new_list

    elif circuit_n==8: #nodeA dele
        B_new,D_new,E_new,F_new  = X_new_list

    elif circuit_n==9: #nodeB dele (top cassette)
        A_new,C_new,D_new,E_new,F_new  = X_new_list
    def check_neighbours(grid,x_pos,y_pos): #returns grid with the neighbouring points
        top_array = [grid[y_pos-1, x_pos-1], grid[y_pos-1,x_pos], grid[y_pos-1,x_pos+1]]
        middle_array = [grid[y_pos, x_pos-1], np.nan, grid[y_pos,x_pos+1]]
        bottom_array = [grid[y_pos+1, x_pos-1], grid[y_pos+1,x_pos], grid[y_pos+1,x_pos+1]]
        neighbours_grid = np.array([top_array,middle_array,bottom_array])
        return neighbours_grid

    def cell_automata_colony(species_list,p_division):
        new_species_list = copy.deepcopy(species_list)
        original_species_list = copy.deepcopy(species_list)
        grid=original_species_list[0]

        for x_pos in np.linspace(1,len(grid)-2,len(grid)-2):
            for y_pos in np.linspace(1,len(grid)-2,len(grid)-2):
                x_pos = int(x_pos)
                y_pos = int(y_pos)
                if grid[x_pos,y_pos]!=0:
                    neighbours_grid = check_neighbours(grid,x_pos,y_pos)
                    if 0 in neighbours_grid:
                        cell_division=np.random.choice([1,0],p=[p_division,1-p_division])
                        if cell_division==1:
                            index_nocells=np.where(np.array(neighbours_grid )== 0)

                            index_newcell_y, index_newcell_x = [np.random.choice(index_nocells[0]),np.random.choice(index_nocells[1])]
                            count=0
                            for species,new_species in zip(original_species_list,new_species_list):
                                new_species_list[count][index_newcell_y+y_pos-1,index_newcell_x+x_pos-1] = species[y_pos,x_pos]/2
                                new_species_list[count][y_pos,x_pos] = species[y_pos,x_pos]/2
                                count+=1


        return new_species_list


    cell_grid = np.where(species_list[0]!=0, 1, species_list[0])

    for n in tqdm(range(N),disable = tqdm_disable):

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
        if n_species==6:
            f_curr_A = eq.dAdt_f(species_list)
            f_curr_B = eq.dBdt_f(species_list)
            f_curr_C = eq.dCdt_f(species_list)
            f_curr_D = eq.dDdt_f(species_list)
            f_curr_E = eq.dEdt_f(species_list)
            f_curr_F = eq.dFdt_f(species_list)

            # recalculate new concentrations

            for i in range(I):
                cells_i = np.where(cell_grid[i] != 0)

                A_new[i, cells_i] = (solve_banded((1, 1), matrix_A_list[0], b_X(0, i) + f_curr_A[i, :] * dt / 2))[cells_i]
                B_new[i, cells_i] = (solve_banded((1, 1), matrix_A_list[1], b_X(1, i) + f_curr_B[i, :] * dt / 2))[cells_i]
                C_new[i, cells_i] = (solve_banded((1, 1), matrix_A_list[2], b_X(2, i) + f_curr_C[i, :] * dt / 2))[cells_i]
                D_new[i, cells_i] = (solve_banded((1, 1), matrix_A_list[3], b_X(3, i) + f_curr_D[i, :] * dt / 2))[cells_i]
                E_new[i, cells_i] = (solve_banded((1, 1), matrix_A_list[4], b_X(4, i) + f_curr_E[i, :] * dt / 2))[cells_i]
                F_new[i, cells_i] = (solve_banded((1, 1), matrix_A_list[5], b_X(5, i) + f_curr_F[i, :] * dt / 2))[cells_i]

            for j in range(J):
                cells_j = np.where(cell_grid[j] != 0)
                A_new[cells_j, j] = (solve_banded((1, 1), matrix_C_list[0], d_X(0, j) + f_curr_A[:, j] * dt / 2))[cells_j]
                B_new[cells_j, j] = (solve_banded((1, 1), matrix_C_list[1], d_X(1, j) + f_curr_B[:, j] * dt / 2))[cells_j]
                C_new[cells_j, j] = (solve_banded((1, 1), matrix_C_list[2], d_X(2, j) + f_curr_C[:, j] * dt / 2))[cells_j]
                D_new[cells_j, j] = (solve_banded((1, 1), matrix_C_list[3], d_X(3, j) + f_curr_D[:, j] * dt / 2))[cells_j]
                E_new[cells_j, j] = (solve_banded((1, 1), matrix_C_list[4], d_X(4, j) + f_curr_E[:, j] * dt / 2))[cells_j]
                F_new[cells_j, j] = (solve_banded((1, 1), matrix_C_list[5], d_X(5, j) + f_curr_F[:, j] * dt / 2))[cells_j]


            X_new_list = [A_new,B_new,C_new,D_new,E_new,F_new]


        elif n_species==2:
            f_curr_A = eq.dAdt_f(species_list)
            f_curr_B = eq.dBdt_f(species_list)

            for i in range(I):
                # cells_i = np.where(cell_grid[i] != 0)
                # A_new[i, cells_i] = (solve_banded((1, 1), matrix_A_list[0], b_X(0, i) + f_curr_A[i, :] * dt / 2))[cells_i]
                A_new[i, :] = (solve_banded((1, 1), matrix_A_list[0], b_X(0, i) + f_curr_A[i, :] * dt / 2))
                B_new[i, :] = (solve_banded((1, 1), matrix_A_list[1], b_X(1, i) + f_curr_B[i, :] * dt / 2))

            for j in range(J):
                # cells_j = np.where(cell_grid[j] != 0)
                A_new[:, j] = (solve_banded((1, 1), matrix_C_list[0], d_X(0, j) + f_curr_A[:, j] * dt / 2))
                B_new[:, j] = (solve_banded((1, 1), matrix_C_list[1], d_X(1, j) + f_curr_B[:, j] * dt / 2))

            X_new_list = [A_new,B_new]

        elif circuit_n==8:
            f_curr_B = eq.dBdt_f(species_list)
            f_curr_D = eq.dDdt_f(species_list)
            f_curr_E = eq.dEdt_f(species_list)
            f_curr_F = eq.dFdt_f(species_list)

            # recalculate new concentrations

            for i in range(I):
                cells_i = np.where(cell_grid[i] != 0)

                B_new[i, cells_i] = (solve_banded((1, 1), matrix_A_list[0], b_X(0, i) + f_curr_B[i, :] * dt / 2))[cells_i]
                D_new[i, cells_i] = (solve_banded((1, 1), matrix_A_list[1], b_X(1, i) + f_curr_D[i, :] * dt / 2))[cells_i]
                E_new[i, cells_i] = (solve_banded((1, 1), matrix_A_list[2], b_X(2, i) + f_curr_E[i, :] * dt / 2))[cells_i]
                F_new[i, cells_i] = (solve_banded((1, 1), matrix_A_list[3], b_X(3, i) + f_curr_F[i, :] * dt / 2))[cells_i]

            for j in range(J):
                cells_j = np.where(cell_grid[j] != 0)
                B_new[cells_j, j] = (solve_banded((1, 1), matrix_C_list[0], d_X(0, j) + f_curr_B[:, j] * dt / 2))[cells_j]
                D_new[cells_j, j] = (solve_banded((1, 1), matrix_C_list[1], d_X(1, j) + f_curr_D[:, j] * dt / 2))[cells_j]
                E_new[cells_j, j] = (solve_banded((1, 1), matrix_C_list[2], d_X(2, j) + f_curr_E[:, j] * dt / 2))[cells_j]
                F_new[cells_j, j] = (solve_banded((1, 1), matrix_C_list[3], d_X(3, j) + f_curr_F[:, j] * dt / 2))[cells_j]


            X_new_list = [B_new,D_new,E_new,F_new]

        elif circuit_n==9: #nodeB dele (top cassette)
            f_curr_A = eq.dAdt_f(species_list)
            f_curr_C = eq.dCdt_f(species_list)
            f_curr_D = eq.dDdt_f(species_list)
            f_curr_E = eq.dEdt_f(species_list)
            f_curr_F = eq.dFdt_f(species_list)

            # recalculate new concentrations

            for i in range(I):
                cells_i = np.where(cell_grid[i] != 0)

                A_new[i, cells_i] = (solve_banded((1, 1), matrix_A_list[0], b_X(0, i) + f_curr_A[i, :] * dt / 2))[cells_i]
                C_new[i, cells_i] = (solve_banded((1, 1), matrix_A_list[1], b_X(1, i) + f_curr_C[i, :] * dt / 2))[cells_i]
                D_new[i, cells_i] = (solve_banded((1, 1), matrix_A_list[2], b_X(2, i) + f_curr_D[i, :] * dt / 2))[cells_i]
                E_new[i, cells_i] = (solve_banded((1, 1), matrix_A_list[3], b_X(3, i) + f_curr_E[i, :] * dt / 2))[cells_i]
                F_new[i, cells_i] = (solve_banded((1, 1), matrix_A_list[4], b_X(4, i) + f_curr_F[i, :] * dt / 2))[cells_i]

            for j in range(J):
                cells_j = np.where(cell_grid[j] != 0)
                A_new[cells_j, j] = (solve_banded((1, 1), matrix_C_list[0], d_X(0, j) + f_curr_A[:, j] * dt / 2))[cells_j]
                C_new[cells_j, j] = (solve_banded((1, 1), matrix_C_list[1], d_X(1, j) + f_curr_C[:, j] * dt / 2))[cells_j]
                D_new[cells_j, j] = (solve_banded((1, 1), matrix_C_list[2], d_X(2, j) + f_curr_D[:, j] * dt / 2))[cells_j]
                E_new[cells_j, j] = (solve_banded((1, 1), matrix_C_list[3], d_X(3, j) + f_curr_E[:, j] * dt / 2))[cells_j]
                F_new[cells_j, j] = (solve_banded((1, 1), matrix_C_list[4], d_X(4, j) + f_curr_F[:, j] * dt / 2))[cells_j]


            X_new_list = [A_new,C_new,D_new,E_new,F_new]
        for species_index in range(n_species):
            X_new_list[species_index] = np.multiply(X_new_list[species_index],cell_grid)
        if n_species==6:
            A_new,B_new,C_new,D_new,E_new,F_new = X_new_list

        if n_species==2:
            A_new,B_new = X_new_list
        if circuit_n==8:
            B_new,D_new,E_new,F_new = X_new_list
        if circuit_n==9:
            A_new,C_new,D_new,E_new,F_new = X_new_list


        hour = n / (N / T)
        if hour % 1 == 0:  #only consider division at unit time (hour)
        # if n % (N / T) == 2: #only consider division at unit time (hour)
            # p_division = derivative_growth[n]*2*x_gridpoints/division_index
            # print(p_division)
            if hour>120:
                p_division=0.01
            X_new_list = cell_automata_colony(X_new_list,p_division)
            cell_grid = np.where(X_new_list[0]!=0, 1, X_new_list[0])

        species_list = X_new_list



    grids = (x_grid, y_grid, t_grid)
    records = X_record_list
    final_concentration = X_new_list


    return records, final_concentration, grids



# GENERAL FUNCTIONS FOR 1D AND 2D
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


def plot_1D_final_concentration(final_concentration,grids,mechanism,shape,filename,n_species, parID,save_figure=False,path=''):
    fig, ax = plt.subplots(n_species, figsize=(4,4))
    x_grid = grids[0]
    labels=['A','B','C','D','E','F']
    for n in range(n_species):
        x_grid = grids[0]

        specie=final_concentration[n]
        ax[n].plot(x_grid, specie, 'k', label = labels[n])
        ax[n].set_xticklabels('')


    fig.subplots_adjust(top = 0.85)
    fig.tight_layout()

    if save_figure == False:
        plt.show()
        plt.close()

    if save_figure == True:
        # plt.savefig(path + '/6eq/numerical_confocal/results/figures/1D/%s/%s/1D_%s.jpeg' % (mechanism,shape,filename))
        plt.savefig(path + '1D_%s.jpeg' % (mechanism,shape,filename))
        plt.close()


def surfpattern(results,grids,par_dict,morphogen = 0):
    results = np.transpose(results[morphogen])
    x_grid = grids[0]
    t_grid = grids[2]
    values = results.reshape(len(t_grid),len(x_grid))
    x, t = np.meshgrid(x_grid, t_grid)
    plt.contourf(x,t,results)
    plt.colorbar()
    plt.title(par_dict,size=8)
    plt.show()
    return x , t , values

#2D PLOTTING

def plot_2D_final_concentration(final_concentration,grids,filename,n_species=6):

    # A,B,C,D,E,F = final_concentration
    x_grid,y_grid,t_grid = grids
    labels=['A','B','C','D','E','F']

    fig, axs = plt.subplots(2,int(n_species/2))#,figsize=(7.5,4))
    ax = axs.flatten()
    black_yellow= [(45/255,45/255,45/255),(255/255,255/255,0/255)]
    black_red = [(45/255,45/255,45/255),(255/255,0/255,0/255)]
    black_green = [(45/255,45/255,45/255),(50/255,205/255,50/255)]
    yellow_cmap = LinearSegmentedColormap.from_list('black_yellow',black_yellow)
    red_cmap = LinearSegmentedColormap.from_list('black_red',black_red)
    green_cmap = LinearSegmentedColormap.from_list('black_green',black_green)
    cmap_list = [yellow_cmap,yellow_cmap,yellow_cmap,green_cmap,red_cmap,green_cmap]
    ims = [0]*n_species
    print(n_species)
    for n in range(n_species):
        specie = final_concentration[n]
        ims[n] = ax[n].pcolormesh(x_grid, y_grid, specie, shading='auto', label = labels[n],cmap=cmap_list[n])
    for ax in axs.flat:
        ax.label_outer()
    #
    count1=0
    morphogens = ('A','B','C','D','E','F')
    for ax in axs.flat:
        ax.set(title=morphogens[count1])
        fig.colorbar(ims[count1], ax=ax)

        count1+=1

    fig.tight_layout()
    # plt.savefig(path + '/results/figures/2D/nonturing_numerical/growing_colony/2D_%s.jpeg' % filename)
    # plt.close()
    plt.show()




#Normalising 2D intensity matrices to 0-255 values to plot as rgb, ploting reg/green contrast


def list_rgb_normalisation(matrix):
    row_n = 0
    NewMatrix = np.zeros(matrix.shape)

    OldMin = np.min(matrix)
    OldMax = np.amax(matrix)
    NewMin = 0
    NewMax = 255
    OldRange = (OldMax - OldMin)
    NewRange = (NewMax - NewMin)

    for value in matrix:
        NewMatrix[row_n] = int((((value - OldMin) * NewRange) / OldRange) + NewMin)
        row_n += 1
    return NewMatrix


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
            NewMatrix[column_n, row_n] = int((((value - OldMin) * NewRange) / OldRange) + NewMin)
            column_n += 1
        row_n += 1
    return NewMatrix, OldMin, OldMax


def plot_redgreen_contrast(final_concentration, mm, mechanism, shape, filename, path, parID=0, dimension='2D', scale_factor=10, save_figure=False):
    green = final_concentration[-1]
    red = final_concentration[-2]
    var_red,var_green = [(np.amax(red) - np.amin(red)), (np.amax(green) - np.amin(green))]

    x_grid = np.linspace(0, mm, len(green))
    if dimension == '1D':
        if any(i>1 for i in [var_red,var_green]): #enough variability to normalise
            # normalised_red = list_rgb_normalisation(red)
            # normalised_green = list_rgb_normalisation(green)

            red *= 255.0 / red.max()
            green *= 255.0 / green.max()


        else:
            print('non normalised')
        normalised_red,normalised_green = red,green

        plt.plot(x_grid, normalised_red, 'r', label='E')
        plt.plot(x_grid, normalised_green, 'g', label='D')
        plt.legend()



    elif dimension == '2D':
        normalised_red, redmin, redmax = matrix_rgb_normalisation(red)
        normalised_green, greenmin, greenmax = matrix_rgb_normalisation(green)
        zeros = np.zeros(normalised_green.shape)
        rgb = np.dstack((normalised_red, normalised_green, zeros))
        rgb = np.rot90(rgb)
        if save_figure != 'results':
            plt.imshow(rgb.astype('uint8'), origin='lower')
            tick_positions = np.arange(0, len(normalised_green), len(normalised_green) / 4)
            tick_labels = np.arange(0, len(normalised_green) / scale_factor,
                                    len(normalised_green) / scale_factor / 4).round(decimals=2)
            plt.xticks(tick_positions, tick_labels)
            plt.yticks(tick_positions, tick_labels)
            plt.ylabel('y axis (mm)', size=16)
            plt.xlabel('x axis (mm)', size=16)
            plt.yticks(size=15)
            plt.xticks(size=15)
            plt.title('parID=' + str(parID), size=14)
            np.set_printoptions(precision=2)
            plt.text(1,1,'mCherry = [%r-%r]'%(np.around(redmin,2),np.around(redmax,2)),c='r')
            plt.text(1,5,'GPF = [%r-%r]'%(np.around(greenmin),np.around(greenmax)),c='g')
            plt.tight_layout()

            if save_figure == True:
                plt.savefig(path + '/%s_%s.jpeg' % (dimension, filename),dpi=2000)
                plt.close()
            else:
                plt.show()


        return rgb




def redgreen_contrast_timeseries(records):
    rgb_timeseries = []
    simulation_time = len(records[0][0][0])
    for time in range (simulation_time):
        red_timeseries,green_timeseries = records[-2],records[-1]
        red = red_timeseries[:,:,time]
        green = green_timeseries[:,:,time]
        normalised_red = matrix_rgb_normalisation(red)[0]
        normalised_green = matrix_rgb_normalisation(green)[0]
        print(normalised_red[100,100]); print(normalised_green[100,100])
        zeros = np.zeros(red.shape)
        rgb = np.dstack((normalised_red,normalised_green,zeros))
        rgb = np.rot90(rgb)
        rgb_timeseries.append(rgb)
    return rgb_timeseries


def save_numerical_results(records,grids,filename,path):

    pickle_out = open(path + '/results/simulation/records_%s.pkl' % filename,"wb")
    pickle.dump(records, pickle_out)
    pickle_out.close()

def show_rgbvideo(timeseries_unstacked,parID):
    time=0

    fig = plt.plot()
    rgb_timeseries=timeseries_unstacked # Read the numpy matrix with images in the rows
    im=plt.imshow(rgb_timeseries[0].astype('uint8'), origin= 'lower')
    for time in range(len(rgb_timeseries)):
        im.set_data(rgb_timeseries[time].astype('uint8'))
        plt.title(parID)
        plt.pause(0.01)
    plt.show()

def timeseries_1D_rgb(simulation, L, J,T,filename,path,save_figure=False):

    red = simulation[4]
    green = simulation[5]
    var_red,var_green = [(np.amax(red[:,-2]) - np.amin(red[:,-2])), (np.amax(green[:,-2]) - np.amin(green[:,-2]))]

    # print('var red=' + str(var_red )+ ', max green=' + str(var_green))
    red *= 255.0/red.max()
    green *= 255.0/green.max()
    i=0
    fig, ax = plt.subplots()

    x = np.linspace(0,L,J)
    line_red, = ax.plot(x, red[:,i],c='r')
    line_green, = ax.plot(x, green[:,i],c='g')
    ax.set_ylim(0,255)


    def animate(i):
        line_red.set_ydata(red[:,i])  # update the data.
        line_green.set_ydata(green[:,i])  # update the data.
        return line_red,line_green,


    ani = animation.FuncAnimation(
        fig, animate, frames=range(T),interval=20, save_count=50)
    # plt.show()
    plt.rcParams['animation.ffmpeg_path'] = '~/Documents/modelling/env1/bin/ffmpeg'
    Writer = animation.writers['pillow']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    if save_figure == True:
        ani.save(path + '//6eq/numerical_confocal/results/videos/1D/general/growing_colony/%s.gif'%filename, writer=writer)
    else:
        plt.show()
