


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
def t_gridpoints_stability(L,J,T):
    dx = float(L)/float(J-1)
    N = T/(0.49*(dx**2)) + 1
    return int(N/T)


def adi_ca(par_dict, initial_condition, L_x, L_y, J, I, T, N, circuit_n, boundary_coef=1, perturbation=0.001,timeresolution='fast',n_species = 9,tqdm_disable=True,p_division=0.41,seed=1):
    np.random.seed(seed)

    # parent_list = [circuit1_eq, circuit2_eq,circuit3_eq,circuit4_eq,circuit5_eq,circuit6_eq,circuit7_eq,circuit8_eq,circuit9_eq]
    eq = circuit11_eq(par_dict)
    d_H = eq.d_H
    diffusion_rate_list = [0,d_H,0,0,0,0,0,0,0]


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

    mB_new,H_new,mD_new,mE_new,mF_new,B_new,D_new,E_new,F_new  = X_new_list

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
        # print(species_list[0])
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
        f_curr_mB,f_curr_H,f_curr_mD,f_curr_mE,f_curr_mF,f_curr_B,f_curr_D,f_curr_E,f_curr_F = [eq.function_list(species_list,wvn=0)[i] for i in range(n_species)]
        # print(f_curr_mE)
        # print(species_list[3])
        # print('afeg')

        for i in range(I):
            cells_i = np.where(cell_grid[i] != 0)
            # print('hrtg')
            # print(cells_i[0])
            mB_new[i, cells_i] = (solve_banded((1, 1), matrix_A_list[0], b_X(0, i) + f_curr_mB[i, :] * dt / 2))[cells_i]
            H_new[i, cells_i] = (solve_banded((1, 1), matrix_A_list[1], b_X(1, i) + f_curr_H[i, :] * dt / 2))[cells_i]
            mD_new[i, cells_i] = (solve_banded((1, 1), matrix_A_list[2], b_X(2, i) + f_curr_mD[i, :] * dt / 2))[cells_i]
            mE_new[i, cells_i] = (solve_banded((1, 1), matrix_A_list[3], b_X(3, i) + f_curr_mE[i, :] * dt / 2))[cells_i]
            mF_new[i, cells_i] = (solve_banded((1, 1), matrix_A_list[4], b_X(4, i) + f_curr_mF[i, :] * dt / 2))[cells_i]
            B_new[i, cells_i] = (solve_banded((1, 1), matrix_A_list[5], b_X(5, i) + f_curr_B[i, :] * dt / 2))[cells_i]
            D_new[i, cells_i] = (solve_banded((1, 1), matrix_A_list[6], b_X(6, i) + f_curr_D[i, :] * dt / 2))[cells_i]
            E_new[i, cells_i] = (solve_banded((1, 1), matrix_A_list[7], b_X(7, i) + f_curr_E[i, :] * dt / 2))[cells_i]
            F_new[i, cells_i] = (solve_banded((1, 1), matrix_A_list[8], b_X(8, i) + f_curr_F[i, :] * dt / 2))[cells_i]
            # print(H_new)
            # print('fasdgh')

        for j in range(J):
            cells_j = np.where(cell_grid[j] != 0)
            mB_new[cells_j, j] = (solve_banded((1, 1), matrix_C_list[0], d_X(0, j) + f_curr_mB[:, j] * dt / 2))[cells_j]
            H_new[cells_j, j] = (solve_banded((1, 1), matrix_C_list[1], d_X(1, j) + f_curr_H[:, j] * dt / 2))[cells_j]
            mD_new[cells_j, j] = (solve_banded((1, 1), matrix_C_list[2], d_X(2, j) + f_curr_mD[:, j] * dt / 2))[cells_j]
            mE_new[cells_j, j] = (solve_banded((1, 1), matrix_C_list[3], d_X(3, j) + f_curr_mE[:, j] * dt / 2))[cells_j]
            mF_new[cells_j, j] = (solve_banded((1, 1), matrix_C_list[4], d_X(4, j) + f_curr_mF[:, j] * dt / 2))[cells_j]
            B_new[cells_j, j] = (solve_banded((1, 1), matrix_C_list[5], d_X(5, j) + f_curr_B[:, j] * dt / 2))[cells_j]
            D_new[cells_j, j] = (solve_banded((1, 1), matrix_C_list[6], d_X(6, j) + f_curr_D[:, j] * dt / 2))[cells_j]
            E_new[cells_j, j] = (solve_banded((1, 1), matrix_C_list[7], d_X(7, j) + f_curr_E[:, j] * dt / 2))[cells_j]
            F_new[cells_j, j] = (solve_banded((1, 1), matrix_C_list[8], d_X(8, j) + f_curr_F[:, j] * dt / 2))[cells_j]


        X_new_list = [mB_new,H_new,mD_new,mE_new,mF_new,B_new,D_new,E_new,F_new]

        for species_index in range(n_species):
            X_new_list[species_index] = np.multiply(X_new_list[species_index],cell_grid)
            mB_new,H_new,mD_new,mE_new,mF_new,B_new,D_new,E_new,F_new = X_new_list
        # print(H_new)

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

    # print(X_record_list[0])

    grids = (x_grid, y_grid, t_grid)
    records = X_record_list
    final_concentration = X_new_list


    return records, final_concentration, grids


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



def plot_redgreen_contrast(final_concentration, mm, mechanism, shape, filename, path, parID=0, dimension='1D', scale_factor=10, save_figure=False):
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
