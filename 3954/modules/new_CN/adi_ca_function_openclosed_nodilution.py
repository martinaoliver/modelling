import numpy


import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import spdiags, diags
import matplotlib.pyplot as plt
from tqdm import tqdm
import copy
from scipy.sparse import linalg
from scipy.linalg import solve_banded
from class_circuit_eq import *

def adi_ca_openclosed_nodilution(par_dict,L_x,L_y,J,I,T,N, circuit_n, n_species,D,tqdm_disable=False, p_division=0.5,stochasticity=0, seed=1,growth='Slow', boundarycoeff=1.5):

    parent_list = [circuit1_eq, circuit2_eq,circuit3_eq,circuit4_eq,circuit5_eq,circuit6_eq,circuit7_eq,circuit8_eq,circuit9_eq, circuit10_eq, circuit11_eq]
    f = parent_list[circuit_n-1](par_dict, stochasticity=stochasticity)

    #spatial variables
    dx = float(L_x)/float(J-1); dy = float(L_y)/float(I-1)
    x_grid = numpy.array([j*dx for j in range(J)]); y_grid = numpy.array([i*dy for i in range(I)])
    diffusing_species =np.nonzero(D)[0]
    nondiffusing_species = np.nonzero(D==0)[0]

    # #time variables
    # def t_gridpoints_stability(T, dx):
    #     N = T/(0.49*(dx**2)) + 1
    #     return int(N)
    # N = t_gridpoints_stability(L_x,dx,T)

    dt = float(T)/float(N-1)
    t_grid = numpy.array([n*dt for n in range(N)])

    alpha = [D[n]*dt/(2.*dx*dx) for n in range(n_species)]


    #Define initial conditions and cell matrix
    U0 = []
    perturbation=0.001
    steadystates=[0.1]*n_species
    np.random.seed(seed)

    cell_matrix = np.zeros(shape=(I,J))
    cell_matrix[int(I/2), int(J/2)] = 1
    for index in range(n_species):
        U0.append(np.random.uniform(low=steadystates[index] - perturbation, high=steadystates[index] + perturbation, size=(I, J)))
    U0 = U0*cell_matrix


    #A matrix (right-hand side of Ax=b)
    def A(alphan):
        bottomdiag = [-alphan for j in range(J-1)]
        centraldiag = [1.+boundarycoeff*alphan]+[1.+2.*alphan for j in range(J-2)]+[1.+boundarycoeff*alphan]
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


    #b vector (left-hand side of Ax=b)
    def b(axis,ij,alphan,Un):
        b_t_stencil = np.array( [0] + [(1-boundarycoeff*alphan)] + [alphan])
        b_c_stencil = np.array( [alphan] + [(1-2*alphan)] + [alphan])
        b_b_stencil = np.array( [alphan] + [(1-boundarycoeff*alphan)] + [0])

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


    def check_neighbours(cell_matrix,y_pos,x_pos): #returns grid with the neighbouring points
        top_array = [cell_matrix[y_pos-1, x_pos-1], cell_matrix[y_pos-1,x_pos], cell_matrix[y_pos-1,x_pos+1]]
        middle_array = [cell_matrix[y_pos, x_pos-1], np.nan, cell_matrix[y_pos,x_pos+1]]
        bottom_array = [cell_matrix[y_pos+1, x_pos-1], cell_matrix[y_pos+1,x_pos], cell_matrix[y_pos+1,x_pos+1]]
        neighbours_cellmatrix = np.array([top_array,middle_array,bottom_array])
        return neighbours_cellmatrix
    def cell_automata_colony(species_list,cell_matrix, p_division):
        new_species_list = copy.deepcopy(species_list)
        original_species_list = copy.deepcopy(species_list)
        cell_matrix_new = copy.deepcopy(cell_matrix)
        for y_pos in np.linspace(1,len(cell_matrix)-2,len(cell_matrix)-2):
            for x_pos in np.linspace(1,len(cell_matrix)-2,len(cell_matrix)-2):
                y_pos = int(y_pos)
                x_pos = int(x_pos)
                if cell_matrix[y_pos, x_pos]!=0:
                    neighbours_cellmatrix = check_neighbours(cell_matrix,y_pos,x_pos)
                    if 0 in neighbours_cellmatrix:
                        cell_division=np.random.choice([1,0],p=[p_division,1-p_division])
                        if cell_division==1:
                            index_nocells=np.where(np.array(neighbours_cellmatrix )== 0)
                            divided_cell_index = np.random.choice(range(len(index_nocells[0])))
                            index_newcell_y, index_newcell_x = (index_nocells[n][divided_cell_index] for n in range(2))
                            for count,species in enumerate(original_species_list):
                                # new_species_list[count][index_newcell_y+y_pos-1,index_newcell_x+x_pos-1] += species[y_pos,x_pos]/2
                                # new_species_list[count][y_pos,x_pos] += species[y_pos,x_pos]/2
                                new_species_list[count][index_newcell_y+y_pos-1,index_newcell_x+x_pos-1] += species[y_pos,x_pos]/2
                                new_species_list[count][y_pos,x_pos] += species[y_pos,x_pos]/2
                            cell_matrix_new[index_newcell_y+y_pos-1,index_newcell_x+x_pos-1]=1


        return new_species_list, cell_matrix_new


    U = copy.deepcopy(U0)
    U_record = []
    for species_index in range(n_species):
        U_record.append(np.zeros([J, I, T])) #DO NOT SIMPLIFY TO U_record = [np.zeros([J, I, T])]*n_species

    t_gridpoints=(N/T)
    tdivider=0.1*t_gridpoints #every 0.1 hour we divide. 
    for ti in tqdm(range(N), disable = tqdm_disable):
        #First step: solve in y direction from n -> n+1/2
        U_half = copy.deepcopy(U)
        f0 = f.dudt_growth(U,cell_matrix)
        for i in range(I):
            for n in diffusing_species:
                U_half[n][:,i] = solve_banded((1, 1), ab_list[n], b('y',i,alpha[n],U[n]) +  f0[n][:,i]*(dt/2)) #CN step in one dimension to get banded(tridiagonal) A matrix
            for n in nondiffusing_species: #Species with D=0, no CN included, basic euler
                U_half[n][:,i] =  U[n][:,i] + f0[n][:,i]*(dt/2)

        #Second step: solve in x direction from n+1/2 -> n+1
        U_new = copy.deepcopy(U_half)
        f1 = f.dudt_growth(U_half,cell_matrix)
        for j in range(J):
            for n in diffusing_species:
                U_new[n][j,:] = solve_banded((1, 1), ab_list[n], b('x',j,alpha[n],U_half[n]) + f1[n][j,:]*(dt/2))
            for n in nondiffusing_species:
                U_new[n][j,:] =  U_half[n][j,:] + f1[n][j,:]*(dt/2)

        hour = ti / (N / T)

        if hour % 1 == 0:  #only consider recording at unit time (hour)
            #append results into top_array for records
            for species_index in range(n_species):
                U_record[species_index][:, :, int(hour)] = U_new[species_index] #issue in this line
            if growth=='Slow': 
                #predict if division occurs based on the p_division, the current cell matrix
                #return new cell matrix and updated concentrations with dilution
                U_new, cell_matrix_new = cell_automata_colony(U_new, cell_matrix, p_division)
                cell_matrix = copy.deepcopy(cell_matrix_new)
        if growth=='Fast':
            if (ti%tdivider==0):
                #predict if division occurs based on the p_division, the current cell matrix
                #return new cell matrix and updated concentrations with dilution
                U_new, cell_matrix_new = cell_automata_colony(U_new, cell_matrix, p_division)
                cell_matrix = copy.deepcopy(cell_matrix_new)

        U = copy.deepcopy(U_new)

    return U_record, U
