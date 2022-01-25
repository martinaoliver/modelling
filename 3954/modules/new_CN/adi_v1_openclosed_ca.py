import numpy
import matplotlib as mpl
mpl.use('tkagg')

import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import spdiags, diags
import matplotlib.pyplot as plt
from tqdm import tqdm
import copy
from scipy.sparse import linalg
from scipy.linalg import solve_banded

n_species=3

#spatial variables
L_x = 8; L_y=L_x
# L_x = 1; L_y=L_x
T =  50
J = L_x*4; I=J
# J = L_x*5; I=J
dx = float(L_x)/float(J-1)
dy = float(L_y)/float(I-1)
x_grid = numpy.array([j*dx for j in range(J)])
y_grid = numpy.array([i*dy for i in range(I)])
D = np.array([ 0.125,0.5,0])
diffusing_species =np.nonzero(D)[0]
nondiffusing_species = np.nonzero(D==0)[0]
#time variables
# N = T*10 #try to have uneven number of gridpoints so that initial plug stays in the middle. Otherwise, initial cell defined as 2 cells
def t_gridpoints_stability(L,dx,T):
    N = T/(0.49*(dx**2)) + 1
    return int(N)
N = t_gridpoints_stability(L_x,dx,T)
N=T*5
dt = float(T)/float(N-1)
t_grid = numpy.array([n*dt for n in range(N)])



#Define initial conditions
U0 = []
perturbation=0.001
steadystates=[0.1,0.1,0.1]
np.random.seed(2)

cell_matrix = np.zeros(shape=(I,J))
cell_matrix[int(I/2), int(J/2)] = 1
cell_matrix
for index in range(n_species):
    U0.append(np.random.uniform(low=steadystates[index] - perturbation, high=steadystates[index] + perturbation, size=(I, J)))
    print(U0)
U0 = U0*cell_matrix
U0
alpha = [D[n]*dt/(2.*dx*dx) for n in diffusing_species]

#A matrix (right-hand side of Ax=b)
def A(alphan):
    bottomdiag = [-alphan for j in range(J-1)]
    centraldiag = [1.+alphan]+[1.+2.*alphan for j in range(J-2)]+[1.+alphan]
    topdiag = [-alphan for j in range(J-1)]
    diagonals = [bottomdiag,centraldiag,topdiag]
    A = diags(diagonals, [ -1, 0,1]).toarray()
#     A = scipy.sparse.csc_matrix(A, dtype=float)
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

#define ODE equations
def f(U, cell_matrix):
    dudt = [0]*n_species
    # dudt[0]= 5*U[0] - 6*U[1] + 1
    # dudt[1] = 6*U[0]- 7*U[1]
    dudt[0]= 5*1
    dudt[1] = 6*1
    dudt[2] =  U[0]*0
    dudt = [eq*cell_matrix for eq in dudt]
    return dudt

def check_neighbours(cell_matrix,y_pos,x_pos): #returns grid with the neighbouring points
    top_array = [cell_matrix[y_pos-1, x_pos-1], cell_matrix[y_pos-1,x_pos], cell_matrix[y_pos-1,x_pos+1]]
    middle_array = [cell_matrix[y_pos, x_pos-1], np.nan, cell_matrix[y_pos,x_pos+1]]
    # print('afe')
    # print(cell_matrix[2,1])
    bottom_array = [cell_matrix[y_pos+1, x_pos-1], cell_matrix[y_pos+1,x_pos], cell_matrix[y_pos+1,x_pos+1]]
    neighbours_cellmatrix = np.array([top_array,middle_array,bottom_array])
    return neighbours_cellmatrix

def cell_automata_colony(species_list,cell_matrix, p_division):
    new_species_list = copy.deepcopy(species_list)
    original_species_list = copy.deepcopy(species_list)
    cell_matrix_new = copy.deepcopy(cell_matrix)
    for y_pos in np.linspace(1,len(cell_matrix)-2,len(cell_matrix)-2):
    # for x_pos in np.linspace(0,len(cell_matrix)-2,len(cell_matrix)-2):
        for x_pos in np.linspace(1,len(cell_matrix)-2,len(cell_matrix)-2):
        # for y_pos in np.linspace(0,len(cell_matrix)-2,len(cell_matrix)-2):
            y_pos = int(y_pos)
            x_pos = int(x_pos)
            # print(y_pos,x_pos,cell_matrix)
            if cell_matrix[y_pos, x_pos]!=0:
                print(y_pos,x_pos)
                neighbours_cellmatrix = check_neighbours(cell_matrix,y_pos,x_pos)
                print(neighbours_cellmatrix)
                if 0 in neighbours_cellmatrix:
                    cell_division=np.random.choice([1,0],p=[p_division,1-p_division])
                    if cell_division==1:
                        index_nocells=np.where(np.array(neighbours_cellmatrix )== 0)
                        print(index_nocells)
                        divided_cell_index = np.random.choice(range(len(index_nocells[0])))
                        index_newcell_y, index_newcell_x = (index_nocells[n][divided_cell_index] for n in range(2))
                        # index_newcell_y , index_newcell_x = 1,2

                        print(index_newcell_y, index_newcell_x)# count=0
                        # for species,new_species in zip(original_species_list,new_species_list):
                        print(index_newcell_y+y_pos-1,index_newcell_x+x_pos-1)
                        for count,species in enumerate(original_species_list):
                            new_species_list[count][index_newcell_y+y_pos-1,index_newcell_x+x_pos-1] += species[y_pos,x_pos]/2
                            new_species_list[count][y_pos,x_pos] += species[y_pos,x_pos]/2
                            # count+=1
                        cell_matrix_new[index_newcell_y+y_pos-1,index_newcell_x+x_pos-1]=1
                        print(cell_matrix_new)


    return new_species_list, cell_matrix_new





U = copy.deepcopy(U0)

ab_list = [diagonal_form(A(alphan)) for alphan in alpha]# = [A(alphan) for alphan in alpha]

for ti in tqdm(range(N)):

    # if ti%10==0:
        # plt.subplot(121)
        # plt.imshow(cell_matrix)#, vmin=0.09, vmax=0.2)
        # plt.subplot(122)
        # plt.imshow(U[0])#, vmin=0.09, vmax=0.2)
        # plt.colorbar()
        # plt.tight_layout()
        # plt.show()

    U_half = copy.deepcopy(U)
    f0 = f(U,cell_matrix)
    for i in range(I):
        for n in diffusing_species:
            # U_half[n][:,i] = solve_banded((1, 1), ab_list[n], b('y',i,alpha[n],U[n])) +  f(U)[n][:,i]*(dt/2)
            # U_half[n][:,i] = solve_banded((1, 1), ab_list[n], b('y',i,alpha[n],U[n]) +  f(U,cell_matrix)[n][:,i]*(dt/2))
            U_half[n][:,i] = solve_banded((1, 1), ab_list[n], b('y',i,alpha[n],U[n]) +  f0[n][:,i]*(dt/2))
        for n in nondiffusing_species:
            # print(n)
            U_half[n][:,i] =  U[n][:,i] + f0[n][:,i]*(dt/2)#f(U)[n][:,i]*(dt/2)
            # U_half[n][:,i] =  U[n][:,i] + 1#f(U)[n][:,i]*(dt/2)

    U_new = copy.deepcopy(U_half)
    f1 = f(U_half,cell_matrix)
    for j in range(J):
        for n in diffusing_species:
            # U_new[n][j,:] = solve_banded((1, 1), ab_list[n], b('x',j,alpha[n],U_half[n])) + f(U_half)[n][j,:]*(dt/2)
            # U_new[n][j,:] = solve_banded((1, 1), ab_list[n], b('x',j,alpha[n],U_half[n]) + f(U_half,cell_matrix)[n][j,:]*(dt/2))
            U_new[n][j,:] = solve_banded((1, 1), ab_list[n], b('x',j,alpha[n],U_half[n]) + f1[n][j,:]*(dt/2))
        for n in nondiffusing_species:
            # print(n)
            U_new[n][j,:] =  U_half[n][j,:] + f1[n][j,:]*(dt/2)#f(U_half)[n][j,:]*(dt/2)
            # U_new[n][j,:] =  U_half[n][j,:] + 1#f(U_half)[n][j,:]*(dt/2)

    hour = ti / (N / T)
    if hour % 1 == 0:  #only consider division at unit time (hour)
        plt.subplot(121)
        plt.imshow(cell_matrix)#, vmin=0.09, vmax=0.2)
        plt.subplot(122)
        plt.imshow(U[0])#, vmin=0.09, vmax=0.2)
        plt.colorbar()
        plt.tight_layout()
        plt.show()
        p_division=1
        U_new, cell_matrix_new = cell_automata_colony(U_new, cell_matrix, p_division)
        plt.imshow(cell_matrix)
        cell_matrix = copy.deepcopy(cell_matrix_new)

    # print(U[2])
    U = copy.deepcopy(U_new)


#
#
# plt.imshow(U[0])
# plt.colorbar()
# plt.show()
#
#
#
# plt.imshow(U[1])
# plt.colorbar()
# plt.show()
#
#
#
#
# plt.imshow(U[2])
# plt.colorbar()
# plt.show()
