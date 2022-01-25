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
L_x = 1; L_y=L_x
T =  10
J = L_x*4; I=J
J = L_x*5; I=J
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
np.random.seed(1)

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
    dudt[0]= 5*U[0] - 6*U[1] + 1
    dudt[1] = 6*U[0]- 7*U[1]
    dudt[2] =  U[0]*0
    dudt = [eq*cell_matrix for eq in dudt]
    return dudt




U = copy.deepcopy(U0)

ab_list = [diagonal_form(A(alphan)) for alphan in alpha]# = [A(alphan) for alphan in alpha]

for ti in tqdm(range(N)):

    # if ti%10==0:
    #     plt.imshow(U[0])#, vmin=0.09, vmax=0.2)
    #     plt.colorbar()
    #     plt.show()

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

    print(U[2])
    U = copy.deepcopy(U_new)


plt.imshow(U[0])
plt.colorbar()
plt.show()



plt.imshow(U[1])
plt.colorbar()
plt.show()




plt.imshow(U[2])
plt.colorbar()
plt.show()
