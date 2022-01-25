import numpy as np
from scipy.sparse import spdiags, eye
from scipy import signal
import matplotlib.pyplot as plt
import scipy
import timeit
from scipy.linalg import solve, lu_factor, lu_solve
import matplotlib as mpl
mpl.use('tkagg')


D=1

#space
D=1; Lx=1;Ly=1 #Diffusion constant; Domain length
Jx = 10; Jy=Jx #Number of cells in x and y
dx = float(Lx/(Jx-1)); dy = float(Ly/(Jy-1)) #space between each cell
ncells = Jx*Jy #number of cells in simulation
#time
T = 100 #total time in units
N = T*10#total number of timesteps
dt = float(T/(N-1)) #time of each timesteps in units


# dt=(dx*dy)/(2*D); #borderline stability of FTCS scheme
print(dt)
alpha=dt*D/(dx*dy);


nx=Jx


def no(coordinates,n=Jx):
    i=coordinates[0]
    j = coordinates[1]
    k = i + (j)*n
    return k

top  = [(0,j) for j in range(Jx)]
bottom  = [(Jx-1,j) for j in range(Jx)]
left =  [(i,0) for i in range(Jy)]
right  = [(i,Jy-1) for i in range(Jy)]
top_no = [no(top_cell) for top_cell in top]
bottom_no = [no(bottom_cell) for bottom_cell in bottom]
left_no = [no(left_cell) for left_cell in left]
right_no = [no(right_cell) for right_cell in right]


interior_no = np.linspace(0,ncells-1,ncells)
boundaries_no = top_no + bottom_no + left_no + right_no
interior_no = [int(x) for x in interior_no if x not in boundaries_no]
# interior_no



ncells = Jx*Jy
# ncells = Jx

data = np.array([-alpha*np.ones(ncells),-alpha*np.ones(ncells), 2*(1+2*alpha)*np.ones(ncells), -alpha*np.ones(ncells),-alpha*np.ones(ncells),])
diags = np.array([-Jx, -1, 0, 1, Jx])
I = eye(ncells).toarray()
A = spdiags(data, diags,ncells,ncells).toarray()
A[boundaries_no] = I[boundaries_no]



# sigma=Lx/1; #stedv of the gaussian - sharp gaussian as initial condition
# u = signal.gaussian(ncells,sigma)*100
# u=u.transpose()
# plt.plot(np.linspace(0,ncells,ncells),u)# xlabel('$x$','Interpreter','latex','FontSize',14);
# plt.title('initial condition')
u = np.zeros(ncells)
center_no = no((int(Jx/2),int(Jx/2)))
u[center_no] = 100
u0 = u

def f(u):
    return 1

from scipy.linalg import lu_factor, lu_solve
from tqdm import tqdm
luA, piv = lu_factor(A)
start = timeit.timeit()
b=[0]*ncells
for m in tqdm(range(N)):
    # plt.plot(np.linspace(0,ncells,ncells),u,label=m)
    for n in interior_no:
        b[n] = 2*(1-2*alpha)*u[n] + alpha*(u[n+1] + u[n-1] + u[n+Jx] + u[n-Jx])
    u = scipy.linalg.lu_solve((luA,piv),b)

    # u = scipy.linalg.lu_solve((luA,piv),b)
end = timeit.timeit()
# plt.legend()
np.linspace(0,ncells,ncells)
# plt.xlim(0,10)
print(end - start)


def plot_matrix(vector,Jx=Jx):
    u_matrix = vector.reshape(Jx,Jx)
    plt.figure()
    plt.title('initial condition')
    plt.imshow(u_matrix)
    plt.colorbar()
    plt.show()

plot_matrix(u,Jx)
plot_matrix(u0,Jx)
