from matplotlib import pyplot as plt 
from scipy.optimize import leastsq
import numpy as np

# >>>>>> Students: complete code as explained in the text

#  Define model
def A(L,x):
    a,N,eps,Koff,Kon=np.absolute(x[0]), np.absolute(x[1]), -np.absolute(x[2]), np.absolute(x[3]), np.absolute(x[4])
    #a,N,eps,Koff,Kon=x
    ratio=np.divide(1+L/Koff,1+L/Kon)
    logratio=np.log(ratio)
    return np.divide(a, 1+np.exp(N*(eps+np.log(np.divide(1+L/Koff,1+L/Kon))))) 

# Define residuals (errors)
def residuals(x, L, data):
    return data-A(L,x)

p0=1, 2, -2, 0.00001, 0.001
fit=leastsq(residuals, p0, args=(dat1[0], dat1[1]))
print("Best fit parameters: ", fit[0])

# <<<<<<<<<
