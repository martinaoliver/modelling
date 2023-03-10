

import matplotlib.pyplot as plt
def noncompetitiveact(X,args):
    b, V, km, n = args
    V = b + Vmax*((X / km) ** (n)) / (1 + (X / km) ** (n))
    # V = b + Vmax*(1) / (1 + (km / X) ** (n))
    return V

Vmax=100
b=0.1
km=0.4
n=1
args = (b, Vmax, km, n)

X = np.linspace(0, 20, 500)
plt.axvline(x=km, color='darkcyan',ls=':', lw=2, label='axvline - full height')
plt.axvline(x=6.03, color='black',ls=':', lw=2, label='axvline - full height')

plt.plot(X, noncompetitiveact(X, args), c='darkcyan')
plt.xlabel('[HSL]')
plt.ylabel('GFP')
plt.show()
