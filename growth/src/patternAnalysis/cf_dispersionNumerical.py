import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt


cf_matrix_scholes = np.array([[1,0],[0,1]])
ax = sns.heatmap(cf_matrix_scholes, annot=True, cmap='Blues',cbar=False)

ax.set_xlabel('\nNumerical solution')
ax.set_ylabel('Dispersion relation\n');

## Ticket labels - List must be in alphabetical order
ax.yaxis.set_ticklabels(['Turing I','Else'])
ax.xaxis.set_ticklabels(['Pattern','No Pattern'])

## Display the visualization of the Confusion Matrix.
plt.tight_layout()
plt.savefig('cf_scholes')
plt.show()




cf_matrix_marti = np.array([[1,0,0,0,0],[0,0,0,0,1],[1,0,1,0,0],[0,1,1,0,0],[0,0,0,1,0],[0,0,0,0,1]])
ax = sns.heatmap(cf_matrix_marti, annot=True, cmap='Blues',cbar=False)

ax.set_xlabel('\nNumerical solution')
ax.set_ylabel('Dispersion relation\n');

## Ticket labels - List must be in alphabetical order
ax.yaxis.set_ticklabels(['Turing I','Turing II', 'Turing I Hopf', 'Hopf instability', 'Complex with no instability', 'Real with no instability'], rotation=0)
ax.xaxis.set_ticklabels(['Stationary spatial wave','Oscillatory spatial wave', 'Temportal oscillations', 'Chaotic inhomogeneities', 'Homogeneous'], rotation=90)
## Display the visualization of the Confusion Matrix.
plt.tight_layout()
plt.savefig('cf_marti')
plt.show()



