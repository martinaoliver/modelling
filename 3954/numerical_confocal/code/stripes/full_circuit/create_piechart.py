import matplotlib.pyplot as plt
import numpy as np
from decimal import *
import matplotlib as mpl
mpl.rcParams['font.size'] = 15
mpl.use('tkagg')
total = 10000
mylabels = ["No Pattern %.2f %%"%(100-1.56), 'Pattern %.2f %%'%(1.56)]
mycolors = ['lightseagreen','darkslategrey']
y = np.array([100-1.56,1.56])


plt.pie(y, labels = mylabels,colors = mycolors)
plt.tight_layout()
plt.savefig('10000_fullcircuit_piechart_patterns_nonoise.png',bbox_inches='tight')
plt.show()


mpl.rcParams['font.size'] = 13
y = np.array([16,10,130])
percentages = (y/156)*100
mylabels = ['Turing Patterns %.2f %%'%percentages[0], 'Boundary induced patterns %.2f %%'%percentages[1],'Growth induced patterns %.2f %%'%percentages[2]]
mycolors = ['tomato','lightsalmon', 'peachpuff']#,'lightcoral']
plt.pie(y, labels = mylabels,colors = mycolors)
plt.tight_layout()
plt.savefig('10000_fullcircuit_piechart_mechanisms_nonoise.png',bbox_inches='tight')
plt.show()


