import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# scholes_cfMatrix = np.array([[0,100],[100,0]])
# ax = sns.heatmap(scholes_cfMatrix, annot=False, cmap='Greens')
# # ax.set_yticklabels(labels=ax.get_yticklabels(), va='center')
# # ax.set_xticklabels(labels=ax.get_xticklabels(), va='center')
# # # ax.set_xlabel('Dispersion relation')
# # # ax.set_ylabel('Numerical solution')
# ax.set_yticks([.5, 2.5])
# # ## Ticket labels - List must be in alphabetical order
# # ax.xaxis.set_ticklabels(['Classical Turing','Else'])
# # ax.yaxis.set_ticklabels(['Classical Turing','Else'])
# # # plt.tight_layout()
# # ## Display the visualization of the Confusion Matrix.
# # plt.show()



import numpy as np; np.random.seed(0)
import seaborn as sns; sns.set_theme()
uniform_data = np.random.rand(10, 12)
ax = sns.heatmap(uniform_data)
plt.show()