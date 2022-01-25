import matplotlib.pyplot as plt
import numpy as np

y = np.array([99.61, 0.39])
mylabels = ["No Turing pattern 99.61%", "Turing pattern 0.39%"]
colors = ['lightseagreen', 'darkslategrey']
plt.pie(y, labels = mylabels, colors=colors)
plt.tight_layout()
plt.show()
