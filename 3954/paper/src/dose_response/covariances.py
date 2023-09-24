import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

# Set the seed for reproducibility
np.random.seed(0)

# Generate random data for x and y
mean = [0, 0]
covariance_negative = [[1, -0.8], [-0.8, 1]]
covariance_zero = [[1, 0], [0, 1]]
covariance_positive = [[1, 0.8], [0.8, 1]]

data_negative = np.random.multivariate_normal(mean, covariance_negative, 10000)
data_zero = np.random.multivariate_normal(mean, covariance_zero, 10000)
data_positive = np.random.multivariate_normal(mean, covariance_positive, 10000)

# Extract x and y from the generated data
x_negative, y_negative = data_negative[:, 0], data_negative[:, 1]
x_zero, y_zero = data_zero[:, 0], data_zero[:, 1]
x_positive, y_positive = data_positive[:, 0], data_positive[:, 1]

# Calculate distances from the mean for each point
distance_negative = np.sqrt((x_negative - mean[0])**2 + (y_negative - mean[1])**2)
distance_zero = np.sqrt((x_zero - mean[0])**2 + (y_zero - mean[1])**2)
distance_positive = np.sqrt((x_positive - mean[0])**2 + (y_positive - mean[1])**2)
#%%
# Create subplots for the three distributions
plt.figure(figsize=(9, 2.5))

plt.subplot(131)
plt.scatter(x_negative, y_negative, c=distance_negative, cmap='viridis', alpha=0.5)
plt.title('Covariance X,Y < 0')
plt.xlabel('X')
plt.ylabel('Y')
plt.colorbar(label='Distance from Mean')

plt.subplot(132)
plt.scatter(x_zero, y_zero, c=distance_zero, cmap='viridis', alpha=0.5)
plt.title('Covariance X,Y = 0')
plt.xlabel('X')
plt.ylabel('Y')
plt.colorbar(label='Distance from Mean')

plt.subplot(133)
plt.scatter(x_positive, y_positive, c=distance_positive, cmap='viridis', alpha=0.5)
plt.title('Covariance X,Y > 0')
plt.xlabel('X')
plt.ylabel('Y')
plt.colorbar(label='Distance from Mean')

plt.tight_layout()
plt.savefig('covariances.pdf')
plt.show()

# %%
