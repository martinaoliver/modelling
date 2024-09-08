
#%%
import numpy as np
import matplotlib.pyplot as plt
import tensorflow as tf

# Load the .npy file
# loaded_array = np.load('dataset_newsim_2x500x500_crop100.npy')
loaded_array = np.load('dataset_oldsim_2x100x500_crop9087.npy')
train_data = loaded_array[:800,:1].transpose(0,2,3,1)
test_data = loaded_array[900:,:1].transpose(0,2,3,1)
train_data = (train_data - np.mean(train_data, axis=(-3,-2), keepdims=True) )/ (np.std(train_data, axis=(-3,-2), keepdims=True) +1)
test_data = (test_data - np.mean(test_data, axis=(-3,-2), keepdims=True) )/ (np.std(test_data, axis=(-3,-2), keepdims=True) +1)
# train_data = (train_data - np.min(train_data, axis=(-3,-2), keepdims=True)) / (np.max(train_data, axis=(-3,-2), keepdims=True) - np.min(train_data, axis=(-2,-3), keepdims=True))
# Now, `loaded_array` is a numpy array containing the data from the .npy file

print(train_data.shape)
train_data = np.repeat(train_data, 5, axis=1)
test_data = np.repeat(test_data, 5, axis=1)
print(train_data.shape)

#%%

n = 10
plt.figure(figsize=(20, 20))
for i in range(n):
  # display original
  ax = plt.subplot(2, n, i + 1)
  plt.imshow(train_data[i])
  plt.title("original")
  plt.gray()
  ax.get_xaxis().set_visible(False)
  ax.get_yaxis().set_visible(False)

#   # display reconstruction
#   ax = plt.subplot(2, n, i + 1 + n)
#   plt.imshow(decoded_imgs[i])
#   plt.title("reconstructed")
#   plt.gray()
#   ax.get_xaxis().set_visible(False)
#   ax.get_yaxis().set_visible(False)
plt.show()

#%%
noise_factor = 0.2
x_train_noisy = train_data + noise_factor * tf.random.normal(shape=train_data.shape) 
x_test_noisy = test_data + noise_factor * tf.random.normal(shape=test_data.shape) 

x_train_noisy = tf.clip_by_value(x_train_noisy, clip_value_min=0., clip_value_max=1.)
x_test_noisy = tf.clip_by_value(x_test_noisy, clip_value_min=0., clip_value_max=1.)
# %%

BATCH_SIZE = 4 
SHUFFLE_BUFFER_SIZE = 100

train_dataset = train_dataset.shuffle(SHUFFLE_BUFFER_SIZE).batch(BATCH_SIZE)
test_dataset = test_dataset.batch(BATCH_SIZE)

#%%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tensorflow as tf
from sklearn.metrics import accuracy_score, precision_score, recall_score
from sklearn.model_selection import train_test_split
from tensorflow.keras import layers, losses
from tensorflow.keras.datasets import fashion_mnist
from tensorflow.keras.models import Model

class Denoise(Model):
  def __init__(self):
    super(Denoise, self).__init__()
    self.encoder = tf.keras.Sequential([
      layers.Input(shape=(500, 500, 1)),
      layers.Conv2D(32, (3, 3), activation='relu', padding='same', strides=2),
      layers.Conv2D(16, (3, 3), activation='relu', padding='same', strides=2),
      layers.Conv2D(8, (5, 5), activation='relu', padding='same', strides=5),
      layers.Flatten(),
      layers.Dense(3)])

    self.decoder = tf.keras.Sequential([
      layers.Input(shape=(3,)),
      layers.Dense(25*25*8, activation='relu'),
      layers.Reshape((25, 25, 8)),
      layers.Conv2DTranspose(8, kernel_size=5, strides=5, activation='relu', padding='same'),
      layers.Conv2DTranspose(16, kernel_size=3, strides=2, activation='relu', padding='same'),
      layers.Conv2DTranspose(32, kernel_size=3, strides=2, activation='relu', padding='same'),
      layers.Conv2D(1, kernel_size=(3, 3), padding='same')])

  def call(self, x):
    encoded = self.encoder(x)
    decoded = self.decoder(encoded)
    return decoded

autoencoder = Denoise()
#%%
optimizer = tf.keras.optimizers.Adam(
    learning_rate=1e-4)
    
autoencoder.compile(optimizer=optimizer, loss=losses.MeanSquaredError())

#%%

autoencoder.fit(train_data, train_data,
                epochs=20,
                shuffle=True,
                
                )
# %%
