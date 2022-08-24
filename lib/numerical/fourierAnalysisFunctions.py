import numpy as np
from scipy.fft import fft, fftfreq
import matplotlib.pyplot as plt
from scipy.stats import entropy


def psEntropyFunction(U0):
    #fourier transform, then power spectrum, remove zero freq, normalize, then entropy
    fft_U0 = fft(U0) #fft of U0
    ps = np.real(fft_U0)**2 + np.imag(fft_U0)**2 #power spectrum of U0
    nonzeroPs = np.round(ps[1:],decimals=3) #remove zero freq
    if np.sum(nonzeroPs)>0: #if there are nonzero freqs normalize
        nonzeroPsNormalized = nonzeroPs/np.sum(nonzeroPs)
    if np.sum(nonzeroPs)==0: #if there are no nonzero freqs, make a uniform distribution
        nonzeroPsNormalized = np.full(len(nonzeroPs),1/len(nonzeroPs))
    psEntropy = entropy(nonzeroPsNormalized)
    return psEntropy

def plotFourier(U0,c='b'):

    fft_U0 = fft(U0)
    ps = np.real(fft_U0)**2 + np.imag(fft_U0)**2
    freq = np.fft.fftfreq(len(U0))
    nonzeroPs = np.round(ps[1:],decimals=3)
    if np.sum(nonzeroPs)>0:
        nonzeroPsNormalized = nonzeroPs/np.sum(nonzeroPs)
    if np.sum(nonzeroPs)==0:
        nonzeroPsNormalized = np.full(len(nonzeroPs),1/len(nonzeroPs))
    plt.plot(freq[1:],nonzeroPsNormalized, c=c)
