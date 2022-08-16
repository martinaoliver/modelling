import numpy as np
from scipy.fft import fft, ifft

def powerspectrumFunction(fft_U):
    #(re^2 + im^2/N  This powerspectrum is normalised to 1
    sumsq = np.real(fft_U)**2 + np.imag(fft_U)**2
    N = np.sum(sumsq)
    powerspectrum = sumsq/N
    return powerspectrum
def entropyFunction(vector):
    H = - np.sum ( [p * np.log(p) for p in vector if p != 0]  )
    return H

def fourierAnalysisFunction(U0):
    #fourier transform, then power spectrum, then entropy
    fft_U = fft(U0)
    ps = powerspectrumFunction(fft_U)
    H = entropyFunction(ps)
    return H
