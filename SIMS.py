import numpy as np
import scipy.special as sp
def SIMS(x, k, D, t):

    h = k/D
    y = sp.erfc(x/(2*np.sqrt(D*t)))-np.exp(h*x+h**2*D*t)*sp.erfc(x/(2*np.sqrt(D*t))+h*np.sqrt(D*t))

    return y