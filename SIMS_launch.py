import numpy as np
from MCMCSIMSig import MCMCSIMSig
#MCMCSIMS tutorial
#A dataset is included for this tutorial. It is pH2O-dependent 1.txt. This SIMS data set was produced by Dr. Roger
#De-Souze that has already been analyzed
#Tutorial is based on pH2O = 220mbar experiment

#Importing Data
H220x = np.loadtxt('H220x.txt', dtype=float, delimiter=',')
H220z = np.loadtxt('H220z.txt', dtype=float, delimiter=',')

#First column is set of sample thicknesses in cm
#Second column is the set of O18 fractions

H220t = 5640

kmin = 1e-15
kmax = 1e-2
k = 1e-6
SIGMAk = .1

Dmin = 1e-15
Dmax = 1e-2
D = 5e-6
SIGMAD = .1
ps = 0.03
N = 5000
nu = 1000
tau = 1001*ps**2
thinfact = 1

MCMCSIMSig(H220x, H220z, H220t, kmin, kmax, k, SIGMAk, Dmin, Dmax, D, SIGMAD, N, nu, tau, thinfact)