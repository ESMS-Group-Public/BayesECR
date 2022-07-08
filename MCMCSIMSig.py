import numpy as np
from LikelihoodECRrect import LikelihoodECRrect
from SIMS import SIMS
from ig import ig
from Metropolis import METROPOLIS

def MCMCSIMSig(x,z,t,kmin,kmax,k,SIGMAk,Dmin,Dmax,D,SIGMAD,N,nu,tau,thinfact):
    """
    
    :param x: depth, in cm
    :param z: tracer site fraction, at distances corresponding to x
    :param t: time associated with the SIMS dataset
    :param kmin: The smallest permissible k*
    :param kmax: The largest permissible k*
    :param k: Initial value of k*
    :param SIGMAk: Standard deviation of the proposal draws on k
    :param Dmin: The smallest permissible D*
    :param Dmax: The largest permissible D*
    :param D: Initial value of D*
    :param SIGMAD: Standard Deviation of the proposal draws on D
    :param N: Number of cycles to run the calibration (suggested = 100000. There is a burn-in of around 5000, so
        the program returns the results from the lastN-5000)
    :param nu: shape parameter of the inverse gamma prior for the observation error variance (psi)
    :param tau: scale parameter of the inverse gamma prior for the observation error variance
    :param thinfact: thinning factor, fraction of the data you want to use in the calibration, number between 0 and 1
    :return: kpost : the posterior k* distribution
            Dpost: the posterior D* distribution
            psipost: the posterior psi distribution
            SUCCESSk = The fraction of accepted k*
            SUCCESSD = The fraction of accepted D*
    """
    print(''' MCMCECRig is free software: you can redistribute it and/or modify 
    it under the terms of the GNU General Public License as published by the Free Software Foundation, either 
    version 3 of the License, or (at your option) any later version\n.''')

    #Preallocation
    kpost = np.zeros((1,N))
    Dpost = np.zeros((1,N))
    psipost = np.zeros((1,N))
    likelies = np.zeros((1,N))

    kpost[0][0] = k
    Dpost[0][0] = D
    psipost[0][0] = tau/(nu+1)

    acceptedk = 0
    acceptedD = 0

    totacck = 0
    totaccD = 0

    y = SIMS(x,k,D,t)

    LOGLIKELY, EPSILON = LikelihoodECRrect(y,z,t,psipost[0][0])

    likelies[0][0] = LOGLIKELY

    # Monte Carlo one-at-a-time iterations
    for i in range(1,N):
        if np.mod(i,100) == 0:
            print('currently on draw {:d} of {:d} total draws\n'.format(i,N))
            print('acceptance rate k =\n')
            print(acceptedk/100)
            print('acceptance rate D =\n')
            print(acceptedD/100)

            if acceptedk/100 < 0.01:
                SIGMAk = SIGMAk/2
            elif acceptedk/100 > 0.1:
                SIGMAk = 1.5*SIGMAk
            if acceptedD/100 < 0.01:
                SIGMAD = SIGMAD/2
            elif acceptedD/100 > 0.1:
                SIGMAD = 1.5*SIGMAD
            totacck = totacck + acceptedk
            totaccD = totaccD + acceptedD
            acceptedD = 0
            acceptedk = 0

        proposedk = np.log10(kmin) - 1
        while proposedk < np.log10(kmin) or proposedk > np.log10(kmax):
            proposedk = np.random.normal(np.log10(kpost[0][i-1]),SIGMAk)

        proposedk = 10**(proposedk)

        proposedy = SIMS(x, proposedk, Dpost[0][i-1],t)

        proposedLikely, propsilon = LikelihoodECRrect(proposedy,z,t,psipost[0][i-1])

        if METROPOLIS(proposedLikely,LOGLIKELY, thinfact):
            kpost[0][i] = proposedk
            LOGLIKELY = proposedLikely
            EPSILON = propsilon
            acceptedk = acceptedk + 1
        else:
            kpost[0][i] = kpost[0][i-1]

        proposedD = np.log10(Dmin) - 1

        while proposedD < np.log10(Dmin) or proposedD > np.log10(Dmax):
            proposedD = np.random.normal(np.log10(Dpost[0][i-1]),SIGMAD)

        proposedD = 10**proposedD

        proposedy = SIMS(x,kpost[0][i-1],proposedD,t)

        proposedLikely,propsilon = LikelihoodECRrect(proposedy,z,t,psipost[0][i-1])

        if METROPOLIS(proposedLikely, LOGLIKELY, thinfact):
            Dpost[0][i] = proposedD
            LOGLIKELY = proposedLikely
            EPSILON = propsilon
            acceptedD = acceptedD + 1
        else:
            Dpost[0][i] = Dpost[0][i - 1]

        psipost[0][i] = ig(nu,tau,EPSILON)
        likelies[0][i] = LOGLIKELY

    successk = totacck/N
    successD = totaccD/N
    weights = (1-thinfact)*likelies
    weights = weights - np.max(weights)
    weights = np.exp(weights)

    return kpost, Dpost, psipost, weights, successk, successD