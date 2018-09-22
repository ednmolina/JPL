"""
CALCULATES THE RANDOM PSF OF THE 83 SUBSAMPLES
"""
import numpy as np

import pylab
pylab.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size':10})
pylab.rc('text', usetex=True)

"""
CALCULATING RANDOM OF E MODE
"""
def setRandomFiles(Ran_Sums, l_RandomJKIndices, WetWex):
    Random_DeltaSigma = np.zeros((83, 24))

    for i in range(83):
        Delta_Sigma = np.zeros((100, 24))
        for j in range(100):
            #Indicies of Random
            index = l_RandomJKIndices[j]
            Need_to_Delete_Indices = np.where(index == i)[0]
            #RANDOM 100 sums.fits
            S = Ran_Sums[j]

            """DELETES CORRESPONDING INSICIES AND SUMS UP DATA"""
            A = np.delete(S["%s" %WetWex], Need_to_Delete_Indices, 0)
            B = np.delete(S["w"], Need_to_Delete_Indices, 0)
            A_Sum = np.sum(A, axis = 0)
            B_Sum = np.sum(B, axis=0)

            Delta_Sigma[j]=((0.5) * (1. / (1. - (0.365) ** 2.)) *A_Sum/B_Sum)

        Random_DeltaSigma[i] = np.average(Delta_Sigma, axis=0)

    print "Randoms Calculated"
    return Random_DeltaSigma
