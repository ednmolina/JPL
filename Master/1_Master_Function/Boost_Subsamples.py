"""
THIS PROGRAM WILL CALCUALTE THE BOOST OF THE 83 SUBSAMPLES
"""
import numpy as np
from astropy.io import fits as pyfits #Import astro py as this
import matplotlib.pyplot as plt
import pylab
pylab.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size':10})
pylab.rc('text', usetex=True)

#Deletes indicies then calculates boost and write boost to a file
def setBoostFiles(LensJKIndices, LensSums,l_RandomSums, l_RandomJKIndices): #RedMapper, gglens_rand_dr8_rand_zmask_redmapper_v5.10_randcat_z0.05-0.60_lgt020.Z_LAMBDA_0.1_0.33.LAMBDA_lt100.fit.lens/dr8_rand_zmask_redmapper_v5.10_randcat_z0.05-0.60_lgt020.Z_LAMBDA_0.1_0.33.LAMBDA_lt100.%s.fit/sums.fits
    Boost_JackKnife = np.zeros((83, 24))
    for i in range(83):
        #Delete data from sums files
        Index = LensJKIndices
        C = LensSums
        Need_to_Delete = np.where(Index == i)[0]  # Contains the indicies of the data points needed to be deleted
        Realizations = np.zeros((100, 24))
        for j in range(100):
            """
            LOAD RANDOM FILES
            """
            R = l_RandomSums[j]

            """
            LOAD THE INDICIES OF THE RANDOM
            """
            IndicesRandom = l_RandomJKIndices[j]

            Need_to_Delete_From_Random = np.where(IndicesRandom == i)[0] #Sets the indices that need to be deleted

            """
            CALCULATE NUMERATOR FOR BOOST
            """
            Boost_Lensed = C['w']
            New_Boost_Lensed = np.delete(Boost_Lensed, Need_to_Delete, 0)
            numerator = np.sum(New_Boost_Lensed, axis = 0)/len(New_Boost_Lensed)

            """
            CALCULATE DENOMENATOR OR BOOST ***LENGTH OF THE RANDOM FILES IS VARIABLE***
            """
            Boost_Random = R['w']
            New_Boost_Random = np.delete(Boost_Random, Need_to_Delete_From_Random, 0)
            denominator = np.sum(New_Boost_Random, axis = 0)/len(New_Boost_Random)
            """
            CALCULATES THE BOOST
            """
            BoostFactor = numerator/denominator
            Realizations[j] = BoostFactor
        Average_Boost = np.average(Realizations, axis=0)  # Collapses to a 1x24 array
        Boost_JackKnife[i] = Average_Boost
    "Boost Calculating Done"
    return Boost_JackKnife, Average_Boost