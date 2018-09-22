import numpy as np
import pylab
import sys
sys.path.append("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/1_Master_Function")
import CovarianceAndCorrelation
from astropy.io import fits as pyfits #Import astro py as this
import cosmology
import matplotlib.pyplot as plt
pylab.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size':10})
pylab.rc('text', usetex=True)

"""
IMPORT FILES
"""
#Uncorrected Delta Sigma
UCDSEM = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutTest/DeltaSigmaFiles/Uncorrected_Calculated_DeltaSigma/UncorrectedEMJKSS.txt")
UCDSBM = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutTest/DeltaSigmaFiles/Uncorrected_Calculated_DeltaSigma/UncorrectedBMJKSS.txt")

#Load the Boost
Boost = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutTest/Boost/Boost.txt")

#Load the Photo-Z
Photo_Z = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutTest/PhotoZ/PhotoZ.txt")

#Load the Randoms
RandomEM = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutTest/Random/RandomEM.txt")
RandomBM = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutTest/Random/RandomBM.txt")

"""
*** PERFORM THE CORRECTIONS ***
"""

"""
Boost Correction
"""
Boost_DSEM = UCDSEM * Boost
#Calculate Covariance
Boost_Cov_EM = CovarianceAndCorrelation.getCovariance(Boost_DSEM)

#Save the correction
np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages/EM_Boost.txt", np.c_[Boost_DSEM])
np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages/EM_BoostCov.txt", np.c_[Boost_Cov_EM])

"""
Photo-Z Correction
"""
Photo_Z_DSEM = Boost_DSEM * Photo_Z
Photo_Z_DSEMCov = CovarianceAndCorrelation.getCovariance(Photo_Z_DSEM)
#Save the correction
np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages/EM_Boost_PhotoZ.txt", np.c_[Photo_Z_DSEM])
np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages/EM_Boost_PhotoZCov.txt", np.c_[Photo_Z_DSEMCov])