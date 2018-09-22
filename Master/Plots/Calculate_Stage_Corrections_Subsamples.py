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
LOAD THE UNCORRECTED SUBSAMPLE DELTA SIGMAS
"""
#Large Red Subsample
L_R_Uncorrected = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleRLarge/DeltaSigmaFiles/Uncorrected_Calculated_DeltaSigma/UncorrectedEMJKSS.txt")
L_R_Boost = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleRLarge/Boost/Boost.txt")
L_R_PhotoZ = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleRLarge/PhotoZ/PhotoZ.txt")
L_R_Random = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleRLarge/Random/RandomEM.txt")

#Small Red Subsample
S_R_Uncorrected = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleRSmall/DeltaSigmaFiles/Uncorrected_Calculated_DeltaSigma/UncorrectedEMJKSS.txt")
S_R_Boost = np.loadtxt("//Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleRSmall/Boost/Boost.txt")
S_R_PhotoZ = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleRSmall/PhotoZ/PhotoZ.txt")
S_R_Random = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleRSmall/Random/RandomEM.txt")

#Large Green Subsample
L_G_Uncorrected = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleGLarge/DeltaSigmaFiles/Uncorrected_Calculated_DeltaSigma/UncorrectedEMJKSS.txt")
L_G_Boost = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleGLarge/Boost/Boost.txt")
L_G_PhotoZ = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleGLarge/PhotoZ/PhotoZ.txt")
L_G_Random = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleGLarge/Random/RandomEM.txt")

#Small Green Subsample
S_G_Uncorrected = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleGSmall/DeltaSigmaFiles/Uncorrected_Calculated_DeltaSigma/UncorrectedEMJKSS.txt")
S_G_Boost = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleGSmall/Boost/Boost.txt")
S_G_PhotoZ = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleGSmall/PhotoZ/PhotoZ.txt")
S_G_Random = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleGSmall/Random/RandomEM.txt")

#Large I Subsample
L_I_Uncorrected = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleILarge/DeltaSigmaFiles/Uncorrected_Calculated_DeltaSigma/UncorrectedEMJKSS.txt")
L_I_Boost = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleILarge/Boost/Boost.txt")
L_I_PhotoZ = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleILarge/PhotoZ/PhotoZ.txt")
L_G_Random = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleILarge/Random/RandomEM.txt")

#Small I Subsample
S_I_Uncorrected = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleISmall/DeltaSigmaFiles/Uncorrected_Calculated_DeltaSigma/UncorrectedEMJKSS.txt")
S_I_Boost = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleISmall/Boost/Boost.txt")
S_I_PhotoZ = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleISmall/PhotoZ/PhotoZ.txt")
S_I_Random = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleISmall/Random/RandomEM.txt")

"""
CORRECT THE SIGNAL
"""
#For Large Red
DSLR_Boost = L_R_Uncorrected*L_R_Boost
DSLR_PhotoZ = L_R_Uncorrected*L_R_Boost*L_R_PhotoZ
np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Large/DSLR_Boost.txt", DSLR_Boost)
np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Large/DSLR_PhotoZ.txt", DSLR_PhotoZ)

"""CALUCLATE COVARIANCE"""
DSLR_BoostCov = CovarianceAndCorrelation.getCovariance(DSLR_Boost)
DSLR_PhotoZCov = CovarianceAndCorrelation.getCovariance(DSLR_PhotoZ)
np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Large/DSLR_BoostCov.txt", DSLR_BoostCov)
np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Large/DSLR_PhotoZCov.txt", DSLR_PhotoZCov)

#For Small Read
DSSR_Boost = S_R_Uncorrected*S_R_Boost
DSSR_PhotoZ = S_R_Uncorrected*S_R_Boost*S_R_PhotoZ
np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Small/DSSR_Boost.txt", DSSR_Boost)
np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Small/DSSR_PhotoZ.txt", DSSR_PhotoZ)

"""CALUCLATE COVARIANCE"""
DSSR_BoostCov = CovarianceAndCorrelation.getCovariance(DSSR_Boost)
DSSR_PhotoZCov = CovarianceAndCorrelation.getCovariance(DSSR_PhotoZ)
np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Small/DSSR_BoostCov.txt", DSSR_BoostCov)
np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Small/DSSR_PhotoZCov.txt", DSSR_PhotoZCov)

print "Red Correction Finished!"

#For Large Green
DSLG_Boost = L_G_Uncorrected*L_G_Boost
DSLG_PhotoZ = L_G_Uncorrected*L_G_Boost*L_G_PhotoZ
np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Large/DSLG_Boost.txt", DSLG_Boost)
np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Large/DSLG_PhotoZ.txt", DSLG_PhotoZ)

"""CALUCLATE COVARIANCE"""
DSLG_BoostCov = CovarianceAndCorrelation.getCovariance(DSLG_Boost)
DSLG_PhotoZCov = CovarianceAndCorrelation.getCovariance(DSLG_PhotoZ)
np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Large/DSLG_BoostCov.txt", DSLG_BoostCov)
np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Large/DSLG_PhotoZCov.txt", DSLG_PhotoZCov)

#For Small Green
DSSG_Boost = S_G_Uncorrected*S_G_Boost
DSSG_PhotoZ = S_G_Uncorrected*S_G_Boost*S_G_PhotoZ
np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Small/DSSG_Boost.txt", DSSG_Boost)
np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Small/DSSG_PhotoZ.txt", DSSG_PhotoZ)

"""CALUCLATE COVARIANCE"""
DSSG_BoostCov = CovarianceAndCorrelation.getCovariance(DSSG_Boost)
DSSG_PhotoZCov = CovarianceAndCorrelation.getCovariance(DSSG_PhotoZ)
np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Small/DSSG_BoostCov.txt", DSSG_BoostCov)
np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Small/DSSG_PhotoZCov.txt", DSSG_PhotoZCov)

print "Green Correction Finished!"

#For Large I
DSLI_Boost = L_I_Uncorrected*L_I_Boost
DSLI_PhotoZ = L_I_Uncorrected*L_I_Boost*L_I_PhotoZ
np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Large/DSLI_Boost.txt", DSLI_Boost)
np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Large/DSLI_PhotoZ.txt", DSLI_PhotoZ)

"""CALUCLATE COVARIANCE"""
DSLI_BoostCov = CovarianceAndCorrelation.getCovariance(DSLI_Boost)
DSLI_PhotoZCov = CovarianceAndCorrelation.getCovariance(DSLI_PhotoZ)
np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Large/DSLI_BoostCov.txt", DSLI_BoostCov)
np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Large/DSLI_PhotoZCov.txt", DSLI_PhotoZCov)

#For Small I
DSSI_Boost = S_I_Uncorrected*S_I_Boost
DSSI_PhotoZ = S_I_Uncorrected*S_I_Boost*S_I_PhotoZ
np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Small/DSSI_Boost.txt", DSSI_Boost)
np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Small/DSSI_PhotoZ.txt", DSSI_PhotoZ)

"""CALUCLATE COVARIANCE"""
DSSI_BoostCov = CovarianceAndCorrelation.getCovariance(DSSI_Boost)
DSSI_PhotoZCov = CovarianceAndCorrelation.getCovariance(DSSI_PhotoZ)
np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Small/DSSI_BoostCov.txt", DSSI_BoostCov)
np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Small/DSSI_PhotoZCov.txt", DSSI_PhotoZCov)

print "I Correction Finished!"

"""
***FOR B MODE***
"""