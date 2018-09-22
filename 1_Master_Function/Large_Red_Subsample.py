"""
FOR THE LARGE RED SUBSAMPLE
"""
import numpy as np
from astropy.io import fits as pyfits #Import astro py as this
import matplotlib.pyplot as plt
import pylab
pylab.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size':10})
pylab.rc('text', usetex=True)

import Master_Function

"""
FOR CALCULATING DELTA SIGMA OF UNCORRECTED PARENT SAMPLE
"""
WeakLensingBool = True
BoostBool = True
RandomBool = True
Photo_ZBool = True
CorrectionBool = True

n_rand = 100
"""
DEFINE PATHS FOR CALCULATING DELTA SIGMA
"""
LensSumsPath = "/Users/emolina/PycharmProjects/JPL/Master_Function/Data_Files/sums.fits"
LensJKIndicesPath = "/Users/emolina/PycharmProjects/JPL/Master_Function/Data_Files/redmapper_dr8_public_v5.10_catalog.value_added.Z_LAMBDA_0.1_0.33.LAMBDA_lt100.fits.jk_indices.npy" #Indices for Boost and Photo_Z
RandomSumsPathBase = "/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Data_Files/gglens_rand_dr8_rand_zmask_redmapper_v5.10_randcat_z0.05-0.60_lgt020.Z_LAMBDA_0.1_0.33.czl3_1.fit.lens/dr8_rand_zmask_redmapper_v5.10_randcat_z0.05-0.60_lgt020.Z_LAMBDA_0.1_0.33.czl3_1.%s.fit/sums.fits"
RandomJKIndicesPathBase = "/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Data_Files/dr8_rand_zmask_redmapper_v5.10_randcat_z0.05-0.60_lgt020.Z_LAMBDA_0.1_0.33.czl3_1.fit.lens/dr8_rand_zmask_redmapper_v5.10_randcat_z0.05-0.60_lgt020.Z_LAMBDA_0.1_0.33.czl3_1.%s.fit.jk_indices.npy"
SubsampleIndicesPath = "/Users/emolina/PycharmProjects/JPL/BCG_E_LIndices.txt"
SubsampleColumn = 0
# define output directory name
OutputDirname = "../Output_Files/OutPutSubsampleRLarge"

"""
SET THE FILES
"""
# load subsample indices
SubsampleIndices = np.loadtxt(SubsampleIndicesPath, dtype=np.int).transpose()[SubsampleColumn]
# load lens data
LensSums = pyfits.getdata(LensSumsPath)
LensSums = LensSums[SubsampleIndices]
# LensSums[indices] <- indices = [0, 4, 5, 6,....]
LensJKIndices = np.load(LensJKIndicesPath)
LensJKIndices = LensJKIndices[SubsampleIndices]
# LensJKIndices = LensJKIndices[indices]

# load random data
l_RandomSums = list()
l_RandomJKIndices = list()
for i in range(n_rand):
    l_RandomSums.append(pyfits.getdata(RandomSumsPathBase % i))
    l_RandomJKIndices.append(np.load(RandomJKIndicesPathBase % i))


"""
DETINE PATH FOR THE COMBINED EVERYTHING SUMS FILE
"""
Combined_EverythingPath = "/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Data_Files/combined_everything.fits"
Combined_Everything = pyfits.getdata(Combined_EverythingPath)

Master_Function.getSignalData(LensSums, LensJKIndices, l_RandomSums, l_RandomJKIndices, OutputDirname, Combined_Everything,  WeakLensingBool, BoostBool, RandomBool, Photo_ZBool, CorrectionBool)

"""
CALCULATE DATA FOR THE CORRECTED SIGNAL
"""

