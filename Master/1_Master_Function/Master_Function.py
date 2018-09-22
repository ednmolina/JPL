"""
MASTER FUNCTION CALCULATES THE:
    WEAK LENSING UNCORRECTED SIGNAL OF THE PARENT SAMPLE (E AND B MODE)
        COVARAINCE OF UNCORRECTED WEAK LENSING SIGNAL (E AND B MODE)
    WEAK LENSING UNCORRECTED SIGNAL OF THE 83 JACK KNIFE SUBSAMPLE
        COVARAINE OF THE UNCORRECTED SIGNAL OF THE 83 JACKNIFE SUBSAMPLES

    BOOST OF THE PARENTS SAMPLE (E MODE)
        COVARIANCE AND CORRELATION OF THE BOOST CORRECTION
    BOOST OF THE 83 JACK KNIFE SUBSAMPLES
        COVARIANCE AND CORRELATION OF THE BOOST OF THE 83 SUBSAMPLES

    RANDOM PSF OF THE PARENT SAMPLES (E AND B MODE)
        COVARIANCE AND CORRELATION OF THE RANDOM PSF CORRECTION (E AND B MODE)
    RANDOM PSF OF THE 83 JACK KNIFE SUBSAMPLEDS
        COVARIANCE AND CORRELATION OF THE 83 RANDOM JACKKNFIE SUBSAMPLES

    CALCULATES THE PHOTOS Z OF THE UNCORRECTED PARENT SAMPLE
        COVARIANCE AND THE CORRELATION OF THE PHOTO Z UNCORRECTED PARENT SAMPLE
    CALCULATES THE PHOTO Z OF THE 83 JACKKNIFE SUBSAMPLES
        COVARAIANCE AND THE CORRELATION OF THE PHOTO Z OF THE 83 JACKKNIFE SUBSAMPLES

    CORRECTED WEAK LENSING SIGNAL (PARENT SAMPLE; E AND B MODE)
        COVARIANCE AND CORRELATION OF THE CORRECTED WEAK LENSING SAMPLES( E AND B MODE)
    CORRECTED WEAK LENSING SIGNAL OF THE 83 JACKKNIFE SUBSAMPLES
        COVARIANCE AND CORRELATION OF THE 83 JACKKNIFE SUBSAMPLES
"""

import numpy as np
import pylab
from astropy.io import fits as pyfits #Import astro py as this
import cosmology
import matplotlib.pyplot as plt
import os
pylab.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size':10})
pylab.rc('text', usetex=True)
#For the Uncorrected Weak Lensing Signal
import DeltaSigma_Uncorrected_EandBMode
#For Calculating Covariance and Corerelation
import CovarianceAndCorrelation
#For Calculating Boost of Subsamples
import Boost_Subsamples
#For Calculating Random
import Random_PSF
#For Photo-Z
import Photo_Z
import Photo_Z2

#Input the locations of the data ex: /Users/emoline/PyCharm...etc...
#BoostAndRandomInciciesPath = /Users/emolina/PycharmProjects/JPL/Master_Function/Data_Files/redmapper_dr8_public_v5.10_catalog.value_added.Z_LAMBDA_0.1_0.33.LAMBDA_lt100.fits.jk_indices.npy

def getSignalData(LensSums, LensJKIndices, l_RandomSums, l_RandomJKIndices, OutputDirname, Combined_Everything,WeakLensingBool = True, BoostBool = True, RandomBool = True, Photo_ZBool = True, CorrectionBool = True, n_rand = 100):

    # setup directory
    if os.path.exists(OutputDirname):
        Input = input("Would you like to overwrite the %s?" %OutputDirname)
        if Input == 'y':
            pass
        elif Input == 'n':
            print "Terminated by user."
            exit(1)
        else:
            print "Please enter y/n"
            exit(1)
    else:
        print "%s does not exit. Create a new directory." % OutputDirname
        os.mkdir(OutputDirname)

    """
    SET THE DIRECTORY NAMES
    """
    #FOR UNCORRECTED DELTA SIGMA
    OutputDirnameUncorrectedDeltaSigma = OutputDirname + "/DeltaSigmaFiles/Uncorrected_Calculated_DeltaSigma"
    #FOR BOOST
    OutputDirnameBoost = OutputDirname + "/Boost"
    #RANDOM SIGNAL
    OutputDirnameRandom = OutputDirname + "/Random"
    #PHOTO-Z CORRECTION
    OutputDirnamePhotoZ = OutputDirname + "/PhotoZ"
    #COVARIANCE AND CORRELATION
    OutputDirnameCovariance = OutputDirname + "/CovarianceAndCorrelation"
    #JACKKNIFE
    OutputDirnameJackKnife = OutputDirname + "/JackKnife"


    if WeakLensingBool == True:
        """
        CALCULATE THE WEAK LENSING SIGNAL FOR THE UNCORRECTED PARENT SAMPLE & NAIVE ERRORS
        """
        DeltaSigma_Uncorrected_EM, DeltaSigma_Uncorrected_EM_nError = DeltaSigma_Uncorrected_EandBMode.getDeltaSigmaEM(LensSums)
        DeltaSigma_Uncorrected_BM, DeltaSigma_Uncorrected_BM_nError = DeltaSigma_Uncorrected_EandBMode.getDeltaSigmaBM(LensSums)

        """
        SAVE THE CALCULATED DATA
        """

        if os.path.exists(OutputDirnameUncorrectedDeltaSigma) is False:
            os.makedirs(OutputDirnameUncorrectedDeltaSigma)
        np.savetxt(OutputDirnameUncorrectedDeltaSigma + "/UncorrectedEM.txt", np.c_[DeltaSigma_Uncorrected_EM, DeltaSigma_Uncorrected_EM_nError])
        np.savetxt(OutputDirnameUncorrectedDeltaSigma + "/UncorrectedBM.txt", np.c_[DeltaSigma_Uncorrected_BM, DeltaSigma_Uncorrected_BM_nError])

        # Calculates Delta Sigma from the files created above (Uses 'wetsigma', 'wexsigma', 'w'); E and B Mode
        DeltaSigma_Uncorrected_EM_JKSS = DeltaSigma_Uncorrected_EandBMode.getDeltaSigmaJKSSEM(LensSums, LensJKIndices)
        DeltaSigma_Uncorrected_BM_JKSS = DeltaSigma_Uncorrected_EandBMode.getDeltaSigmaJKSSBM(LensSums, LensJKIndices)

        #Calculate Covariance of Delta Sigma (Uncorrected)
        DeltaSigma_Uncorrected_EM_JKCov = CovarianceAndCorrelation.getCovariance(DeltaSigma_Uncorrected_EM_JKSS)
        DeltaSigma_Uncorrected_BM_JKCov = CovarianceAndCorrelation.getCovariance(DeltaSigma_Uncorrected_BM_JKSS)

        np.savetxt(OutputDirnameUncorrectedDeltaSigma + "/UncorrectedEMJKSS.txt", DeltaSigma_Uncorrected_EM_JKSS)
        np.savetxt(OutputDirnameUncorrectedDeltaSigma + "/UncorrectedBMJKSS.txt", DeltaSigma_Uncorrected_BM_JKSS)
        np.savetxt(OutputDirnameUncorrectedDeltaSigma + "/UncorrectedEMJKCov.txt", DeltaSigma_Uncorrected_EM_JKCov)
        np.savetxt(OutputDirnameUncorrectedDeltaSigma + "/UncorrectedBMJKCov.txt", DeltaSigma_Uncorrected_BM_JKCov)
    else:
        print "Delta Sigma Calculation Skipped"

    if BoostBool == True:
        """
        CALCULATE THE BOOST
        """
        Boost_JackKnife, Average_Boost = Boost_Subsamples.setBoostFiles(LensJKIndices, LensSums, l_RandomSums, l_RandomJKIndices)

        """ REMEMBER TO ADD THIS
        CALCULATE COVARIANCE OF BOOST
        """

        if os.path.exists(OutputDirnameBoost) is False:
            os.makedirs(OutputDirnameBoost)
        np.savetxt(OutputDirnameBoost + "/Boost.txt", np.c_[Boost_JackKnife])

    else:
        print "Boost Skipped"

    if RandomBool == True:
        """
        CALCULATE RANDOM PSF
        """
        Random_ArrayEM = Random_PSF.setRandomFiles(l_RandomSums, l_RandomJKIndices, "wetsigma")
        Random_ArrayBM = Random_PSF.setRandomFiles(l_RandomSums, l_RandomJKIndices, "wexsigma")

        """
        SAVE THE DATA
        """

        if os.path.exists(OutputDirnameRandom) is False:
            os.makedirs(OutputDirnameRandom)
        np.savetxt(OutputDirnameRandom + "/RandomEM.txt", np.c_[Random_ArrayEM])
        np.savetxt(OutputDirnameRandom + "/RandomBM.txt", np.c_[Random_ArrayBM])

    else:
        print "Random Skipped"

    if Photo_ZBool == True:
        """
        CALCULATE PHOTO-Z; AVG AND STDDEV
        """
        if os.path.exists(OutputDirnamePhotoZ) is False:
            os.makedirs(OutputDirnamePhotoZ)

        #Set Photo-Z Fits Files
        Photo_Z.setFitsFiles(LensJKIndices, Combined_Everything, OutputDirnamePhotoZ)
        Photo_Z83 = Photo_Z2.getPhotoZCorrection(OutputDirnamePhotoZ)
        #Calculate Photo-S
        Photo_Z_Array = np.array([Photo_Z83] * 24).transpose()  # Repeats the column 24 times >> axis = 1; without transpose axis = 0

        np.savetxt(OutputDirnamePhotoZ + "/PhotoZ.txt", np.c_[Photo_Z_Array])
    else:
        print "Photo-Skipped"

    if CorrectionBool == True:
        """
        LOAD THE CALCUALTED DATA; UNCORRECTED DELTA SIGMA, BOOST, RANDOM, AND PHOTO-Z
        """
        DeltaSigma_Uncorrected_EM_JKSS = np.loadtxt(OutputDirnameUncorrectedDeltaSigma + "/UncorrectedEMJKSS.txt")
        DeltaSigma_Uncorrected_BM_JKSS =np.loadtxt(OutputDirnameUncorrectedDeltaSigma + "/UncorrectedBMJKSS.txt")

        Boost_JackKnife = np.loadtxt(OutputDirnameBoost + "/Boost.txt")

        Random_ArrayEM = np.loadtxt(OutputDirnameRandom + "/RandomEM.txt")
        Random_ArrayBM = np.loadtxt(OutputDirnameRandom + "/RandomBM.txt")

        Photo_Z_Array = np.loadtxt(OutputDirnamePhotoZ + "/PhotoZ.txt")

        Corrected_EM = (Boost_JackKnife*Photo_Z_Array*DeltaSigma_Uncorrected_EM_JKSS) - (Photo_Z_Array*Random_ArrayEM)
        Corrected_BM = DeltaSigma_Uncorrected_BM_JKSS - Random_ArrayBM

        OutputDirnameDeltaSigma = OutputDirname + "/CorrectedDeltaSigma"
        if os.path.exists(OutputDirnameDeltaSigma) is False:
            os.makedirs(OutputDirnameDeltaSigma)
        np.savetxt(OutputDirnameDeltaSigma + "/CorrectedEM.txt", np.c_[Corrected_EM])
        np.savetxt(OutputDirnameDeltaSigma + "/CorrectedBM.txt", np.c_[Corrected_BM])

        """
        CALCULATE COVARIANCE AND CORRELATION OF THE DATA
        """
        Covariance_CorrectedDeltaSigmaEM = CovarianceAndCorrelation.getCovariance(Corrected_EM)
        Covariance_CorrectedDeltaSigmaBM = CovarianceAndCorrelation.getCovariance(Corrected_BM)

        Correlation_CorrectedDeltaSigmaEM = CovarianceAndCorrelation.getCorrelation(Covariance_CorrectedDeltaSigmaEM)
        Correlation_CorrectedDeltaSigmaBM = CovarianceAndCorrelation.getCorrelation(Covariance_CorrectedDeltaSigmaBM)

        if os.path.exists(OutputDirnameCovariance) is False:
            os.makedirs(OutputDirnameCovariance)
        np.savetxt(OutputDirnameCovariance + "/Covariance_CorrectedDeltaSigmaEM.txt", np.c_[Covariance_CorrectedDeltaSigmaEM])
        np.savetxt(OutputDirnameCovariance + "/Covariance_CorrectedDeltaSigmaBM.txt", np.c_[Covariance_CorrectedDeltaSigmaBM])
        np.savetxt(OutputDirnameCovariance + "/Correlation_CorrectedDeltaSigmaEM.txt", np.c_[Correlation_CorrectedDeltaSigmaEM])
        np.savetxt(OutputDirnameCovariance + "/Correlation_CorrectedDeltaSigmaBM.txt", np.c_[Correlation_CorrectedDeltaSigmaBM])

    else:
        print "Files never corrected"