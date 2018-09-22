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
#For Calculating JackKnife
import JackKnife_Calculator
#For Photo-Z
import Photo_Z
import Photo_Z2

#For Sorting the Ellipticity
import Ellipticity_1_Sort
import Ellipticity_2_CalculateDeltaSigma
import Ellipticity_3_DeltaSigmaRandom
import Ellipticity_4_Correct_DeltaSigma

#Input the locations of the data ex: /Users/emoline/PyCharm...etc...
#BoostAndRandomInciciesPath = /Users/emolina/PycharmProjects/JPL/Master_Function/Data_Files/redmapper_dr8_public_v5.10_catalog.value_added.Z_LAMBDA_0.1_0.33.LAMBDA_lt100.fits.jk_indices.npy

def getSignalData(SumsFilePath, DatFilePath, BoostAndPhotoZIndiciesPath, CombinedEverythingFitsPath, RedMapperData, RedMapperOpen, High, Low, WeakLensingBool, BoostBool, RandomBool, Photo_ZBool, CorrectionBool, PlotBool, SortBool, EllipticityBool, CorrectEllipticityBool, PlotEllipticity):

    """
    LOAD THE DATA
    """
    #Loads the sums file
    Sums = pyfits.getdata(SumsFilePath)
    #Load the dat file containing data for setting the x-axis
    f = np.loadtxt(DatFilePath).transpose()
    x = f[0]

    if WeakLensingBool == True:
        """
        CALCULATE THE WEAK LENSING SIGNAL FOR THE UNCORRECTED PARENT SAMPLE & NAIVE ERRORS
        """
        #DeltaSigma_Uncorrected_EM, DeltaSigma_Uncorrected_EM_nError = DeltaSigma_Uncorrected_EandBMode.getDeltaSigmaEM('wetsigma', 'w', 'wetsigma_sq', Sums)
        #DeltaSigma_Uncorrected_BM, DeltaSigma_Uncorrected_BM_nError = DeltaSigma_Uncorrected_EandBMode.getDeltaSigmaBM('wexsigma', 'w', 'wexsigma_sq', Sums)

        """
        CALCULATE THE WEAK LENSING SIGNAL
        """
        #Parent Sample; Deletes indicies and saves 'wetsigma', 'wexsigma', and 'w' to a file
        DeltaSigma_Uncorrected_EandBMode.getDeltaSigmaSSEM(SumsFilePath, BoostAndPhotoZIndiciesPath)
        DeltaSigma_Uncorrected_EandBMode.getDeltaSigmaSSBM(SumsFilePath, BoostAndPhotoZIndiciesPath)

        #Calculates Delta Sigma from the files created above (Uses 'wetsigma', 'wexsigma', 'w')
        DeltaSigma_Uncorrected_EM_SS = DeltaSigma_Uncorrected_EandBMode.getDeltaSigmaFSSEM()
        DeltaSigma_Uncorrected_BM_SS = DeltaSigma_Uncorrected_EandBMode.getDeltaSigmaFSSBM()

        """
        SAVE THE CALCULATED DATA
        """
        np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/DeltaSigmaFiles/Uncorrected_Calculated_DeltaSigma/UncorrectedEM.txt", np.c_[DeltaSigma_Uncorrected_EM_SS])
        np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/DeltaSigmaFiles/Uncorrected_Calculated_DeltaSigma/UncorrectedBM.txt", np.c_[DeltaSigma_Uncorrected_BM_SS])

    else:
        print "Delta Sigma Calculation Skipped"

    if BoostBool == True:
        """
        CALCULATE THE BOOST OF 83 SUBSAMPLES; FOR THE UNCORRECTED LENSING SIGNAL
        """
        BoostArray, BoostAverage =  Boost_Subsamples.setBoostFiles(BoostAndPhotoZIndiciesPath, SumsFilePath) #Boost Average can also be called
        print BoostArray
        """CALCULATE JACKKNIFE OF THE BOOST CORRECTION"""
        #Boost_Subsample_JackKnife = JackKnife_Calculator.calcJackKnife(BoostArray)
        #Covariance_Boost_Subsample = CovarianceAndCorrelation.getCovariance(BoostArray)
        #Correlation_Boost_Subsample = CovarianceAndCorrelation.getCorrelation(Covariance_Boost_Subsample)

        """
        SAVE THE CALCULATED DATA
        """
        np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Boost/Boost_Array.txt", np.c_[BoostArray])
        #np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Boost/Covariance_Boost.txt", np.c_[Covariance_Boost_Subsample])
        #np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Boost/Correlation_Boost.txt", np.c_[Correlation_Boost_Subsample])
    else:
        print "Boost Skipped"

    if RandomBool == True:
        """
        CALCULATE RANDOM PSF
        """
        Random_ArrayEM = Random_PSF.setRandomFilesEM()
        Random_ArrayBM = Random_PSF.setRandomFilesBM()


        """CALCULATE JACKKNIFE, COVARIANCE, AND CORRELATION OF THE RANDOM CORRECTION; E MODE AND B MODE"""
        Random_Subsample_JackKnifeEM = JackKnife_Calculator.calcJackKnife(Random_ArrayEM)
        #Covariance_RandomEM = CovarianceAndCorrelation.getCovariance(Random_ArrayEM)
        #Correlation_RandomEM = CovarianceAndCorrelation.getCorrelation(Covariance_RandomEM)

        Random_Subsample_JackKnifeBM = JackKnife_Calculator.calcJackKnife(Random_ArrayBM)
        #Covariance_RandomBM = CovarianceAndCorrelation.getCovariance(Random_ArrayBM)
        #Correlation_RandomBM = CovarianceAndCorrelation.getCorrelation(Covariance_RandomBM)

        """
        SAVE THE DATA
        """

        #Saves the Random Arrays for E and B Mode
        np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Random/Random_EM.txt", np.c_[Random_ArrayEM])
        np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Random/Random_BM.txt", np.c_[Random_ArrayBM])
        """
        #Save the Random Covariance and Correlation
        np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Random/Covariance_Random_EM.txt", np.c_[Covariance_RandomEM])
        np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Random/Correlation_Random_EM.txt", np.c_[Correlation_RandomEM])

        np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Random/Covariance_Random_BM.txt", np.c_[Covariance_RandomBM])
        np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Random/Correlation_Random_BM.txt", np.c_[Correlation_RandomBM])
        """
    else:
        print "Random Skipped"

    if Photo_ZBool == True:
        """
        CALCULATE PHOTO-Z; AVG AND STDDEV
        """
        #Set Photo-Z Fits Files
        Photo_Z.setFitsFiles(BoostAndPhotoZIndiciesPath, CombinedEverythingFitsPath)
        Photo_Z83 = Photo_Z2.getPhotoZCorrection()
        #Calculate Photo-S
        Photo_Z_Array = np.array([Photo_Z83] * 24).transpose()  # Repeats the column 24 times >> axis = 1; without transpose axis = 0

        Photo_Z_Average = np.average(Photo_Z83)
        Photo_Z_StdDev = np.std(Photo_Z83)

        np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Photo_Z/Photo_Z_Array.txt", np.c_[Photo_Z_Array])
    else:
        print "Photo-Skipped"

    if CorrectionBool == True:
        """
        LOAD THE DATA
        """
        Boost_Array = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Boost/Boost_Array.txt")
        PhotoZ_Array = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Photo_Z/Photo_Z_Array.txt")
        Random_Array_EM = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Random/Random_EM.txt")
        Random_Array_BM = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Random/Random_BM.txt")

        Photo_ZArray = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Photo_Z/Photo_Z_Array.txt")

        DeltaSigma_EMSS = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/DeltaSigmaFiles/Uncorrected_Calculated_DeltaSigma/UncorrectedEM.txt")
        DeltaSigma_BMSS = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/DeltaSigmaFiles/Uncorrected_Calculated_DeltaSigma/UncorrectedBM.txt")

        """
        APPLY THE CORRECTIONS TO EMODE THEN BMODE
        """
        Corrected_EM = (Boost_Array * PhotoZ_Array * DeltaSigma_EMSS) - (Photo_ZArray * Random_Array_EM)
        Corrected_BM = DeltaSigma_BMSS - Random_Array_BM

        """SAVE THE DATA"""
        np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Corrected_Signals/Corrected_EM.txt", np.c_[Corrected_EM])
        np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Corrected_Signals/Corrected_BM.txt", np.c_[Corrected_BM])


        """
        CALCULATE THE COVARIANCE OF THE DATA
        """
        #Calculates the covariance and correlation fo the coreted E mode and B mode
        Covariance_CorrectedEM = CovarianceAndCorrelation.getCovariance(Corrected_EM)
        Correlation_CorrectedEM = CovarianceAndCorrelation.getCorrelation(Covariance_CorrectedEM)

        Covariance_CorrectedBM = CovarianceAndCorrelation.getCovariance(Corrected_BM)
        Correlation_CorrectedBM = CovarianceAndCorrelation.getCorrelation(Covariance_CorrectedBM)

        """SAVE THE COVARAINCE AND CORRELATION OF THE DATA"""
        np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Covariance_Correlation/Covariance_CorrectedEM.txt", np.c_[Covariance_CorrectedEM])
        np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Covariance_Correlation/Correlation_CorrectedEM.txt", np.c_[Correlation_CorrectedEM])

        np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Covariance_Correlation/Covariance_CorrectedBM.txt", np.c_[Covariance_CorrectedBM])
        np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Covariance_Correlation/Correlation_CorrectedBM.txt", np.c_[Correlation_CorrectedBM])



        ErrorBarEM = np.sqrt(np.diagonal(Covariance_CorrectedEM))
        ErrorBarBM = np.sqrt(np.diagonal(Covariance_CorrectedBM))


        EM_Correct_Avd = np.average(Corrected_EM, axis = 0)
        BM_Correct_Avd = np.average(Corrected_BM, axis=0)

        plt.figure(1)
        plt.errorbar(x, EM_Correct_Avd, yerr = ErrorBarEM)
        plt.semilogx()
        plt.semilogy()

        plt.figure(2)
        plt.errorbar(x, x*BM_Correct_Avd, yerr = ErrorBarBM*x)
        plt.semilogx()

        plt.show()
    else:
        print "Files never corrected"

    if PlotBool == True:
        #***NEED TO LOAD THE X, COVARIANCE, DELTASIGMA (CORRECTED), AND CORRELATION***
        plt.figure(1)
        plt.errorbar(x, EM_Correct_Avd)
    else:
        print "Weak lensing data never plotted"

    if SortBool == True:
        """
        SORT THE REDMAPPER DATA
        """
        Ellipticity_1_Sort.SortRedMapper(RedMapperData, RedMapperOpen, High, Low)
    else:
        print "Data Sorting Skipped. Either because it was already calculated or by accident."

    if EllipticityBool == True:
        """
        CALCULATE THE DELTA SIGMA OF THE ELLIPTICITY
        """
        Ellipticity_2_CalculateDeltaSigma.calculateDeltaSigma(SumsFilePath, DatFilePath)
        Ellipticity_3_DeltaSigmaRandom.calculateDeltaSigmaEllipticityRandom()
    else:
        print "Delta Sigma of the Ellipticity Skipped"

    if CorrectEllipticityBool == True:
        """
        CORRECT THE ELLIPTICITY DATA
        """
        Ellipticity_4_Correct_DeltaSigma.calculatePhotoZ()
        Ellipticity_4_Correct_DeltaSigma.correctDeltaSigmaE()
        print "Signals Corrected!"

    else:
        print "Ellipticity Never Corrected"

    if PlotEllipticity == True:
        """
        PLOT THE ELLIPTICITY DATA (HISTOGRAMS....ETC)
        """
    else:
        print "Ellipticity Data never corrected"

