"""
CALCULATES THE DELTA SIGMA, AND NAIVE ERROR OF THE E AND B MODE
"""
import numpy as np
from astropy.io import fits as pyfits #Import astro py as this
import matplotlib.pyplot as plt
import pylab
pylab.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size':10})
pylab.rc('text', usetex=True)

def getDeltaSigma(A , B, Error = None):
    a = np.sum(A, axis=0)
    b = np.sum(B, axis=0)
    DeltaSigma = (0.5) * (1. / (1 - (0.365) ** 2)) * (a / b)
    if Error is None:
        return DeltaSigma
    else:
        error = np.sum(Error, axis=0)
        DeltaSigma_Error = (0.5) * (1. / (1 - (0.365) ** 2)) * (error**0.5/ b)
        return DeltaSigma, DeltaSigma_Error

def getDeltaSigmaEM(sums):
    return getDeltaSigma(sums['wetsigma'], sums['w'], sums['wetsigma_sq'])

def getDeltaSigmaBM(sums):
    return getDeltaSigma(sums['wexsigma'], sums['w'], sums['wexsigma_sq'])

#For  the 83 subsamples
def getDeltaSigmaJKSSEM(sums, jk_indices):
    DeltaSigma_Array = np.zeros((83, 24))
    for i in range(83):
        WetSigma = sums["wetsigma"]
        W = sums["w"]

        Need_to_Delete = np.where(jk_indices == i)[0] #Contains the indicies of the data points needed to be deleted
        New_WetSigma = np.delete(WetSigma, Need_to_Delete, 0) #Indicies that need to be corresponded with the sums.fits file
        New_W = np.delete(W, Need_to_Delete, 0)

        DeltaSigma_Array[i] = getDeltaSigma(New_WetSigma , New_W)
    print "JKSSEM Calculated"
    return DeltaSigma_Array

def getDeltaSigmaJKSSBM(sums, jk_indices):
    DeltaSigma_Array = np.zeros((83, 24))
    for i in range(83):
        WetSigma = sums["wexsigma"]
        W = sums["w"]

        Need_to_Delete = np.where(jk_indices == i)[0] #Contains the indicies of the data points needed to be deleted
        New_WetSigma = np.delete(WetSigma, Need_to_Delete, 0) #Indicies that need to be corresponded with the sums.fits file
        New_W = np.delete(W, Need_to_Delete, 0)

        DeltaSigma_Array[i] = getDeltaSigma(New_WetSigma , New_W)
    print "JKSSBM Calculated"
    return DeltaSigma_Array

def getDeltaSigmaSSBM(SumsFilePath, BoostAndPhotoZIndiciesPath):
    DeltaSigma = ((83, 24))
    for i in range(83):
        index = np.load("%s" %BoostAndPhotoZIndiciesPath)
        S = pyfits.getdata("%s" %SumsFilePath)
        WexSigma = S["wexsigma"]
        W = S["w"]

        Need_to_Delete = np.where(index == i)[0] #Contains the indicies of the data points needed to be deleted
        New_WexSigma = np.delete(WexSigma, Need_to_Delete, 0)#Indicies that need to be corresponded with the sums.fits file
        New_W = np.delete(W, Need_to_Delete, 0)
        np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/DeltaSigmaFiles/BWexSigma%s" % i,np.c_[New_WexSigma])
        np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/DeltaSigmaFiles/BW%s" % i,np.c_[New_W])
    print "SSBM Calculated"

"""
CALCULATES DELTA SIGMA
"""
#Arrays will store all the summed up values
def getDeltaSigmaFSSEM():
    DeltaSigma_Array = np.zeros((83, 24))
    for i in range(83):
        #Loads the data
        WetSigma = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/DeltaSigmaFiles/EWetSigma%s" % i)
        W = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/DeltaSigmaFiles/EW%s" % i)
        a = np.sum(WetSigma, axis = 0) #Sums up all the wetsigma
        b = np.sum(W, axis = 0) #Sums up all the w
        #Calculate Delta Sigma And Error and stores into DeltaSigma_Array
        DeltaSigma_Array[i] = (0.5) * (1. / (1 - (0.365) ** 2)) * (a / b)
    print "FSSEM calculated"
    return DeltaSigma_Array

def getDeltaSigmaFSSBM():
    DeltaSigma_Array = np.zeros((83, 24))
    for i in range(83):
        #Loads the data
        WexSigma = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/DeltaSigmaFiles/BWexSigma%s" % i)
        W = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/DeltaSigmaFiles/BW%s" % i)
        a = np.sum(WexSigma, axis = 0) #Sums up all the wetsigma
        b = np.sum(W, axis = 0) #Sums up all the w
        #Calculate Delta Sigma And Error and stores into DeltaSigma_Array
        DeltaSigma_Array[i] = (0.5) * (1. / (1 - (0.365) ** 2)) * (a / b)
    print "FSSBM Calculated"
    return DeltaSigma_Array