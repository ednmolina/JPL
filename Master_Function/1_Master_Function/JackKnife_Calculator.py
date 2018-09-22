"""
CALCULATES THE JACKKNFIE OF THE 83 SUBSAMPLES
"""
import numpy as np
from astropy.io import fits as pyfits #Import astro py as this
import matplotlib.pyplot as plt
import pylab
pylab.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size':10})
pylab.rc('text', usetex=True)

#Function that calculates JackKnife for Boost
def calcJackKnife(DataFile):
    #Calculate the Average of the Data File used to calculate JackKnife Error (Variation; sqrt of this gives Standard Deviation)
    DataFile_Avg = np.average(DataFile, axis = 0)
    #n = len(DataFile.transpose())
    n = 83
    JackKnife = np.zeros((len(DataFile), 24))
    for i in range(n):
        JackKnife[i] = ((n-1.) / (n)) * (DataFile[i]-DataFile_Avg)**2
    JackKnife_Sum = np.sum(JackKnife, axis = 0)
    return  JackKnife
