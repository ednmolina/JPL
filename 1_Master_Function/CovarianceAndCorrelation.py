"""
CALCULATES COVARIANCE AND CORRELATION
"""
import numpy as np
from astropy.io import fits as pyfits #Import astro py as this
import matplotlib.pyplot as plt
import pylab
pylab.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size':10})
pylab.rc('text', usetex=True)

def getCovariance(AllData):
    #Calculates Average
    AllData_Avg = np.average(AllData, axis = 0)
    n = len(AllData)
    AllData_Minus_Average = np.zeros(AllData.shape)
    for i in range(len(AllData)):
        AllData_Minus_Average[i] = AllData[i] - AllData_Avg
    Covariance = (n-1.)/(n)*np.dot(AllData_Minus_Average.transpose(), AllData_Minus_Average)
    return Covariance

#Calculates Correlation

def getCorrelation(Covariance):
    Correlation = np.zeros(Covariance.shape)
    for i in range(len(Covariance)):
        for j in range(len(Covariance)):
            Correlation[i, j] = Covariance[i, j]/np.sqrt(Covariance[i, i]*Covariance[j, j])
    return Correlation