"""
2
PROGRAM CALCULATES THE DELTA SIGMA T HE DELTA SIGMA OF THE LARGE AND SMALL ELLIPTICITY SUBSAMPLES
"""

import numpy as np
from astropy.io import fits as pyfits #Import astro py as this
import matplotlib.pyplot as plt
import pylab
pylab.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size':10})
pylab.rc('text', usetex=True)

#Calculates the Delta Sigma, Error, and Numberator for Boost Calculation
def getDeltaSigma(name_col_A , name_col_B , name_col_error , Array_indicies , sums_file):
    a = np.sum(sums_file['%s' %name_col_A][Array_indicies] , axis = 0)
    b = np.sum(sums_file['%s' %name_col_B][Array_indicies], axis=0)
    DeltaSigma = (0.5) * (1. / (1 - (0.365) ** 2)) * (a / b)

    error = np.sum(sums_file['%s' %name_col_error][Array_indicies], axis = 0)
    error = ((0.5) / (1 - (0.365) ** 2) * (error**.5)) / b

    numerator = np.sum(sums_file["%s" %name_col_B], axis = 0)/len(sums_file["%s" %name_col_B])
    return DeltaSigma,error, numerator

#Calculates the Redshift and Richness
def getSubsampleColumn(col_name , Array_indicies , lens):
    col = lens['%s' % col_name][Array_indicies]
    return col

"""EDIT SO THAT IT TAKES IN THE SUMS & DATT FILE ARRAY"""
def calculateDeltaSigma(SumsFilePath, DatFilePath):
    c = pyfits.getdata("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Ellipticity/combined_everything.fits")

    # Gets the optical richness (lambda), redshift (z-lambda), BCG ellipticity (deVAB[r,g,i])
    richness = c["LAMBDA"]
    redshift = c["Z_LAMBDA"]
    BCG_r = c["deVAB_r"]
    BCG_g = c["deVAB_g"]
    BCG_i = c["deVAB_i"]
    RA = c["RA"]
    DEC = c["DEC"]

    """
    CALCULATE THE MEDIAN BCG_R_G_I
    """
    median_BCG_r = np.median(BCG_r)
    median_BCG_g = np.median(BCG_g)
    median_BCG_i = np.median(BCG_i)

    # Sorts the BCG R
    S_large_er = np.where(BCG_r > median_BCG_r)[0]
    S_small_er = np.where(BCG_r < median_BCG_r)[0]

    # Sorts the BCG G
    S_large_eg = np.where(BCG_g > median_BCG_g)[0]
    S_small_eg = np.where(BCG_g < median_BCG_g)[0]

    # Sorts the BCG I
    S_large_ei = np.where(BCG_i > median_BCG_i)[0]
    S_small_ei = np.where(BCG_i < median_BCG_i)[0]

    """
    SAVE INDICIES TO FILES
    """
    np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Ellipticity/BCG_Ellipticity_LIndices.txt", np.c_[S_large_er, S_large_eg, S_large_ei], fmt='%i')
    np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Ellipticity/BCG_Ellipticity_SIndices.txt", np.c_[S_small_er, S_small_eg, S_small_ei], fmt='%i')

    # Gets data from the sums.fits files
    S = pyfits.getdata("%s" %SumsFilePath)
    # Gets data for the x-axis
    f = np.loadtxt("%s" %DatFilePath).transpose()
    x = f[0]

    """
    CALCULATES DELTA SIGMA
    """
    # For the large sample
    DeltaSigma_L_R, error_L_r, numerator_L_r = getDeltaSigma('wetsigma', 'w', 'wetsigma_sq', S_large_er, S)
    DeltaSigma_L_G, error_L_g, numerator_L_g = getDeltaSigma('wetsigma', 'w', 'wetsigma_sq', S_large_eg, S)
    DeltaSigma_L_I, error_L_i, numerator_L_i = getDeltaSigma('wetsigma', 'w', 'wetsigma_sq', S_large_ei, S)

    # For the small sample
    DeltaSigma_S_R, error_S_r, numerator_S_r = getDeltaSigma('wetsigma', 'w', 'wetsigma_sq', S_small_er, S)
    DeltaSigma_S_G, error_S_g, numerator_S_g = getDeltaSigma('wetsigma', 'w', 'wetsigma_sq', S_small_er, S)
    DeltaSigma_S_I, error_S_i, numerator_S_i = getDeltaSigma('wetsigma', 'w', 'wetsigma_sq', S_small_er, S)

    """
    SAVE DELTA SIGMA TO A FILE
    """
    # New Parenthasis is a new column in the text file
    np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Ellipticity/Calculated_Delta_Sigma/BCG_E_DeltaSigma.txt", np.c_[x, DeltaSigma_L_R, DeltaSigma_L_G, DeltaSigma_L_I, DeltaSigma_S_R, DeltaSigma_S_G, DeltaSigma_S_I])
    np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Ellipticity/Calculated_Delta_Sigma/BCG_E_DeltaSigma_Error.txt", np.c_[error_L_r, error_L_g, error_L_i, error_S_r, error_S_g, error_S_i])
    np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Ellipticity/Calculated_Delta_Sigma/BCG_E_BoostNumerator.txt", np.c_[numerator_L_r, numerator_L_g, numerator_L_i, numerator_S_r, numerator_S_g, numerator_S_i])

    """
    GET DATA FOR THE HISTORGRAMS
    """
    # DATA FOR THE RICHNESS (Data then Median)
    richness_L_r = getSubsampleColumn("LAMBDA", S_large_er, c)
    richness_S_r = getSubsampleColumn("LAMBDA", S_small_er, c)
    richness_median_LR = np.median(richness_L_r)
    richness_median_SR = np.median(richness_S_r)

    richness_L_g = getSubsampleColumn("LAMBDA", S_large_eg, c)
    richness_S_g = getSubsampleColumn("LAMBDA", S_small_eg, c)
    richness_median_LG = np.median(richness_L_g)
    richness_median_SG = np.median(richness_S_g)

    richness_L_i = getSubsampleColumn("LAMBDA", S_large_ei, c)
    richness_S_i = getSubsampleColumn("LAMBDA", S_small_ei, c)
    richness_median_LI = np.median(richness_L_i)
    richness_median_SI = np.median(richness_S_i)

    # DATA FOR THE REDSHIFT (Data then Median)
    redshift_L_r = getSubsampleColumn("Z_LAMBDA", S_large_er, c)
    redshift_S_r = getSubsampleColumn("Z_LAMBDA", S_small_er, c)
    redshift_median_LR = np.median(redshift_L_r)
    redshift_median_SR = np.median(redshift_S_r)

    redshift_L_g = getSubsampleColumn("Z_LAMBDA", S_large_eg, c)
    redshift_S_g = getSubsampleColumn("Z_LAMBDA", S_small_eg, c)
    redshift_median_LG = np.median(redshift_L_g)
    redshift_median_SG = np.median(redshift_S_g)

    redshift_L_i = getSubsampleColumn("Z_LAMBDA", S_large_ei, c)
    redshift_S_i = getSubsampleColumn("Z_LAMBDA", S_small_ei, c)
    redshift_median_LI = np.median(redshift_L_i)
    redshift_median_SI = np.median(redshift_S_i)

    """
    SAVE THE HISTORGRAM DATA TO FILE
    """
    np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Ellipticity/Calculated_Delta_Sigma/BCG_E_Richness_Redshift.txt", np.c_[richness_L_r, richness_L_g, richness_L_i, richness_S_r, richness_S_g, richness_S_i, redshift_L_r, redshift_L_g, redshift_L_i, redshift_S_r, redshift_S_g, redshift_S_i])

    """
    SAVE THE MEDIAN S TO A FILE
    """
    np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Ellipticity/Calculated_Delta_Sigma/BCG_E_Median.txt", np.c_[richness_median_LR, richness_median_LG, richness_median_LI, richness_median_SR, richness_median_SG, richness_median_SI, redshift_median_LR, redshift_median_LG, redshift_median_LI, redshift_median_SR, redshift_median_SG, redshift_median_SI])