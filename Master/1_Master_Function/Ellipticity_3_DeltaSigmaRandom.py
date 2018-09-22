"""
CALCULATE THE DELTA SIGMA OF THE 100 RANDOM SAMPLES FOR THE LARGE AND SMALL SUBSAMPLES
"""

import numpy as np
from astropy.io import fits as pyfits #Import astro py as this

def getDeltaSigma(A , B):
    DeltaSigma = (0.5)*(1./(1-(0.365)**2))*(A/B)
    return DeltaSigma

def calculateDeltaSigmaEllipticityRandom():
    for i in range(100):
        """
        CALCULATES DELTA SIGMA FOR LARGE SUBSAMPLES
        """
        # Load the catalogs
        d_L = pyfits.getdata("/Users/emolina/PycharmProjects/JPL/Master_Function/Data_Files/gglens_rand_dr8_rand_zmask_redmapper_v5.10_randcat_z0.05-0.60_lgt020.Z_LAMBDA_0.1_0.33.LAMBDA_lt100.fit.lens/dr8_rand_zmask_redmapper_v5.10_randcat_z0.05-0.60_lgt020.Z_LAMBDA_0.1_0.33.LAMBDA_lt100.%s.fit/sums.fits" % i)
        f_L = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Data_Files/gglens_rand_dr8_rand_zmask_redmapper_v5.10_randcat_z0.05-0.60_lgt020.Z_LAMBDA_0.1_0.33.LAMBDA_lt100.fit.lens/dr8_rand_zmask_redmapper_v5.10_randcat_z0.05-0.60_lgt020.Z_LAMBDA_0.1_0.33.LAMBDA_lt100.%s.fit/deltaSigma.dat" % i).transpose()

        # Get X-Axis for the graph
        x = f_L[0]

        # Calculates Delta Sigma
        A_L = np.sum(d_L["wetsigma"][range(len(d_L["wetsigma"]))], axis=0)
        B_L = np.sum(d_L["w"][range(len(d_L["w"]))], axis=0)
        DeltaSigma_L = getDeltaSigma(A_L, B_L)

        """
        CALCULATES DELTA SIGMA FOR SMALL SUBSAMPLES
        """
        # Load the catalogs
        d_S = pyfits.getdata("/Users/emolina/PycharmProjects/JPL/Master_Function/Data_Files/gglens_rand_dr8_rand_zmask_redmapper_v5.10_randcat_z0.05-0.60_lgt020.Z_LAMBDA_0.1_0.33.LAMBDA_lt100.fit.lens/dr8_rand_zmask_redmapper_v5.10_randcat_z0.05-0.60_lgt020.Z_LAMBDA_0.1_0.33.LAMBDA_lt100.%s.fit/sums.fits" % i)
        f_S = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Data_Files/gglens_rand_dr8_rand_zmask_redmapper_v5.10_randcat_z0.05-0.60_lgt020.Z_LAMBDA_0.1_0.33.LAMBDA_lt100.fit.lens/dr8_rand_zmask_redmapper_v5.10_randcat_z0.05-0.60_lgt020.Z_LAMBDA_0.1_0.33.LAMBDA_lt100.%s.fit/deltaSigma.dat" % i).transpose()

        # Get X-Axis for the graph
        x = f_S[0]

        A_S = np.sum(d_S["wetsigma"][range(len(d_L["wetsigma"]))], axis=0)
        B_S = np.sum(d_S["w"][range(len(d_L["w"]))], axis=0)
        DeltaSigma_S = getDeltaSigma(A_S, B_S)

        print "Calculating DeltaSigma of Large and Small", i

        """
        SAVES THE DATA INTO NEW FILES
        """
        np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Ellipticity/Calculated_Delta_Sigma/Random_Large_Subsample/DeltaSigmaL%s.txt" % i, np.c_[x, DeltaSigma_L])
        np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Ellipticity/Calculated_Delta_Sigma/Random_Small_Subsample/DeltaSigmaS%s.txt" % i, np.c_[x, DeltaSigma_S])
