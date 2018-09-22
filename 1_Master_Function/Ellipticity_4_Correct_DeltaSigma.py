"""
4
CORRECTED THE DELTA SIGMA OF THE ELLIPTICITY SUBSAMPLES
"""
import numpy as np
import pylab
from astropy.io import fits as pyfits #Import astro py as this
import cosmology
import matplotlib.pyplot as plt
pylab.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size':10})
pylab.rc('text', usetex=True)

"""
CALCULATE THE PHTOO-Z FIRST
"""

#Calculates the Redshift and Richness
def getSubsampleColumn(col_name , Array_indicies , lens):
    col = lens['%s' %col_name][Array_indicies]
    return col

def run(lens):
    d = np.genfromtxt("/Users/emolina/PycharmProjects/JPL/getbias.finalcat.all_nocosmos.out", dtype = [ ("z_lens", np.float64), ("col2", np.int32), ("bz", np.float64), ("sum_wk", np.float64), ("inv_Sigma_c", np.float64)])
    uni = cosmology.universe(Omega_m0 = 0.25, Omega_l0 = 0.75, H0=100.)
    w = uni.DA(d["z_lens"]) **(-2.) * (1.+d["z_lens"])**(-2.) * d["sum_wk"]

    # make bins f
    bin_edges = np.linspace(0.005, 0.595, 60)
    assign = np.digitize(lens["Z_LAMBDA"], bin_edges)
    n_lens = [np.sum(assign == i) for i in range(1, 60)]
    bz = np.sum(n_lens*w*d["bz"])/np.sum(n_lens*w)
    return bz

#ACTUALLY CALCULATE PHOTO-Z
def calculatePhotoZ():
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
    IMPORT INDICIES
    """
    L = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Ellipticity/BCG_Ellipticity_LIndices.txt", int).transpose()
    S = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Ellipticity/BCG_Ellipticity_SIndices.txt", int).transpose()

    """
    SET THE DATA TO ARRAYS
    """
    richness_LR = getSubsampleColumn("LAMBDA", L[0], c)
    richness_LG = getSubsampleColumn("LAMBDA", L[1], c)
    richness_LI = getSubsampleColumn("LAMBDA", L[2], c)

    redshift_SR = getSubsampleColumn("Z_LAMBDA", S[0], c)
    redshift_SG = getSubsampleColumn("Z_LAMBDA", S[1], c)
    redshift_SI = getSubsampleColumn("Z_LAMBDA", S[2], c)

    """
    CREATES NEW FITS FILES
    """
    # RED BCG
    new_cols = pyfits.ColDefs([pyfits.Column(name='LAMBDA', format='D', array=richness_LR),
                               pyfits.Column(name='Z_LAMBDA', format='D', array=redshift_SR)])
    hdu = pyfits.BinTableHDU.from_columns(new_cols)
    hdu.writeto('/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Ellipticity/Photo_Z_Fits/BCG_Ellip_R.fits')

    # Green BCG
    new_cols = pyfits.ColDefs([pyfits.Column(name='LAMBDA', format='D', array=richness_LG),
                               pyfits.Column(name='Z_LAMBDA', format='D', array=redshift_SG)])
    hdu = pyfits.BinTableHDU.from_columns(new_cols)
    hdu.writeto('/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Ellipticity/Photo_Z_Fits/BCG_Ellip_G.fits')

    # IR BCG
    new_cols = pyfits.ColDefs([pyfits.Column(name='LAMBDA', format='D', array=richness_LI),
                               pyfits.Column(name='Z_LAMBDA', format='D', array=redshift_SI)])
    hdu = pyfits.BinTableHDU.from_columns(new_cols)
    hdu.writeto('/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Ellipticity/Photo_Z_Fits/BCG_Ellip_I.fits')

    """
    CALCULATE THE PHOTO-Z OF THE THREE SAMPLES (RED, GREEN, INFERRED)
    """
    # Calculate Photo-Z Red
    lens = pyfits.getdata('/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Ellipticity/Photo_Z_Fits/BCG_Ellip_R.fits')
    bz = run(lens)
    correction_factorBCGR = 1. / (1. + bz)

    # Calculate Photo-Z Green
    lens = pyfits.getdata('/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Ellipticity/Photo_Z_Fits/BCG_Ellip_G.fits')
    bz = run(lens)
    correction_factorBCGG = 1. / (1. + bz)

    # Calculate Photo-Z I
    lens = pyfits.getdata('/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Ellipticity/Photo_Z_Fits/BCG_Ellip_I.fits')
    bz = run(lens)
    correction_factorBCGI = 1. / (1. + bz)

    """USE THESE VALUES TO CORRECT THE SIGNAL NEXT"""
    np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Ellipticity/Photo_Z_Fits/PhotoZ_BCG_RGI.txt", np.c_[correction_factorBCGR, correction_factorBCGG, correction_factorBCGI])

"""
CORRECT THE DETLA SIGMA
"""
def correctDeltaSigmaE():
    """
    IMPORT THE DELTA SIGMA AND ERROR (FOR ERROR BAR) FROM THE FILES
    """
    DeltaSigma = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Ellipticity/Calculated_Delta_Sigma/BCG_E_DeltaSigma.txt").transpose()
    Error = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Ellipticity/Calculated_Delta_Sigma/BCG_E_DeltaSigma_Error.txt").transpose()
    PhotoZ = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Ellipticity/Photo_Z_Fits/PhotoZ_BCG_RGI.txt").transpose()
    """
    CALCULATE BOOST SIGNALS
    """
    # Calculates numerator for boost
    numerators = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Ellipticity/Calculated_Delta_Sigma/BCG_E_BoostNumerator.txt").transpose()
    # Numerator for Large Boost
    numerator_LR = numerators[0]
    numerator_LG = numerators[1]
    numerator_LI = numerators[2]

    # Numerator for Small Boost
    numerator_SR = numerators[3]
    numerator_SG = numerators[4]
    numerator_SI = numerators[5]

    # Stores Boost of Subsamples
    realtns_LR = np.zeros((100, 24))
    realtns_LG = np.zeros((100, 24))
    realtns_LI = np.zeros((100, 24))
    realtns_SR = np.zeros((100, 24))
    realtns_SG = np.zeros((100, 24))
    realtns_SI = np.zeros((100, 24))

    # Stores large subsample
    all_L = np.zeros((100, 24))
    # Stores small subsample
    all_S = np.zeros((100, 24))

    # Calculate Denomenator
    for i in range(100):
        # For the large subsample
        d_L = pyfits.getdata("/Users/emolina/PycharmProjects/JPL/Master_Function/Data_Files/gglens_rand_dr8_rand_zmask_redmapper_v5.10_randcat_z0.05-0.60_lgt020.Z_LAMBDA_0.1_0.33.LAMBDA_lt100.fit.lens/dr8_rand_zmask_redmapper_v5.10_randcat_z0.05-0.60_lgt020.Z_LAMBDA_0.1_0.33.LAMBDA_lt100.%s.fit/sums.fits" % i)
        f_L = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Data_Files/gglens_rand_dr8_rand_zmask_redmapper_v5.10_randcat_z0.05-0.60_lgt020.Z_LAMBDA_0.1_0.33.LAMBDA_lt100.fit.lens/dr8_rand_zmask_redmapper_v5.10_randcat_z0.05-0.60_lgt020.Z_LAMBDA_0.1_0.33.LAMBDA_lt100.%s.fit/deltaSigma.dat" % i).transpose()

        "Boosts"
        denomenatorLR = np.sum(d_L["w"], axis=0) / len(d_L["w"])
        BoostLR = numerator_LR / denomenatorLR
        realtns_LR[i] = BoostLR

        denomenatorLG = np.sum(d_L["w"], axis=0) / len(d_L["w"])
        BoostLG = numerator_LG / denomenatorLG
        realtns_LG[i] = BoostLG

        denomenatorLI = np.sum(d_L["w"], axis=0) / len(d_L["w"])
        BoostLI = numerator_LI / denomenatorLI
        realtns_LI[i] = BoostLI

        # For the small subsamples
        d_S = pyfits.getdata("/Users/emolina/PycharmProjects/JPL/Master_Function/Data_Files/gglens_rand_dr8_rand_zmask_redmapper_v5.10_randcat_z0.05-0.60_lgt020.Z_LAMBDA_0.1_0.33.LAMBDA_lt100.fit.lens/dr8_rand_zmask_redmapper_v5.10_randcat_z0.05-0.60_lgt020.Z_LAMBDA_0.1_0.33.LAMBDA_lt100.%s.fit/sums.fits" % i)
        f_S = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Data_Files/gglens_rand_dr8_rand_zmask_redmapper_v5.10_randcat_z0.05-0.60_lgt020.Z_LAMBDA_0.1_0.33.LAMBDA_lt100.fit.lens/dr8_rand_zmask_redmapper_v5.10_randcat_z0.05-0.60_lgt020.Z_LAMBDA_0.1_0.33.LAMBDA_lt100.%s.fit/deltaSigma.dat" % i).transpose()

        "Boosts"
        denomenatorSR = np.sum(d_S["w"], axis=0) / len(d_S["w"])
        BoostSR = numerator_SR / denomenatorSR
        realtns_SR[i] = BoostSR

        denomenatorSG = np.sum(d_S["w"], axis=0) / len(d_S["w"])
        BoostSG = numerator_SG / denomenatorSG
        realtns_SG[i] = BoostSG

        denomenatorSI = np.sum(d_S["w"], axis=0) / len(d_S["w"])
        BoostSI = numerator_SI / denomenatorSI
        realtns_SI[i] = BoostSI

        """
        PSF CORRECTION
        """
        """
        LOAD CALCULATED RANDOM DELTA SIGMAS FOR LARGE SUBSAMPLE
        """
        DSR_L = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Ellipticity/Calculated_Delta_Sigma/Random_Large_Subsample/DeltaSigmaL%s.txt" % i).transpose()
        DSR_S = np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Ellipticity/Calculated_Delta_Sigma/Random_Small_Subsample/DeltaSigmaS%s.txt" % i).transpose()
        # Store Delta Sigma into array
        all_L[i] = DSR_L[1]  # Delta Sigma Random Large
        all_S[i] = DSR_S[1]
        print "Calculating Boost For Ellipticity ", i

        """
        CALCULATING THE AVERAGE (BOOST)
        """
        BoostCorr_LR = np.sum(realtns_LR, axis=0) / 100
        BoostCorr_LG = np.sum(realtns_LG, axis=0) / 100
        BoostCorr_LI = np.sum(realtns_LI, axis=0) / 100

        BoostCorr_SR = np.sum(realtns_SR, axis=0) / 100
        BoostCorr_SG = np.sum(realtns_SG, axis=0) / 100
        BoostCorr_SI = np.sum(realtns_SI, axis=0) / 100

        """
        CALCULATE AVERAGE OF PSF
        """
        avg_PSF_L = np.average(all_L, axis=0)

        avg_PSF_S = np.average(all_S, axis=0)

        """
        CORRECTED SIGNALS AND ERROR
        """
        # Only Bosst Correction Applied
        Corr_BST_LR = DeltaSigma[1] * BoostCorr_LR
        Corr_BST_LG = DeltaSigma[2] * BoostCorr_LG
        Corr_BST_LI = DeltaSigma[3] * BoostCorr_LI

        Corr_BST_SR = DeltaSigma[4] * BoostCorr_SR
        Corr_BST_SG = DeltaSigma[5] * BoostCorr_SG
        Corr_BST_SI = DeltaSigma[6] * BoostCorr_SI

        # Boost and Photo-Z
        """***EDIT THIS SECTION. REPLACE PHTOTO-Z WITH VARIABLES***"""
        """GET THE PHOTO Z"""
        RedZ = PhotoZ[0]
        GreenZ = PhotoZ[1]
        IZ = PhotoZ[2]

        Corr_BSTPZ_LR = DeltaSigma[1] * BoostCorr_LR * RedZ
        Corr_BSTPZ_LG = DeltaSigma[2] * BoostCorr_LG * GreenZ
        Corr_BSTPZ_LI = DeltaSigma[3] * BoostCorr_LI * IZ

        Corr_BSTPZ_SR = DeltaSigma[4] * BoostCorr_SR * RedZ
        Corr_BSTPZ_SG = DeltaSigma[5] * BoostCorr_SG * GreenZ
        Corr_BSTPZ_SI = DeltaSigma[6] * BoostCorr_SI * IZ

        # Corrects the signals
        Corr_LR = DeltaSigma[1] * BoostCorr_LR * RedZ - avg_PSF_L
        Corr_LG = DeltaSigma[2] * BoostCorr_LG * GreenZ - avg_PSF_L
        Corr_LI = DeltaSigma[3] * BoostCorr_LI * IZ - avg_PSF_L

        Corr_SR = DeltaSigma[4] * BoostCorr_SR * RedZ - avg_PSF_S
        Corr_SG = DeltaSigma[5] * BoostCorr_SG * GreenZ - avg_PSF_S
        Corr_SI = DeltaSigma[6] * BoostCorr_SI * IZ - avg_PSF_S

        # Corrects the error
        Error_LR = Error[0] * RedZ
        Error_LG = Error[1] * GreenZ
        Error_LI = Error[2] * IZ

        Error_SR = Error[3] * RedZ
        Error_SG = Error[4] * GreenZ
        Error_SI = Error[5] * IZ

        """
        SAVES THE DATA ONTO A FILE
        """
        # Save Corrected Delta Sigma onto a File
        np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Ellipticity/Calculated_Delta_Sigma/Corrected_DeltaSigma/BCG_E_Corrected_DeltaSigma.txt", np.c_[DeltaSigma[0], Corr_LR, Corr_LG, Corr_LI, Corr_SR, Corr_SG, Corr_SI, Error_LR, Error_LG, Error_LI, Error_SR, Error_SG, Error_SI])

        # Save Stages of Correction (Boost and Photo-Z onto a file)
        np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Ellipticity/Calculated_Delta_Sigma/Corrected_DeltaSigma/BCG_E_StageCorrection_DeltaSigma.txt", np.c_[Corr_BST_LR, Corr_BST_LG, Corr_BST_LI, Corr_BST_SR, Corr_BST_SG, Corr_BST_SI, Corr_BSTPZ_LR, Corr_BSTPZ_LG, Corr_BSTPZ_LI, Corr_BSTPZ_SR, Corr_BSTPZ_SG, Corr_BSTPZ_SI])

        # Save the Boost and the Random Signal to a file
        np.savetxt("/Users/emolina/PycharmProjects/JPL/Master_Function/Output_Files/Ellipticity/Calculated_Delta_Sigma/Corrected_DeltaSigma/BCG_E_Boost_Random.txt", np.c_[BoostCorr_LR, BoostCorr_LG, BoostCorr_LI, BoostCorr_SR, BoostCorr_SG, BoostCorr_SI, avg_PSF_L, avg_PSF_S])

        print "Done Correcting Ellipticity!"