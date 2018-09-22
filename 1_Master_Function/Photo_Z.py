"""
CALCULATE PHOTO-Z
"""
import numpy as np
from astropy.io import fits as pyfits #Import astro py as this
import matplotlib.pyplot as plt
import os
import pylab
pylab.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size':10})
pylab.rc('text', usetex=True)



def setFitsFiles(LensJKIndices, Combined_Everything, PhotoZDir):
    for i in range(83):
        index = LensJKIndices
        C = Combined_Everything

        """LOAD THE DATA FROM THE C FILE (FITS FILE CONTAINING LAMBA AND Z LAMBDA)"""
        Lambda = C["LAMBDA"]
        Z_Lambda = C["Z_LAMBDA"]

        Need_to_Delete = np.where(index == i)[0]  # Contains the indicies of the data points needed to be deleted
        New_Lambda = np.delete(Lambda, Need_to_Delete, 0)  # Indicies that need to be corresponded with the sums.fits file
        New_Z_Lambda = np.delete(Z_Lambda, Need_to_Delete, 0)

        new_cols = pyfits.ColDefs([pyfits.Column(name='LAMBDA', format='D', array=New_Lambda),
                                   pyfits.Column(name='Z_LAMBDA', format='D', array=New_Z_Lambda)])
        hdu = pyfits.BinTableHDU.from_columns(new_cols)

        hdu.writeto(PhotoZDir + "/Photo_Z%s.fits" %i)