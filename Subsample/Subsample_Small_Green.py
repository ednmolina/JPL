"""
MAKE FITS FILE FOR THE SMALL GREEN SUBSAMPLE
"""
import numpy as np
from astropy.io import fits as pyfits

RedMaPPerPath = "/Users/emolina/PycharmProjects/JPL/redmapper_dr8_public_v5.10_catalog.value_added.Z_LAMBDA_0.1_0.33.LAMBDA_lt100-2.fits"
SubsampleIndicesPath = "/Users/emolina/PycharmProjects/JPL/BCG_E_SIndices.txt"
SubsampleColumn = 1
# define output directory name
OutputFilename = "redmapper_dr8_public_v5.10_catalog.value_added.Z_LAMBDA_0.1_0.33.LAMBDA_lt100-2" + ".subsample_G_S" ".fits"

"""
SET THE FILES
"""
# load subsample indices
SubsampleIndices = np.loadtxt(SubsampleIndicesPath, dtype=np.int).transpose()[SubsampleColumn]
# load lens data
RedMaPPer = pyfits.getdata(RedMaPPerPath)
pyfits.writeto(OutputFilename, RedMaPPer[SubsampleIndices])