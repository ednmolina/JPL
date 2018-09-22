import numpy as np
from astropy.io import fits as pyfits #Import astro py as this
import matplotlib.pyplot as plt
import cosmology
import pylab
pylab.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size':10})
pylab.rc('text', usetex=True)

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

def getPhotoZCorrection(dirname):
    Photo_Z = np.zeros((83))
    for i in range(83):
        #Calculate Photo-Z
        lens = pyfits.getdata(dirname + '/Photo_Z%s.fits' %i)
        bz = run(lens)
        correction_factor = 1. / (1. + bz)
        Photo_Z[i] = correction_factor
        print "Calculating PhotoZ: ", i
    return Photo_Z