import numpy as np
from astropy.io import fits as pyfits #Import astro py as this
import matplotlib.pyplot as plt
import pylab
pylab.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size':10})
pylab.rc('text', usetex=True)

"""
LOAD THE CALCULATED DATA
"""
# loads data from deltaSigma
f = np.loadtxt("/Users/edenmolina/PycharmProjects/JPL/deltaSigma.dat").transpose()

#X-Axis
x = f[0]

CorrectedEM = np.average(np.loadtxt("/Users/edenmolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutTest/CorrectedDeltaSigma/CorrectedEM.txt"), axis = 0)
CorrectedBM = np.average(np.loadtxt("/Users/edenmolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutTest/CorrectedDeltaSigma/CorrectedBM.txt"), axis = 0)


"""LOAD THE ERRORS"""
Error_EM = np.sqrt(np.diagonal(np.loadtxt("/Users/edenmolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutTest/CovarianceAndCorrelation/Covariance_CorrectedDeltaSigmaEM.txt")))
Error_BM = np.sqrt(np.diagonal(np.loadtxt("/Users/edenmolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutTest/CovarianceAndCorrelation/Covariance_CorrectedDeltaSigmaBM.txt")))

Correlation_EM = np.loadtxt("/Users/edenmolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutTest/CovarianceAndCorrelation/Correlation_CorrectedDeltaSigmaEM.txt")
Correlation_BM = np.loadtxt("/Users/edenmolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutTest/CovarianceAndCorrelation/Correlation_CorrectedDeltaSigmaBM.txt")

Correlation_Org = np.loadtxt("/Users/edenmolina/PycharmProjects/JPL/BCG_E_JackKnifeCorrection_ParentSample_V2.0/Files/Correlation.txt")

JackKNifeER = np.sqrt(np.sum(np.loadtxt("/Users/edenmolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutTest/JackKnife/CorrectedDSEM_JackKnife.txt"), axis = 0))
JackKnifeBR = np.sqrt(np.sum(np.loadtxt("/Users/edenmolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutTest/JackKnife/CorrectedDSBM_JackKnife.txt"), axis = 0))
"""
PLOT THE DATA
"""
plt.figure(1)
plt.errorbar(x, CorrectedEM, yerr = Correlation_EM)
plt.title("E-Mode")
plt.xlabel("$R [h^{-1} Mpc$]")
plt.ylabel("$R [h^{-1} Mpc$] x $\Delta$ $\Sigma$ [$h M_{\odot}/pc^{2}$]")
plt.semilogx()
plt.semilogy()

plt.figure(2)
plt.errorbar(x, CorrectedBM*x, yerr = Error_BM*x)
plt.title("B-Mode")
plt.xlabel("$R [h^{-1} Mpc$]")
plt.ylabel("$R [h^{-1} Mpc$] x $\Delta$ $\Sigma$ [$h M_{\odot}/pc^{2}$]")
plt.semilogx()

plt.figure(3)
plt.imshow(Correlation_EM, origin = 'lower left', interpolation = 'none')
plt.title("Correlation of Corrected E-Mode")
plt.colorbar()

plt.figure(4)
plt.imshow(Correlation_BM, origin = 'lower left', interpolation = 'none')
plt.title("Correlation of Corrected B-Mode")
plt.colorbar()

plt.show()