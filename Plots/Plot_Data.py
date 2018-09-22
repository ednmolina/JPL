"""
PLOTS PARENT SIGNAL
"""
import numpy as np
import pylab
from astropy.io import fits as pyfits #Import astro py as this
import cosmology
import matplotlib.pyplot as plt
pylab.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size':10})
pylab.rc('text', usetex=True)

"""
IMPROT THE DATA
"""
# loads data from deltaSigma
f = np.loadtxt("/Users/emolina/PycharmProjects/JPL/deltaSigma.dat").transpose()
#X-Axis
x = f[0]

#For the uncorrected Signal
UncorrectedDSEM = np.average(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutTest/DeltaSigmaFiles/Uncorrected_Calculated_DeltaSigma/UncorrectedEMJKSS.txt"), axis = 0)
UncorrectedDSBM = np.average(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutTest/DeltaSigmaFiles/Uncorrected_Calculated_DeltaSigma/UncorrectedBMJKSS.txt"), axis = 0)

UncorrectedDSEM_Error = np.sqrt(np.diag(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutTest/DeltaSigmaFiles/Uncorrected_Calculated_DeltaSigma/UncorrectedEMJKCov.txt")))
UncorrectedDSBM_Error = np.sqrt(np.diag(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutTest/DeltaSigmaFiles/Uncorrected_Calculated_DeltaSigma/UncorrectedBMJKCov.txt")))

#Boost Corrected Signal
Boost_DSEM = np.average(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages/EM_Boost.txt"), axis = 0)

Boost_DSEM_Error = np.sqrt(np.diagonal(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages/EM_BoostCov.txt")))

#Photo_Z Corrected Signal
Photo_Z_DSEM = np.average(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages/EM_Boost_PhotoZ.txt"), axis = 0)

Photo_Z_DSEM_Error = np.sqrt(np.diagonal(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages/EM_Boost_PhotoZCov.txt")))

#Corrected Signal
Corrected_DSEM = np.average(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutTest/CorrectedDeltaSigma/CorrectedEM.txt"), axis = 0)
Corrected_DSBM = np.average(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutTest/CorrectedDeltaSigma/CorrectedBM.txt"), axis = 0)

Corrected_DSEM_Erorr = np.sqrt(np.diagonal(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutTest/CovarianceAndCorrelation/Covariance_CorrectedDeltaSigmaEM.txt")))
Corrected_DSBM_Erorr = np.sqrt(np.diagonal(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutTest/CovarianceAndCorrelation/Covariance_CorrectedDeltaSigmaBM.txt")))

"""CALCULATE X**2 (CHI  SQUARED) FOR B MODE"""
chi = np.sum(Corrected_DSBM**2/Corrected_DSBM_Erorr**2)
print chi
"""
PLOT THE DATA
"""
plt.figure(1, figsize=(16, 20))

"""PLOTS THE UNORRECTED SIGNAL"""
plt.subplot(2, 1, 1)
plt.errorbar(x, UncorrectedDSEM, yerr = UncorrectedDSEM_Error, label = "Uncorrected Signal", linestyle = ':', color = 'b')
plt.errorbar(x, Boost_DSEM, yerr = Boost_DSEM_Error, label = "Boost Correction", linestyle = '--', color = 'b')
plt.errorbar(x, Photo_Z_DSEM, yerr = Photo_Z_DSEM_Error, label = "Boost and Photo-Z Correction", linestyle = '-.', color = 'b')
plt.errorbar(x, Corrected_DSEM, yerr = Corrected_DSEM_Erorr, label = "Boost, Photo-Z, and Random Correction")
plt.title("E Mode")
plt.xlabel(r"$R [h^{-1} \mathrm{Mpc}]$")
plt.ylabel(r"$\Delta \Sigma [h \mathrm{M}_{\odot}/\mathrm{pc^{2}}]$")
plt.semilogx()
plt.semilogy()
plt.xlim([x[2],10**1.8])
plt.legend()

plt.subplot(2, 1, 2)
plt.errorbar(x, UncorrectedDSBM*x, yerr = UncorrectedDSBM_Error*x, label = "Uncorrected Signal", linestyle = ':', color = 'b')
plt.errorbar(x, Corrected_DSBM*x, yerr = Corrected_DSBM_Erorr*x, label = "Random Correction", color = 'b')
plt.title("B Mode")
plt.xlabel(r"$R [h^{-1} \mathrm{Mpc}]$")
plt.ylabel(r"$R [h^{-1} \mathrm{Mpc}] \times \Delta \Sigma_\times [h \mathrm{M_{\odot}}/\mathrm{pc^{2}}]$")
plt.semilogx()
plt.xlim([x[2],10**1.8])
plt.axhline(y=0, color = 'k')
plt.legend()
plt.annotate("$\chi^2$/DoF =  %s/24" %round(chi, 2), xy =(23, 6.25), xytext=(23, 6.25))

plt.savefig("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/StageCorrectionPlotParent.png", dpi = 520, bbox_inches='tight')
print "Plot Saved!"
