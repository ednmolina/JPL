"""
PLOTS THE UNCORRECTED SIGNALS FOR THE KEYNOTE
"""
import numpy as np
import pylab
from astropy.io import fits as pyfits #Import astro py as this
import cosmology
import matplotlib.pyplot as plt
pylab.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size':12})
pylab.rc('text', usetex=True)

"""
IMPROT THE DATA
"""
# loads data from deltaSigma
f = np.loadtxt("/Users/emolina/PycharmProjects/JPL/deltaSigma.dat").transpose()
#X-Axis
x = f[0]


#Gets the Delta Sigmas
d = np.loadtxt("/Users/emolina/PycharmProjects/JPL/BCG_E_DeltaSigma.txt").transpose()
e = np.loadtxt("/Users/emolina/PycharmProjects/JPL/BCG_E_DeltaSigma_Error.txt").transpose()

DeltaSigma_L_R = d[1]
DeltaSigma_L_G = d[2]
DeltaSigma_L_I = d[3]

DeltaSigma_S_R = d[4]
DeltaSigma_S_G = d[5]
DeltaSigma_S_I = d[6]

error_L_r = e[0]
error_L_g = e[1]
error_L_i = e[2]

error_S_r = e[3]
error_S_g = e[4]
error_S_i = e[5]

"""
#For the uncorrected Signal
UncorrectedDSEM = np.average(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutTest/DeltaSigmaFiles/Uncorrected_Calculated_DeltaSigma/UncorrectedEMJKSS.txt"), axis = 0)
UncorrectedDSBM = np.average(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutTest/DeltaSigmaFiles/Uncorrected_Calculated_DeltaSigma/UncorrectedBMJKSS.txt"), axis = 0)

UncorrectedDSEM_Error = np.sqrt(np.diag(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutTest/DeltaSigmaFiles/Uncorrected_Calculated_DeltaSigma/UncorrectedEMJKCov.txt")))
UncorrectedDSBM_Error = np.sqrt(np.diag(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutTest/DeltaSigmaFiles/Uncorrected_Calculated_DeltaSigma/UncorrectedBMJKCov.txt")))

plt.figure(1)
plt.errorbar(x, UncorrectedDSEM, yerr = UncorrectedDSEM_Error)
plt.semilogx()
plt.semilogy()
plt.title("E Mode")
plt.xlabel(r"$R [h^{-1} \mathrm{Mpc}]$")
plt.ylabel(r"$\Delta \Sigma [h \mathrm{M}_{\odot}/\mathrm{pc^{2}}]$")
plt.xlim([x[2],10**1.8])

plt.figure(2)
plt.errorbar(x, UncorrectedDSBM, yerr = UncorrectedDSBM_Error)
plt.semilogx()
plt.axhline(y=0, color = 'k')
plt.title("B Mode")
plt.xlabel(r"$R [h^{-1} \mathrm{Mpc}]$")
plt.ylabel(r"$R [h^{-1} \mathrm{Mpc}] \times \Delta \Sigma_\times [h \mathrm{M_{\odot}}/\mathrm{pc^{2}}]$")
plt.xlim([x[2],10**1.8])"""
plt.figure(1, figsize=(22, 7))
plt.title("$\Delta \Sigma \mathrm{of the Uncorrected Subsamples} $")

plt.subplot(1,3,1)
plt.title("$e_R$")
plt.errorbar(x , DeltaSigma_L_R , yerr = error_L_r , label = '$e_R$ Large')
plt.errorbar(x , DeltaSigma_S_R , yerr = error_S_r , label = '$e_R$ Small')
plt.semilogx()
plt.semilogy()
plt.xlim([x[2],10**1.8])
plt.legend(loc=3)

plt.subplot(1,3,2)
plt.title("$e_G$")
plt.errorbar(x , DeltaSigma_L_G , yerr = error_L_g , label = '$e_G$ Large')
plt.errorbar(x , DeltaSigma_S_G , yerr = error_S_g , label = '$e_G$ Small')
plt.semilogx()
plt.semilogy()
plt.xlim([x[2],10**1.8])
plt.legend(loc=3)

plt.subplot(1,3,3)
plt.title("$e_I$")
plt.errorbar(x , DeltaSigma_L_I , yerr = error_L_i , label = '$e_I$ Large')
plt.errorbar(x , DeltaSigma_S_I , yerr = error_S_i , label = '$e_I$ Small')
plt.semilogx()
plt.semilogy()
plt.xlim([x[2],10**1.8])
plt.legend(loc=3)
plt.savefig("/Users/emolina/Desktop/Keynote Plots/UncorrectedEM.png", dpi = 326)

"""
PLOT THE BOOST OVER RADIAL BIN R OF THE PARENT SAMPLE
"""
#Boost Corrected Signal
Boost_DSEM = np.average(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutTest/Boost/Boost.txt"), axis = 0)

plt.figure(2)

plt.plot(x,Boost_DSEM)
plt.title("Boost Factor over Radial Bins")
plt.xlim([x[2],10**1.8])
plt.semilogx()
plt.xlabel(r"$R [h^{-1} \mathrm{Mpc}]$")
plt.ylabel("Boost Factor")
plt.savefig("/Users/emolina/Desktop/Keynote Plots/Boost Factor.png", dpi = 326)
plt.show()
