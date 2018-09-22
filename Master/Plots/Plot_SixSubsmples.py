"""
PLOTS THE SUBSAMPLES CALCULATED FROM THE MASTER FILES
LARGE AND SMALL DELTA SIGMAS
"""
import numpy as np
from astropy.io import fits as pyfits #Import astro py as this
import cosmology
import matplotlib.pyplot as plt
import os
import pylab
pylab.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size':10})
pylab.rc('text', usetex=True)

# loads data from deltaSigma
f = np.loadtxt("/Users/emolina/PycharmProjects/JPL/deltaSigma.dat").transpose()
#X-Axis
x = f[0]


"""
LOAD THE UNCORRECTED DATA
"""
UC_L_R = np.average(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleRLarge/DeltaSigmaFiles/Uncorrected_Calculated_DeltaSigma/UncorrectedEMJKSS.txt"), axis = 0)
UC_L_R_Error = np.sqrt(np.diagonal(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleRLarge/DeltaSigmaFiles/Uncorrected_Calculated_DeltaSigma/UncorrectedEMJKCov.txt")))

UC_L_G = np.average(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleGLarge/DeltaSigmaFiles/Uncorrected_Calculated_DeltaSigma/UncorrectedEMJKSS.txt"), axis = 0)
UC_L_G_Error = np.sqrt(np.diagonal(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleGLarge/DeltaSigmaFiles/Uncorrected_Calculated_DeltaSigma/UncorrectedEMJKCov.txt")))

UC_L_I = np.average(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleILarge/DeltaSigmaFiles/Uncorrected_Calculated_DeltaSigma/UncorrectedEMJKSS.txt"), axis = 0)
UC_L_I_Error = np.sqrt(np.diagonal(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleILarge/DeltaSigmaFiles/Uncorrected_Calculated_DeltaSigma/UncorrectedEMJKCov.txt")))

UC_S_R = np.average(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleRSmall/DeltaSigmaFiles/Uncorrected_Calculated_DeltaSigma/UncorrectedEMJKSS.txt"), axis = 0)
UC_S_R_Error = np.sqrt(np.diagonal(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleRSmall/DeltaSigmaFiles/Uncorrected_Calculated_DeltaSigma/UncorrectedEMJKCov.txt")))

UC_S_G = np.average(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleGSmall/DeltaSigmaFiles/Uncorrected_Calculated_DeltaSigma/UncorrectedEMJKSS.txt"), axis = 0)
UC_S_G_Error = np.sqrt(np.diagonal(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleGSmall/DeltaSigmaFiles/Uncorrected_Calculated_DeltaSigma/UncorrectedEMJKCov.txt")))

UC_S_I = np.average(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleISmall/DeltaSigmaFiles/Uncorrected_Calculated_DeltaSigma/UncorrectedEMJKSS.txt"), axis = 0)
UC_S_I_Error = np.sqrt(np.diagonal(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleISmall/DeltaSigmaFiles/Uncorrected_Calculated_DeltaSigma/UncorrectedEMJKCov.txt")))

"""FOR B MODE"""
UC_L_R_BM = np.average(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleRLarge/DeltaSigmaFiles/Uncorrected_Calculated_DeltaSigma/UncorrectedBMJKSS.txt"), axis = 0)
UC_L_R_BM_Error = np.sqrt(np.diagonal(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleRLarge/DeltaSigmaFiles/Uncorrected_Calculated_DeltaSigma/UncorrectedBMJKCov.txt")))

UC_L_G_BM = np.average(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleGLarge/DeltaSigmaFiles/Uncorrected_Calculated_DeltaSigma/UncorrectedBMJKSS.txt"), axis = 0)
UC_L_G_BM_Error = np.sqrt(np.diagonal(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleGLarge/DeltaSigmaFiles/Uncorrected_Calculated_DeltaSigma/UncorrectedBMJKCov.txt")))

UC_L_I_BM = np.average(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleILarge/DeltaSigmaFiles/Uncorrected_Calculated_DeltaSigma/UncorrectedBMJKSS.txt"), axis = 0)
UC_L_I_BM_Error = np.sqrt(np.diagonal(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleILarge/DeltaSigmaFiles/Uncorrected_Calculated_DeltaSigma/UncorrectedBMJKCov.txt")))

UC_S_R_BM = np.average(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleRSmall/DeltaSigmaFiles/Uncorrected_Calculated_DeltaSigma/UncorrectedBMJKSS.txt"), axis = 0)
UC_S_R_BM_Error = np.sqrt(np.diagonal(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleRSmall/DeltaSigmaFiles/Uncorrected_Calculated_DeltaSigma/UncorrectedBMJKCov.txt")))

UC_S_G_BM = np.average(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleGSmall/DeltaSigmaFiles/Uncorrected_Calculated_DeltaSigma/UncorrectedBMJKSS.txt"), axis = 0)
UC_S_G_BM_Error = np.sqrt(np.diagonal(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleGSmall/DeltaSigmaFiles/Uncorrected_Calculated_DeltaSigma/UncorrectedBMJKCov.txt")))

UC_S_I_BM = np.average(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleISmall/DeltaSigmaFiles/Uncorrected_Calculated_DeltaSigma/UncorrectedBMJKSS.txt"), axis = 0)
UC_S_I_BM_Error = np.sqrt(np.diagonal(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleISmall/DeltaSigmaFiles/Uncorrected_Calculated_DeltaSigma/UncorrectedBMJKCov.txt")))


"""
LOAD THE STAGE CORRECTIONS
"""
#FOR BOOST
Boost_L_R = np.average(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Large/DSLR_Boost.txt"), axis = 0)
Boost_L_R_Error = np.sqrt(np.diagonal(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Large/DSLR_BoostCov.txt")))

Boost_L_G = np.average(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Large/DSLG_Boost.txt"), axis = 0)
Boost_L_G_Error = np.sqrt(np.diagonal(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Large/DSLG_BoostCov.txt")))

Boost_L_I = np.average(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Large/DSLI_Boost.txt"), axis = 0)
Boost_L_I_Error = np.sqrt(np.diagonal(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Large/DSLI_BoostCov.txt")))

Boost_S_R = np.average(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Small/DSSR_Boost.txt"), axis = 0)
Boost_S_R_Error = np.sqrt(np.diagonal(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Small/DSSR_BoostCov.txt")))

Boost_S_G = np.average(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Small/DSSG_Boost.txt"), axis = 0)
Boost_S_G_Error = np.sqrt(np.diagonal(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Small/DSSG_BoostCov.txt")))

Boost_S_I = np.average(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Small/DSSI_Boost.txt"), axis = 0)
Boost_S_I_Error = np.sqrt(np.diagonal(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Small/DSSI_BoostCov.txt")))

#FOR PHOTO-Z
PhotoZ_L_R = np.average(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Large/DSLR_PhotoZ.txt"), axis = 0)
PhotoZ_L_R_Error = np.sqrt(np.diagonal(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Large/DSLR_PhotoZCov.txt")))

PhotoZ_L_G = np.average(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Large/DSLG_PhotoZ.txt"), axis = 0)
PhotoZ_L_G_Error = np.sqrt(np.diagonal(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Large/DSLG_PhotoZCov.txt")))

PhotoZ_L_I = np.average(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Large/DSLI_PhotoZ.txt"), axis = 0)
PhotoZ_L_I_Error = np.sqrt(np.diagonal(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Large/DSLI_PhotoZCov.txt")))

PhotoZ_S_R = np.average(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Small/DSSR_PhotoZ.txt"), axis = 0)
PhotoZ_S_R_Error = np.sqrt(np.diagonal(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Small/DSSR_PhotoZCov.txt")))

PhotoZ_S_G = np.average(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Small/DSSG_PhotoZ.txt"), axis = 0)
PhotoZ_S_G_Error = np.sqrt(np.diagonal(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Small/DSSG_PhotoZCov.txt")))

PhotoZ_S_I = np.average(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Small/DSSI_PhotoZ.txt"), axis = 0)
PhotoZ_S_I_Error = np.sqrt(np.diagonal(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/Correction_Stages_Subsamples/Small/DSSI_PhotoZCov.txt")))

"""
LOAD THE CORRECTED DATA
"""
L_R = np.average(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleRLarge/CorrectedDeltaSigma/CorrectedEM.txt"), axis = 0)
L_R_Error = np.sqrt(np.diagonal(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleRLarge/CovarianceAndCorrelation/Covariance_CorrectedDeltaSigmaEM.txt")))

L_G = np.average(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleGLarge/CorrectedDeltaSigma/CorrectedEM.txt"), axis = 0)
L_G_Error = np.sqrt(np.diagonal(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleGLarge/CovarianceAndCorrelation/Covariance_CorrectedDeltaSigmaEM.txt")))

L_I = np.average(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleILarge/CorrectedDeltaSigma/CorrectedEM.txt"), axis = 0)
L_I_Error = np.sqrt(np.diagonal(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleILarge/CovarianceAndCorrelation/Covariance_CorrectedDeltaSigmaEM.txt")))

S_R = np.average(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleRSmall/CorrectedDeltaSigma/CorrectedEM.txt"), axis = 0)
S_R_Error = np.sqrt(np.diagonal(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleRSmall/CovarianceAndCorrelation/Covariance_CorrectedDeltaSigmaEM.txt")))

S_G = np.average(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleGSmall/CorrectedDeltaSigma/CorrectedEM.txt"), axis = 0)
S_G_Error = np.sqrt(np.diagonal(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleGSmall/CovarianceAndCorrelation/Covariance_CorrectedDeltaSigmaEM.txt")))

S_I = np.average(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleISmall/CorrectedDeltaSigma/CorrectedEM.txt"), axis = 0)
S_I_Error = np.sqrt(np.diagonal(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleISmall/CovarianceAndCorrelation/Covariance_CorrectedDeltaSigmaEM.txt")))

"""FOR B MODE"""
L_R_BM = np.average(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleRLarge/CorrectedDeltaSigma/CorrectedBM.txt"), axis = 0)
L_R_BM_Error = np.sqrt(np.diagonal(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleRLarge/CovarianceAndCorrelation/Covariance_CorrectedDeltaSigmaBM.txt")))
chi_LRBM = np.sum(L_R_BM**2/L_R_BM_Error**2)
print chi_LRBM

L_G_BM = np.average(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleGLarge/CorrectedDeltaSigma/CorrectedBM.txt"), axis = 0)
L_G_BM_Error = np.sqrt(np.diagonal(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleGLarge/CovarianceAndCorrelation/Covariance_CorrectedDeltaSigmaBM.txt")))
chi_LGBM = np.sum(L_G_BM**2/L_G_BM_Error**2)
print chi_LGBM

L_I_BM = np.average(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleILarge/CorrectedDeltaSigma/CorrectedBM.txt"), axis = 0)
L_I_BM_Error = np.sqrt(np.diagonal(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleILarge/CovarianceAndCorrelation/Covariance_CorrectedDeltaSigmaBM.txt")))
chi_LIBM = np.sum(L_I_BM**2/L_I_BM_Error**2)
print chi_LIBM

S_R_BM = np.average(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleRSmall/CorrectedDeltaSigma/CorrectedBM.txt"), axis = 0)
S_R_BM_Error = np.sqrt(np.diagonal(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleRSmall/CovarianceAndCorrelation/Covariance_CorrectedDeltaSigmaBM.txt")))
chi_SRBM = np.sum(S_R_BM**2/S_R_BM_Error**2)
print chi_SRBM

S_G_BM = np.average(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleGSmall/CorrectedDeltaSigma/CorrectedBM.txt"), axis = 0)
S_G_BM_Error = np.sqrt(np.diagonal(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleGSmall/CovarianceAndCorrelation/Covariance_CorrectedDeltaSigmaBM.txt")))
chi_SGBM = np.sum(S_G_BM**2/S_G_BM_Error**2)
print chi_SGBM

S_I_BM = np.average(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleISmall/CorrectedDeltaSigma/CorrectedBM.txt"), axis = 0)
S_I_BM_Error = np.sqrt(np.diagonal(np.loadtxt("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Output_Files/OutPutSubsampleISmall/CovarianceAndCorrelation/Covariance_CorrectedDeltaSigmaBM.txt")))
chi_SIBM = np.sum(S_I_BM**2/S_I_BM_Error**2)
print chi_SIBM

"""
PLOT THE SUBSAMPLES
"""
plt.figure(1, figsize=(22.5, 20))

plt.subplot(4, 3, 1)
plt.title("Large $e_G$")
plt.errorbar(x, UC_L_G, yerr=UC_L_G_Error, label="Uncorrected", color='b', linestyle=':')
plt.errorbar(x, Boost_L_G, yerr=Boost_L_G_Error, label="Boost Correction", color='b', linestyle='--')
plt.errorbar(x, PhotoZ_L_G, yerr=PhotoZ_L_G_Error, label="Boost and photo-z Correction", color='b', linestyle='-.')
plt.errorbar(x, L_G, yerr=L_G_Error, label="Boost, photo-z, and Random Correction", color='b')
plt.xlabel(r"$R [h^{-1} \mathrm{Mpc}]$")
plt.ylabel(r"$\Delta \Sigma [h \mathrm{M}_{\odot}/\mathrm{pc^{2}}]$")
plt.semilogx()
plt.yscale("log", nonposy='clip')
plt.xlim([x[2],10**1.8])
plt.legend(loc=3)


plt.subplot(4, 3, 2)
plt.title("Large $e_R$")
plt.errorbar(x, UC_L_R, yerr=UC_L_R_Error, label = "Uncorrected", color = 'b', linestyle = ':')
plt.errorbar(x, Boost_L_R, yerr=Boost_L_R_Error, label = "Boost Correction",color = 'b', linestyle = '--')
plt.errorbar(x, PhotoZ_L_R, yerr=PhotoZ_L_R_Error, label = "Boost and photo-z Correction", color = 'b', linestyle = '-.')
plt.errorbar(x, L_R, yerr=L_R_Error, label = "Boost, photo-z, and Random Correction", color = 'b')
plt.xlabel(r"$R [h^{-1} \mathrm{Mpc}]$")
plt.ylabel(r"$\Delta \Sigma [h \mathrm{M}_{\odot}/\mathrm{pc^{2}}]$")
plt.semilogx()
plt.yscale("log", nonposy='clip')
plt.xlim([x[2],10**1.8])
plt.legend(loc = 3)

plt.subplot(4, 3, 3)
plt.title("Large $e_I$")
plt.errorbar(x, UC_L_I, yerr=UC_L_I_Error, label = "Uncorrected", color = 'b', linestyle = ':')
plt.errorbar(x, Boost_L_I, yerr=Boost_L_I_Error, label = "Boost Correction", color = 'b', linestyle = '--')
plt.errorbar(x, PhotoZ_L_I, yerr=PhotoZ_L_I_Error, label = "Boost and photo-z Correction", color = 'b', linestyle = '-.')
plt.errorbar(x, L_I, yerr=L_I_Error, label = "Boost, photo-z, and Random Correction", color = 'b')
plt.xlabel(r"$R [h^{-1} \mathrm{Mpc}]$")
plt.ylabel(r"$\Delta \Sigma [h \mathrm{M}_{\odot}/\mathrm{pc^{2}}]$")
plt.semilogx()
plt.yscale("log", nonposy='clip')
plt.xlim([x[2],10**1.8])
plt.legend(loc = 3)


plt.subplot(4, 3, 4)
plt.title("Large $e_G$")
plt.errorbar(x, x * UC_L_G_BM, yerr=x * UC_L_G_BM_Error, label="Uncorrected", linestyle=':', color = 'b')
plt.errorbar(x, x * L_G_BM, yerr=x * L_G_BM_Error, label="Random Correction", color = 'b')
plt.xlabel(r"$R [h^{-1} \mathrm{Mpc}]$")
plt.ylabel(r"$R [h^{-1} \mathrm{Mpc}] \times \Delta \Sigma_\times [h \mathrm{M_{\odot}}/\mathrm{pc^{2}}]$")
plt.ylim([-10, 10])
plt.annotate("$\chi^2$ / DoF =  %s/24" % round(chi_LGBM, 2), xy=(0.1, 0.9), xycoords="axes fraction")
plt.semilogx()
plt.axhline(y=0, color='k')
plt.xlim([x[2],10**1.8])
plt.legend(loc=3)

plt.subplot(4, 3, 5)
plt.title("Large $e_R$")
plt.errorbar(x, x*UC_L_R_BM, yerr=x*UC_L_R_BM_Error, label = "Uncorrected",   linestyle = ':', color = 'b')
plt.errorbar(x, x*L_R_BM, yerr=x*L_R_BM_Error, label = "Random Correction", color = 'b')
plt.xlabel(r"$R [h^{-1} \mathrm{Mpc}]$")
plt.ylabel(r"$R [h^{-1} \mathrm{Mpc}] \times \Delta \Sigma_\times [h \mathrm{M}_{\odot}/\mathrm{pc^{2}}]$")
plt.ylim([-10, 10])
plt.annotate("$\chi^2$ / DoF =  %s/24" %round(chi_LRBM, 2), xy =(0.1, 0.9), xycoords="axes fraction")
plt.semilogx()
plt.axhline(y=0, color = 'k')
plt.xlim([x[2],10**1.8])
plt.legend(loc = 3)


plt.subplot(4, 3, 6)
plt.title("Large $e_I$")
plt.errorbar(x, x*UC_L_I_BM, yerr=x*UC_L_I_BM_Error, label = "Uncorrected",  linestyle = ':', color = 'b')
plt.errorbar(x, x*L_I_BM, yerr=x*L_I_BM_Error, label = "Random Correction", color = 'b')
plt.xlabel(r"$R [h^{-1} \mathrm{Mpc}]$")
plt.ylabel(r"$R [h^{-1} \mathrm{Mpc}] \times \Delta \Sigma_\times [h \mathrm{M}_{\odot}/\mathrm{pc^{2}}]$")
plt.ylim([-10, 10])
plt.annotate("$\chi^2$ / DoF =  %s/24" %round(chi_LIBM, 2), xy =(0.1, 0.9), xycoords="axes fraction")
plt.semilogx()
plt.axhline(y=0, color = 'k')
plt.xlim([x[2],10**1.8])
plt.legend(loc = 3)


plt.subplot(4, 3, 7)
plt.title("Small $e_G$")
plt.errorbar(x, UC_S_G, yerr=UC_S_G_Error, label="Uncorrected", color='b', linestyle=':')
plt.errorbar(x, Boost_S_G, yerr=Boost_S_G_Error, label="Boost Correction", color='b', linestyle='--')
plt.errorbar(x, PhotoZ_S_G, yerr=PhotoZ_S_G_Error, label="Boost and photo-z Correction", color='b', linestyle='-.')
plt.errorbar(x, S_G, yerr=S_G_Error, label="Boost, photo-z, and Random Correction", color='b')
plt.xlabel(r"$R [h^{-1} \mathrm{Mpc}]$")
plt.ylabel(r"$\Delta \Sigma [h \mathrm{M}_{\odot}/\mathrm{pc^{2}}]$")
plt.semilogx()
plt.yscale("log", nonposy='clip')
plt.xlim([x[2],10**1.8])
plt.legend(loc=3)


plt.subplot(4, 3, 8)
plt.title("Small $e_R$")
plt.errorbar(x, UC_S_R, yerr=UC_S_R_Error, label = "Uncorrected", color = 'b', linestyle = ':')
plt.errorbar(x, Boost_S_R, yerr=Boost_S_R_Error, label = "Boost Correction", color = 'b', linestyle = '--')
plt.errorbar(x, PhotoZ_S_R, yerr=PhotoZ_S_R_Error, label = "Boost and photo-z Correction", color = 'b', linestyle = '-.')
plt.errorbar(x, S_R, yerr=S_R_Error, label = "Boost, photo-z, and Random Correction", color = 'b')
plt.xlabel(r"$R [h^{-1} \mathrm{Mpc}]$")
plt.ylabel(r"$\Delta \Sigma [h \mathrm{M}_{\odot}/\mathrm{pc^{2}}]$")
plt.semilogx()
plt.yscale("log", nonposy='clip')
plt.xlim([x[2],10**1.8])
plt.legend(loc = 3)

plt.subplot(4, 3, 9)
plt.title("Small $e_I$")
plt.errorbar(x, UC_S_I, yerr=UC_S_I_Error, label = "Uncorrected", color = 'b', linestyle = ':')
plt.errorbar(x, Boost_S_I, yerr=Boost_S_I_Error, label = "Boost Correction", color = 'b', linestyle = '--')
plt.errorbar(x, PhotoZ_S_I, yerr=PhotoZ_S_I_Error, label = "Boost and photo-z Correction", color = 'b', linestyle = '-.')
plt.errorbar(x, S_I, yerr=S_I_Error, label = "Boost, photo-z, and Random Correction", color = 'b')
plt.xlabel(r"$R [h^{-1} \mathrm{Mpc}]$")
plt.ylabel(r"$\Delta \Sigma [h \mathrm{M}_{\odot}/\mathrm{pc^{2}}]$")
plt.semilogx()
plt.yscale("log", nonposy='clip')
plt.xlim([x[2],10**1.8])
plt.legend(loc = 3)


plt.subplot(4, 3, 10)
plt.title("Small $e_G$")
plt.errorbar(x, x * UC_S_G_BM, yerr=x * UC_S_G_BM_Error, label="Uncorrected", linestyle=':', color = 'b')
plt.errorbar(x, x * S_G_BM, yerr=x * S_G_BM_Error, label="Random Correction", color = 'b')
plt.xlabel(r"$R [h^{-1} \mathrm{Mpc}]$")
plt.ylabel(r"$R [h^{-1} \mathrm{Mpc}] \times \Delta \Sigma_\times [h \mathrm{M}_{\odot}/\mathrm{pc^{2}}]$")
plt.ylim([-10, 10])
plt.annotate("$\chi^2$ / DoF =  %s/24" % round(chi_SGBM, 2), xy=(0.1, 0.9), xycoords="axes fraction")
plt.semilogx()
plt.axhline(y=0, color='k')
plt.xlim([x[2],10**1.8])
plt.legend(loc=3)


plt.subplot(4, 3, 11)
plt.title("Small $e_R$")
plt.errorbar(x, x*UC_S_R_BM, yerr=x*UC_S_R_BM_Error, label = "Uncorrected",  linestyle = ':', color = 'b')
plt.errorbar(x, x*S_R_BM, yerr=x*S_R_BM_Error, label = "Random Correction", color = 'b')
plt.xlabel(r"$R [h^{-1} \mathrm{Mpc}]$")
plt.ylabel(r"$R [h^{-1} \mathrm{Mpc}] \times \Delta \Sigma_\times [h \mathrm{M}_{\odot}/\mathrm{pc^{2}}]$")
plt.ylim([-10, 10])
plt.annotate("$\chi^2$ / DoF =  %s/24" %round(chi_SRBM, 2), xy =(0.1, 0.9), xycoords="axes fraction")
plt.semilogx()
plt.axhline(y=0, color = 'k')
plt.xlim([x[2],10**1.8])
plt.legend(loc = 3)

plt.subplot(4, 3, 12)
plt.title("Small $e_I$")
plt.errorbar(x, x*UC_S_I_BM, yerr=x*UC_S_I_BM_Error, label = "Uncorrected",  linestyle = ':', color = 'b')
plt.errorbar(x, x*S_I_BM, yerr=x*S_I_BM_Error, label = "Random Correction", color = 'b')
plt.xlabel(r"$R [h^{-1} \mathrm{Mpc}]$")
plt.ylabel(r"$R [h^{-1} \mathrm{Mpc}] \times \Delta \Sigma_\times [h \mathrm{M}_{\odot}/\mathrm{pc^{2}}]$")
plt.ylim([-10, 10])
plt.annotate("$\chi^2$ / DoF =  %s/24" %round(chi_SIBM, 2), xy =(0.1, 0.9), xycoords="axes fraction")
plt.semilogx()
plt.axhline(y=0, color = 'k')
plt.xlim([x[2],10**1.8])
plt.legend(loc = 3)

plt.savefig("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/StageSixSubSamplePlot_BMEM.png", dpi = 520, bbox_inches='tight')
print "BM EM MODE Plot Saved!"

"""
PLOTTING JUST CORRECTED SIGNALS
"""
plt.figure(3, figsize=(22.5, 5))

plt.subplot(1, 3, 1)
plt.title("$e_G$")
plt.errorbar(x, L_G, yerr=L_G_Error, label = "Large $e_G$")
plt.errorbar(x, S_G, yerr=S_G_Error, label = "Small $e_G$")
plt.xlabel("$R [h^{-1} Mpc$]")
plt.ylabel(r"$\Delta \Sigma [h \mathrm{M}_{\odot}/\mathrm{pc^{2}}]$")
plt.semilogx()
plt.yscale("log", nonposy='clip')
plt.xlim([x[2],10**1.8])
plt.legend(loc = 3)

plt.subplot(1, 3, 2)
plt.title("$e_R$")
plt.errorbar(x, L_R, yerr=L_R_Error, label = "Large $e_R$")
plt.errorbar(x, S_R, yerr=S_R_Error, label = "Small $e_R$")
plt.xlabel("$R [h^{-1} Mpc$]")
plt.ylabel(r"$\Delta \Sigma [h \mathrm{M}_{\odot}/\mathrm{pc^{2}}]$")
plt.semilogx()
#plt.ylim([10**-100, 10**3])
plt.yscale("log", nonposy='clip')
plt.xlim([x[2],10**1.8])
plt.legend(loc = 3)

plt.subplot(1, 3, 3)
plt.title("$e_I$")
plt.errorbar(x, L_I, yerr=L_I_Error, label = "Large $e_I$")
plt.errorbar(x, S_I, yerr=S_I_Error, label = "Small $e_I$")
plt.xlabel(r"$R [h^{-1} \mathrm{Mpc}]$")
plt.ylabel(r"$\Delta \Sigma [h \mathrm{M}_{\odot}/\mathrm{pc^{2}}]$")
plt.semilogx()
plt.yscale("log", nonposy='clip')
plt.xlim([x[2],10**1.8])
plt.legend(loc = 3)

plt.savefig("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/SixSubSamplePlotCorrectedEM.png", dpi = 520, bbox_inches='tight')
print "E MODE Plot Saved!"
plt.show()

"""FOR B MODE"""
"""
plt.figure(4, figsize=(22.5, 10))

plt.subplot(1, 3, 1)
plt.title("$e_R$")
plt.errorbar(x, x*L_R_BM, yerr=x*L_R_BM_Error, label = "Corrected Large $e_R$")
plt.errorbar(x, x*S_R_BM, yerr=x*S_R_BM_Error, label = "Corrected Small $e_R$")
plt.xlabel("$R [h^{-1} Mpc$]")
plt.ylabel("$R [h^{-1} Mpc$] x $\Delta$ $\Sigma$ [$h M_{\odot}/pc^{2}$]")
plt.ylim([-10, 10])
plt.semilogx()
plt.axhline(y=0, color = 'k')
plt.legend(loc = 3)

plt.subplot(1, 3, 2)
plt.title("$e_G$")
plt.errorbar(x, x*L_G_BM, yerr=x*L_G_BM_Error, label = "Corrected Large $e_G$")
plt.errorbar(x, x*S_G_BM, yerr=x*S_G_BM_Error, label = "Corrected Small $e_G$")
plt.xlabel("$R [h^{-1} Mpc$]")
plt.ylabel("$R [h^{-1} Mpc$] x $\Delta$ $\Sigma$ [$h M_{\odot}/pc^{2}$]")
plt.ylim([-10, 10])
plt.semilogx()
plt.axhline(y=0, color = 'k')
plt.legend(loc = 3)

plt.subplot(1, 3, 3)
plt.title("$e_I$")
plt.errorbar(x, x*L_I_BM, yerr=x*L_I_BM_Error, label = "Corrected Large $e_I$")
plt.errorbar(x, x*S_I_BM, yerr=x*S_I_BM_Error, label = "Corrected Large $e_I$")
plt.xlabel("$R [h^{-1} Mpc$]")
plt.ylabel("$R [h^{-1} Mpc$] x $\Delta$ $\Sigma$ [$h M_{\odot}/pc^{2}$]")
plt.ylim([-10, 10])
plt.semilogx()
plt.axhline(y=0, color = 'k')
plt.legend(loc = 3)

plt.savefig("/Users/emolina/PycharmProjects/JPL/Master_FunctionCopy/Plots/SixSubSamplePlotCorrectedBM.png", dpi = 520, bbox_inches='tight')
print "B MODE Plot Saved!
"""