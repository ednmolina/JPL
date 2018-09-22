"""
PLOTS THE 6 SUBSAMPLES (WP) AS CALCULATED ON ZODIAC
"""

import numpy as np
import pylab
from astropy.io import fits as pyfits #Import astro py as this
import cosmology
import matplotlib.pyplot as plt
import os
pylab.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size':10})
pylab.rc('text', usetex=True)

"""
LOAD THE CORRECTED DATA
"""

L_R = np.loadtxt("/Users/edenmolina/PycharmProjects/JPL/Subsample/Subsamples_Calculated_Zodiac_Data/Large/Red/wp.dat").transpose() #Array 2, 20
L_R_Error = np.sqrt(np.diagonal(np.loadtxt("/Users/edenmolina/PycharmProjects/JPL/Subsample/Subsamples_Calculated_Zodiac_Data/Large/Red/wp.cov.dat")))

L_G = np.loadtxt("/Users/edenmolina/PycharmProjects/JPL/Subsample/Subsamples_Calculated_Zodiac_Data/Large/Green/wp.dat").transpose() #Array 2, 20
L_G_Error = np.sqrt(np.diagonal(np.loadtxt("/Users/edenmolina/PycharmProjects/JPL/Subsample/Subsamples_Calculated_Zodiac_Data/Large/Green/wp.cov.dat")))

L_I = np.loadtxt("/Users/edenmolina/PycharmProjects/JPL/Subsample/Subsamples_Calculated_Zodiac_Data/Large/I/wp.dat").transpose() #Array 2, 20
L_I_Error = np.sqrt(np.diagonal(np.loadtxt("/Users/edenmolina/PycharmProjects/JPL/Subsample/Subsamples_Calculated_Zodiac_Data/Large/I/wp.cov.dat")))

S_R = np.loadtxt("/Users/edenmolina/PycharmProjects/JPL/Subsample/Subsamples_Calculated_Zodiac_Data/Small/Red/wp.dat").transpose() #Array 2, 20
S_R_Error = np.sqrt(np.diagonal(np.loadtxt("/Users/edenmolina/PycharmProjects/JPL/Subsample/Subsamples_Calculated_Zodiac_Data/Small/Red/wp.cov.dat")))

S_G = np.loadtxt("/Users/edenmolina/PycharmProjects/JPL/Subsample/Subsamples_Calculated_Zodiac_Data/Small/Green/wp.dat").transpose() #Array 2, 20
S_G_Error = np.sqrt(np.diagonal(np.loadtxt("/Users/edenmolina/PycharmProjects/JPL/Subsample/Subsamples_Calculated_Zodiac_Data/Small/Green/wp.cov.dat")))

S_I = np.loadtxt("/Users/edenmolina/PycharmProjects/JPL/Subsample/Subsamples_Calculated_Zodiac_Data/Small/I/wp.dat").transpose() #Array 2, 20
S_I_Error = np.sqrt(np.diagonal(np.loadtxt("/Users/edenmolina/PycharmProjects/JPL/Subsample/Subsamples_Calculated_Zodiac_Data/Small/I/wp.cov.dat")))

print L_R[0]
"""
PLOT THE SUBSAMPLES
"""

plt.figure(1, figsize=(22.5, 10))

plt.subplot(2, 3, 1)
plt.title("Red Subsample")
plt.errorbar(L_R[0], L_R[1], yerr=L_R_Error, label = "Large Subsample")
plt.xlabel("$R [h^{-1} Mpc$]")
plt.ylabel("$w_p$")
plt.semilogx()
plt.yscale("log", nonposy='clip')
plt.xlim([10**.25, 10**2])
plt.legend(loc = 3)

plt.subplot(2, 3, 2)
plt.title("Green Subsample")
plt.errorbar(L_G[0], L_G[1], yerr=L_G_Error, label = "Large Subsample")
plt.xlabel("$R [h^{-1} Mpc$]")
plt.ylabel("$w_p$")
plt.semilogx()
plt.yscale("log", nonposy='clip')
plt.xlim([10**.25, 10**2])
plt.legend(loc = 3)

plt.subplot(2, 3, 3)
plt.title("I Subsample")
plt.errorbar(L_I[0], L_I[1], yerr=L_I_Error, label = "Large Subsample")
plt.xlabel("$R [h^{-1} Mpc$]")
plt.ylabel("$w_p$")
plt.semilogx()
plt.yscale("log", nonposy='clip')
plt.xlim([10**.25, 10**2])
plt.legend(loc = 3)

plt.subplot(2, 3, 4)
plt.title("Red Subsample")
plt.errorbar(S_R[0], S_R[1], yerr=S_R_Error, label = "Small Subsample")
plt.xlabel("$R [h^{-1} Mpc$]")
plt.ylabel("$w_p$")
plt.semilogx()
plt.yscale("log", nonposy='clip')
plt.xlim([10**.25, 10**2])
plt.legend(loc = 3)

plt.subplot(2, 3, 5)
plt.title("Green Subsample")
plt.errorbar(S_G[0], S_G[1], yerr=S_G_Error, label = "Corrected Small Subsample")
plt.xlabel("$R [h^{-1} Mpc$]")
plt.ylabel("$w_p$")
plt.semilogx()
plt.yscale("log", nonposy='clip')
plt.xlim([10**.25, 10**2])
plt.legend(loc = 3)

plt.subplot(2, 3, 6)
plt.title("I Subsample")
plt.errorbar(S_I[0], S_I[1], yerr=S_I_Error, label = "Corrected Small Subsample")
plt.xlabel("$R [h^{-1} Mpc$]")
plt.ylabel("$w_p$")
plt.semilogx()
plt.yscale("log", nonposy='clip')
plt.xlim([10**.25, 10**2])
plt.legend(loc = 3)


plt.savefig("/Users/edenmolina/PycharmProjects/JPL/Subsample/Plots/SixSubSamplePlot_From_Zodiac.png", dpi = 520, bbox_inches='tight')
print "Figure 1!"


plt.figure(2, figsize=(22.5, 5))

plt.subplot(1, 3, 1)
plt.title("$e_G$")
plt.errorbar(L_G[0], L_G[1], yerr=L_G_Error, label="Large $e_G$")
plt.errorbar(S_G[0], S_G[1], yerr=S_G_Error, label="Small $e_G$")
plt.xlabel(r"$R [h^{-1} \mathrm{Mpc}]$")
plt.ylabel("$w_p$")
plt.semilogx()
plt.xlim([10**.25, 10**2])
plt.ylim([10 ** 0, 10 ** 3])
plt.yscale("log", nonposy='clip')
plt.legend(loc=3)


plt.subplot(1, 3, 2)
plt.title("$e_R$")
plt.errorbar(L_R[0], L_R[1], yerr=L_R_Error, label = "Large $e_R$")
plt.errorbar(S_R[0], S_R[1], yerr=S_R_Error, label = "Small $e_R$")
plt.xlabel(r"$R [h^{-1} \mathrm{Mpc}]$")
plt.ylabel("$w_p$")
plt.semilogx()
plt.xlim([10**.25, 10**2])
plt.ylim([10 ** 0, 10 ** 3])
plt.yscale("log", nonposy='clip')
plt.legend(loc = 3)

plt.subplot(1, 3, 3)
plt.title("$e_I$")
plt.errorbar(L_I[0], L_I[1], yerr=L_I_Error, label = "Large $e_I$")
plt.errorbar(S_I[0], S_I[1], yerr=S_I_Error, label = "Small $e_I$")
plt.xlabel(r"$R [h^{-1} \mathrm{Mpc}]$")
plt.ylabel("$w_p$")
plt.semilogx()
plt.xlim([10**.25, 10**2])
plt.ylim([10 ** 0, 10 ** 3])
plt.yscale("log", nonposy='clip')
plt.legend(loc = 3)

plt.savefig("/Users/edenmolina/PycharmProjects/JPL/Subsample/Plots/SixSubSamplePlot_From_Zodiac_Combined.png", dpi = 520, bbox_inches='tight')
print "Figure 2!"
plt.show()