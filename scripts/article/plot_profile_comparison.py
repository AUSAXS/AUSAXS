import numpy as np
import matplotlib.pyplot as plt
import sys

params = {
    'legend.fontsize': 22,
    'figure.figsize': (10, 8),
    'axes.labelsize': 22,
    'axes.titlesize': 22,
    'xtick.labelsize': 18,
    'ytick.labelsize': 18,
    'lines.markersize': 10,
    'lines.linewidth': 3,
}
plt.rcParams.update(params)

path = sys.argv[1]
if path[-1] != '/':
    path += '/'
title = path.split('/')[-3]

crysol_aa =         np.loadtxt(path + 'crysol_aa.dat', skiprows=1)
foxs_aa =           np.loadtxt(path + 'foxs_aa.dat', skiprows=1)
pepsi_aa =          np.loadtxt(path + 'pepsi_aa.dat', skiprows=1)
ausaxs_aa =         np.loadtxt(path + 'ausaxs_aa.dat', skiprows=1)

crysol_xx =         np.loadtxt(path + 'crysol_xx.dat', skiprows=1)
foxs_xx =           np.loadtxt(path + 'foxs_xx.dat', skiprows=1)
pepsi_xx =          np.loadtxt(path + 'pepsi_xx.dat', skiprows=1)
ausaxs_xx =         np.loadtxt(path + 'ausaxs_xx.dat', skiprows=1)
ausaxs_crysol_xx =  np.loadtxt(path + 'ausaxs_crysol_xx.dat', skiprows=1)
ausaxs_foxs_xx =    np.loadtxt(path + 'ausaxs_foxs_xx.dat', skiprows=1)
ausaxs_pepsi_xx =   np.loadtxt(path + 'ausaxs_pepsi_xx.dat', skiprows=1)

crysol_aa_diff =    np.loadtxt(path + 'crysol_aa_diff.dat', skiprows=1)
foxs_aa_diff =      np.loadtxt(path + 'foxs_aa_diff.dat', skiprows=1)
pepsi_aa_diff =     np.loadtxt(path + 'pepsi_aa_diff.dat', skiprows=1)

crysol_xx_diff =    np.loadtxt(path + 'crysol_xx_diff.dat', skiprows=1)
foxs_xx_diff =      np.loadtxt(path + 'foxs_xx_diff.dat', skiprows=1)
pepsi_xx_diff =     np.loadtxt(path + 'pepsi_xx_diff.dat', skiprows=1)

crysol_aa_diff[:, 1] /= crysol_aa_diff[0, 1]
foxs_aa_diff[:, 1] /= foxs_aa_diff[0, 1]
pepsi_aa_diff[:, 1] /= pepsi_aa_diff[0, 1]

crysol_xx_diff[:, 1] /= crysol_xx_diff[0, 1]
foxs_xx_diff[:, 1] /= foxs_xx_diff[0, 1]
pepsi_xx_diff[:, 1] /= pepsi_xx_diff[0, 1]

fig, ax = plt.subplots(2, 1, figsize=(10, 8), sharex=True, gridspec_kw={'height_ratios': [3, 2], 'hspace': 0.1})
plt.sca(ax[0])
plt.plot(crysol_aa[:, 0], crysol_aa[:, 1], label='CRYSOL', color='tab:cyan')
plt.plot(foxs_aa[:, 0], foxs_aa[:, 1], label='FoXS', color='tab:orange')
plt.plot(pepsi_aa[:, 0], pepsi_aa[:, 1], label='Pepsi-SAXS', color='tab:blue')
plt.plot(ausaxs_aa[:, 0], ausaxs_aa[:, 1], label='AUSAXS', color='k')
plt.semilogy()
plt.ylabel('Intensity')
plt.legend()
plt.title(title + ' $I_{AA}$ profiles')
plt.sca(ax[1])
plt.axhline(1, color='k')
plt.plot(crysol_aa_diff[:, 0], crysol_aa_diff[:, 1], color='tab:cyan')
plt.plot(foxs_aa_diff[:, 0], foxs_aa_diff[:, 1], color='tab:orange')
plt.plot(pepsi_aa_diff[:, 0], pepsi_aa_diff[:, 1], color='tab:blue')
plt.ylim(0.95, 1.05)
plt.xlabel('$q$ [A⁻¹]')
plt.ylabel('mimic ratio')
plt.savefig(path + 'aa_comparison.png', dpi=600)

fig, ax = plt.subplots(2, 1, figsize=(10, 8), sharex=True, gridspec_kw={'height_ratios': [3, 2], 'hspace': 0.1})
plt.sca(ax[0])
plt.plot(crysol_xx[:, 0], crysol_xx[:, 1], label='CRYSOL', color='tab:cyan')
plt.plot(foxs_xx[:, 0], foxs_xx[:, 1], label='FoXS', color='tab:orange')
plt.plot(pepsi_xx[:, 0], pepsi_xx[:, 1], label='Pepsi-SAXS', color='tab:blue')
plt.plot(ausaxs_xx[:, 0], ausaxs_xx[:, 1], label='AUSAXS', color='k')
plt.plot(ausaxs_crysol_xx[:, 0], ausaxs_crysol_xx[:, 1], label='AUSAXS$_{mimics}$', color='k', lw=1, ls='--')
plt.plot(ausaxs_foxs_xx[:, 0], ausaxs_foxs_xx[:, 1], color='k', lw=1, ls='--')
plt.plot(ausaxs_pepsi_xx[:, 0], ausaxs_pepsi_xx[:, 1], color='k', lw=1, ls='--')
plt.semilogy()
plt.ylabel('Intensity')
plt.legend()
plt.title(title + ' $I_{XX}$ profiles')
plt.sca(ax[1])
plt.axhline(1, color='k')
plt.plot(crysol_xx_diff[:, 0], crysol_xx_diff[:, 1], color='tab:cyan')
plt.plot(foxs_xx_diff[:, 0], foxs_xx_diff[:, 1], color='tab:orange')
plt.plot(pepsi_xx_diff[:, 0], pepsi_xx_diff[:, 1], color='tab:blue')
plt.ylim(0.95, 1.05)
plt.xlabel('$q$ [A⁻¹]$')
plt.ylabel('mimic ratio')
plt.savefig(path + 'xx_comparison.png', dpi=600)