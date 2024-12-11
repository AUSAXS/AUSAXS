import matplotlib.pyplot as plt
import numpy as np
from enum import Enum

params = {
    'legend.fontsize': 20,
    'figure.figsize': (10, 8),
    'axes.labelsize': 20,
    'axes.titlesize': 20,
    'xtick.labelsize': 18,
    'ytick.labelsize': 18,
    'lines.markersize': 10,
    'lines.linewidth': 2,
}
plt.rcParams.update(params)

rm = 1.62
use_original_expressions = True
class FormFactor(Enum):
    CH = 1
    NH = 2

def get_traube_volume(ff):
    match (ff):
        case FormFactor.CH:
            return 21.59
        case FormFactor.NH:
            return 7.64
        case _:
            raise ValueError("Unknown form factor")

def get_schaefer_volume(ff):
    match (ff):
        # Voronoi
        # case FormFactor.CH:
        #     return 12.430
        # case FormFactor.NH:
        #     return 14.944

        # Minimum fluctuation
        case FormFactor.CH:
            return 11.640
        case FormFactor.NH:
            return 2.181
        case _:
            raise ValueError("Unknown form factor")
        
def get_radius(vol):
    return np.power(3*vol/(4*np.pi), 1/3)

def V_crysol(r):
    return 4*np.pi/3*r**3

def V_gauss(r):
    return np.power(np.pi, 3/2)*r**3

def ff_exact(q, r0, Va):
    V = Va*np.power(r0/rm, 3)
    return V*np.exp(-q*q*np.power(V, 2./3)/(4*np.pi))

def ff_crysol(q, r0, Va):
    # original
    if use_original_expressions:
        Gq = (r0/rm)**3 * np.exp(-np.power(4*np.pi/3, 3./2) * np.pi * q**2 * (r0**2 - rm**2))
        return Gq * Va * np.exp(-np.pi * q**2 * np.power(Va, 2./3)/(4*np.pi*np.pi))

    # likely implementation
    Gq = (r0/rm)**3 * np.exp(-np.power(4*np.pi/3, 2./3) * np.pi * q**2 * (r0**2 - rm**2) / (4*np.pi*np.pi))
    return Gq * Va * np.exp(-np.pi * q**2 * np.power(Va, 2./3)/(4*np.pi*np.pi))

def ff_pepsi(q, r0, Va):
    x = (r0 - rm)/rm

    # original
    if use_original_expressions:
        Gq = (1 + x*(3 - 2*np.pi*np.power(V_crysol(rm), 2./3)))
        return Gq * Va * np.exp(-np.pi*q**2*np.power(Va, 2./3)/(4*np.pi*np.pi))

    # likely implementation
    Gq = 1 + x*(3 - 2*np.pi*np.power(V_crysol(rm), 2./3)/(4*np.pi*np.pi))
    return Gq * Va * np.exp(-np.pi*q**2*np.power(Va, 2./3)/(4*np.pi*np.pi))

    # correct derivation
    # Gq = 1 - x*(3 - np.power(Va, 2./3)*q*q/(2*np.pi))
    # return Gq * Va * np.exp(-q*q*np.power(Va, 2./3)/(4*np.pi))

def ff_foxs(q, r0, ff):
    match (ff):
        case FormFactor.CH:
            scaling = 7.21106
        case FormFactor.NH:
            scaling = 2.55176
        case _:
            raise ValueError("Unknown form factor")

    # original from article
    if use_original_expressions:
        c1 = r0/rm
        Gq = c1**3 * np.exp(-np.power(4*np.pi/3, 3./2)*rm*rm*q*q*(c1*c1 - 1)/(4*np.pi))
        return scaling*Gq

    # original from source code
    c1 = r0/rm
    ff = scaling*np.exp(-0.23/2*q*q)
    coeff = -rm*rm*(c1*c1 - 1)/(4*np.pi)
    x = coeff*q*q
    Gq = c1**3 * np.exp(x)
    return ff*Gq

def plot_ff(ff: FormFactor):
    exact = {}
    crysol = {}
    pepsi = {}
    foxs = {}
    tag = ["min", "mid", "max"]
    q_axis = np.linspace(0, 1, 100)
    for l, t in zip([0.865, 1, 1.265], tag):
        crysol[t] = [ff_crysol(q, l*rm, get_traube_volume(ff)) for q in q_axis]

    for l, t in zip([0.92, 1, 1.08], tag):
        exact[t] = [ff_exact(q, l*rm, get_schaefer_volume(ff)) for q in q_axis]

    for l, t in zip([0.95, 1, 1.05], tag):
        pepsi[t] = [ff_pepsi(q, l*rm, get_traube_volume(ff)) for q in q_axis]
        foxs[t] =  [ ff_foxs(q, l*rm, ff) for q in q_axis]

    plt.plot(q_axis, exact["mid"], "k-", label='AUSAXS', lw=2)
    plt.plot(q_axis, exact["min"], "k--")
    plt.plot(q_axis, exact["max"], "k--")

    step = 10
    for i in range(0, len(crysol["mid"]), step):
        if i == 0:
            plt.plot(q_axis[i:i+step+1], crysol["mid"][i:i+step+1], "c-", label="CRYSOL", lw=2)
            continue
        elif i == 10:
            plt.plot(q_axis[i:i+step+1], pepsi["mid"][i:i+step+1], "b-", label='Pepsi-SAXS', lw=2)
            continue

        if i/step % 2 == 0:
            plt.plot(q_axis[i:i+step+1], crysol["mid"][i:i+step+1], "c-", lw=2)
        else:
            plt.plot(q_axis[i:i+step+1], pepsi["mid"][i:i+step+1], "b-", lw=2)
    plt.plot(q_axis, crysol["min"], "c--")
    plt.plot(q_axis, crysol["max"], "c--")

    plt.plot(q_axis, pepsi["min"], "b--")
    plt.plot(q_axis, pepsi["max"], "b--")

    plt.plot(q_axis, foxs["mid"], "-", color="tab:orange", label='FoXS', lw=2)
    plt.plot(q_axis, foxs["min"], "--", color="tab:orange")
    plt.plot(q_axis, foxs["max"], "--", color="tab:orange")

    plt.xlabel('$q\ (Å^{-1})$')
    plt.title(ff.name + ' excluded volume form factor')
    plt.grid(True)

### ORIGINAL PLOT ###
fig, ax = plt.subplots(1, 2, figsize=(16, 8), gridspec_kw={'wspace': 0.01})
plt.sca(ax[0])
plot_ff(FormFactor.CH)
plt.ylabel('Amplitude')
plt.ylim(-25, 70)

plt.sca(ax[1])
plt.gca().yaxis.tick_right()
plot_ff(FormFactor.NH)
plt.ylim(-10, 30)
plt.legend(loc='upper right', bbox_to_anchor=(0.5, 1.1), ncol=2)
plt.savefig("ff_ranges_original.png", dpi=300)
plt.show()

### CORRECTED PLOT ###
use_original_expressions = False
fig, ax = plt.subplots(1, 2, figsize=(16, 7), gridspec_kw={'wspace': 0.01})
plt.sca(ax[0])
plot_ff(FormFactor.CH)
plt.ylabel('Amplitude')

plt.sca(ax[1])
plt.gca().yaxis.tick_right()
plot_ff(FormFactor.NH)
plt.legend(loc='upper right', bbox_to_anchor=(0.5, 0.83), ncol=2)
plt.savefig("ff_ranges_corrected.png", dpi=300)
plt.show()