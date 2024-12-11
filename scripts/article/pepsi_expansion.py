import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

def pepsi_expansion_checker():
    dr, A, q = sp.symbols('dr A q')
    Vi = A * (1 - dr)**3
    g = Vi * sp.exp(-sp.pi * q**2 * sp.cbrt(Vi)**2)
    g2 = g.series(dr, x0=0, n=2).removeO()
    g3 = g.series(dr, x0=0, n=3).removeO()

    vals = {"A": 1, "dr": 0.05}
    q_axis = np.linspace(0, 1, 100)
    g_y =       [ g.subs({A: vals["A"], dr: vals["dr"], q: q_val}) for q_val in q_axis]
    g2_y =      [g2.subs({A: vals["A"], dr: vals["dr"], q: q_val}) for q_val in q_axis]
    g3_y =      [g3.subs({A: vals["A"], dr: vals["dr"], q: q_val}) for q_val in q_axis]

    fig, ax = plt.subplots(2, 1, height_ratios=[2, 1], figsize=(10, 6))
    plt.sca(ax[0])
    plt.plot(q_axis, g_y, "k", label='exact')
    plt.plot(q_axis, g2_y, "--", label='first order')
    plt.plot(q_axis, g3_y, "--", label='second order')
    plt.xlabel('q')
    plt.ylabel('Value')
    plt.title('Original function g and its series expansion')
    plt.legend()
    plt.grid(True)

    plt.sca(ax[1])
    plt.axhline(0, color='k', linestyle='--')
    plt.plot(q_axis, [100 * (g2_y[i] - g_y[i]) / g_y[i] for i in range(len(q_axis))])
    plt.plot(q_axis, [100 * (g3_y[i] - g_y[i]) / g_y[i] for i in range(len(q_axis))])
    plt.xlabel('q')
    plt.ylabel('Deviation [%]')

    plt.show()
#pepsi_expansion_checker()

def pepsi_comparison():
    rm = 1.62
    rw = 1.58
    x, A, q = sp.symbols('x A q')
    V = A
    Vdr = A * (1 - x)**3

    g_exact = Vdr*sp.exp(-sp.pi * q**2 * sp.cbrt(Vdr)**2)
    g_pepsi = V*sp.exp(-sp.pi * q**2 * sp.cbrt(V)**2) * (1 + x*(3 - sp.cbrt(4*sp.pi/3)**2 * 2*sp.pi*rm*rm))
    g_pepsi_corrected = V*sp.exp(-sp.pi * q**2 * sp.cbrt(V)**2) * (1 - x*(3 - 2*sp.pi*sp.cbrt(V)**2 * q*q))
    g_pepsi_approx_2 = g_exact.series(x, x0=0, n=2).removeO()

    vals = {"A": 4*sp.pi/3*rw**3, "x": 0.05}
    q_axis = np.linspace(0, 1, 100)
    g_exact_y = [g_exact.subs({A: vals["A"], x: vals["x"], q: q_val}) for q_val in q_axis]
    g_pepsi_y = [g_pepsi.subs({A: vals["A"], x: vals["x"], q: q_val}) for q_val in q_axis]
    g_pepsi_corrected_y = [g_pepsi_corrected.subs({A: vals["A"], x: vals["x"], q: q_val}) for q_val in q_axis]
    g_pepsi_approx_2_y = [g_pepsi_approx_2.subs({A: vals["A"], x: vals["x"], q: q_val}) for q_val in q_axis]

    fig, ax = plt.subplots(2, 1, height_ratios=[2, 1], figsize=(10, 6))
    plt.sca(ax[0])
    plt.plot(q_axis, g_exact_y, "k", label='exact')
    plt.plot(q_axis, g_pepsi_y, "r", label='pepsi')
    plt.plot(q_axis, g_pepsi_corrected_y, "--", lw=2, label='pepsi corrected', color='g')
    plt.plot(q_axis, g_pepsi_approx_2_y, ":",  lw=2, label='pepsi approx 2', color='b')
    plt.xlabel('q')
    plt.ylabel('Value')
    # plt.semilogy()
    plt.title('Comparing crysol and pepsi functions')
    plt.legend()
    plt.grid(True)

    plt.sca(ax[1])
    plt.axhline(1, color='k', linestyle='--')
    plt.plot(q_axis, [g_pepsi_y[i] / g_exact_y[i] for i in range(len(q_axis))], "r")
    plt.plot(q_axis, [g_pepsi_approx_2_y[i] / g_exact_y[i] for i in range(len(q_axis))], "b")
    plt.xlabel('q')
    plt.ylabel('Deviation factor')
    plt.show()
pepsi_comparison()

def fit_range():
    rm = 1.62
    rw_val = 1.92
    r0, rw, q = sp.symbols('r0 rw q')
    V = 4*np.pi/3 * rw**3
    V_exact = np.power(np.pi, 3./2)*rw**3
    dr = r0 - rm
    c1 = r0/rm
    Gq_crysol = (r0/rm)**3 * sp.exp(-sp.sqrt(4*sp.pi/3)**3 * sp.pi * q*q * (r0**2 - rm**2)/(4*sp.pi*sp.pi))
    Gq_exact  = (r0/rm)**3 * sp.exp(-(rw/rm)**2 * (r0**2 - rm**2)*q*q / 4)
    Gq_approx  = (r0/rm)**3 * sp.exp(-(r0**2 - rm**2)*q*q / 4)
    # Gq_foxs =        c1**3 * sp.exp(-sp.sqrt(4*sp.pi/3)**3*q*q*rm*rm*(c1*c1 - 1)/(4*sp.pi))
    Gq_foxs =        c1**3 * sp.exp(-q*q*rm*rm*(c1*c1 - 1)/(4*sp.pi))
    Gq_pepsi = (1 + dr*(3 - 2*sp.pi*sp.cbrt(4*sp.pi/3)**2 * rm*rm))

    g_exact = Gq_exact*V_exact*sp.exp(-q**2 * sp.cbrt(V_exact)**2/(4*sp.pi))
    g_approx = Gq_approx*V_exact*sp.exp(-q**2 * sp.cbrt(V_exact)**2/(4*sp.pi))
    g_crysol = Gq_crysol*V*sp.exp(-sp.pi * q**2 * sp.cbrt(V)**2/(4*sp.pi*sp.pi))
    g_pepsi = Gq_pepsi*V*sp.exp(-sp.pi*q**2*sp.cbrt(V)**2/(4*sp.pi*sp.pi))
    g_foxs = Gq_foxs*5.49096*sp.exp(-0.23/2*q*q)

    exact = {}
    approx = {}
    crysol = {}
    pepsi = {}
    foxs = {}
    tag = {"min", "max"}
    for l, t in zip([0.865, 1.265], tag):
        vals = {"rw": rw_val, "r0": l*rm}
        q_axis = np.linspace(0, 1, 100)
        crysol[t] = [g_crysol.subs({rw: vals["rw"], r0: vals["r0"], q: q_val}) for q_val in q_axis]

    for l, t in zip([0.95, 1.05], tag):
        vals = {"rw": rw_val, "r0": l*rm}
        q_axis = np.linspace(0, 1, 100)
        exact[t] = [g_exact.subs({rw: vals["rw"], r0: vals["r0"], q: q_val}) for q_val in q_axis]
        approx[t] = [g_approx.subs({rw: vals["rw"], r0: vals["r0"], q: q_val}) for q_val in q_axis]
        pepsi[t] = [g_pepsi.subs({rw: vals["rw"], r0: vals["r0"], q: q_val}) for q_val in q_axis]   
        foxs[t] = [g_foxs.subs({rw: vals["rw"], r0: vals["r0"], q: q_val}) for q_val in q_axis]

    plt.figure(figsize=(10, 6))
    plt.plot(q_axis, exact["min"], "k--", label='exact')
    plt.plot(q_axis, exact["max"], "k--")
    plt.plot(q_axis, approx["min"], "g--", label='approx')
    plt.plot(q_axis, approx["max"], "g--")
    plt.plot(q_axis, crysol["min"], "c", label='crysol')
    plt.plot(q_axis, crysol["max"], "c")
    plt.plot(q_axis, pepsi["min"], "b", label='pepsi')
    plt.plot(q_axis, pepsi["max"], "b")
    plt.plot(q_axis, foxs["min"], "tab:orange", label='foxs')
    plt.plot(q_axis, foxs["max"], "tab:orange")
    plt.xlabel('q')
    plt.ylabel('Value')
    plt.title('Range of fitting')
    plt.legend()
    plt.grid(True)
    plt.show()

    # plt.sca(ax[1])
    # plt.axhline(1, color='k', linestyle='--')
    # plt.plot(q_axis, [g_crysol_y[i] / g_exact_y[i] for i in range(len(q_axis))])
    # plt.plot(q_axis, [g_pepsi_y[i] / g_exact_y[i] for i in range(len(q_axis))])
    # plt.xlabel('q')
    # plt.ylabel('Deviation factor')
    # plt.show()
fit_range()