import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

def pepsi_expansion_checker():
    dr, A, q = sp.symbols('dr A q')
    Vi = A * (1 - dr)**3
    g = Vi * sp.exp(-sp.pi * q**2 * sp.cbrt(Vi)**2)
    g2 = g.series(dr, x0=0, n=2).removeO()
    g3 = g.series(dr, x0=0, n=3).removeO()

    vals = {"A": 1, "dr": 0.1}
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

def pepsi_comparison():
    rm = 1.62
    r0, rw, q = sp.symbols('r0 rw q')
    dr = r0 - rm
    V_crysol = 4*np.pi/3 * rw**3
    V_pepsi = 4*np.pi/3 * rw**3 * (1 - dr)**3

    Gq = (r0/rm)**3 * sp.exp(-sp.sqrt(4*sp.pi/3)**3 * sp.pi * q*q * (r0**2 - rm**2))
    g_crysol = Gq*V_crysol*sp.exp(-sp.pi * q**2 * sp.cbrt(V_crysol)**2)
    # g_pepsi_12 = V_crysol*sp.exp(-sp.pi*q**2*sp.cbrt(V_crysol)**2)*(1 + dr*(3 - sp.cbrt(4*sp.pi/3)**2 * 2*q*q*sp.pi*rm*rm))
    g_pepsi_12 = V_crysol*sp.exp(-sp.pi*q**2*sp.cbrt(V_crysol)**2)*(1 + dr*(3 - sp.cbrt(4*sp.pi/3)**2 * 2*q*q*sp.pi*rm*rm))

    g_pepsi_11 = V_crysol*sp.exp(-sp.pi*q**2*sp.cbrt(V_crysol)**2)*(1 + dr*(3 - 2*sp.pi*q*q*sp.cbrt(V_crysol)**2))
    g_pepsi_approx = V_pepsi*sp.exp(-sp.pi*q**2*sp.cbrt(V_pepsi)**2)
    g_pepsi_approx_2 = g_pepsi_approx.series(dr, x0=0, n=2).removeO()
    # g_pepsi_approx_3 = g_pepsi_approx.series(dr, x0=0, n=3).removeO()

    vals = {"rw": 1.58, "r0": 0.95*rm}
    q_axis = np.linspace(0, 1, 100)
    g_crysol_y   = [  g_crysol.subs({rw: vals["rw"], r0: vals["r0"], q: q_val}) for q_val in q_axis]
    g_pepsi_12_y = [g_pepsi_12.subs({rw: vals["rw"], r0: vals["r0"], q: q_val}) for q_val in q_axis]
    g_pepsi_11_y = [g_pepsi_11.subs({rw: vals["rw"], r0: vals["r0"], q: q_val}) for q_val in q_axis]
    g_pepsi_approx_y   = [  g_pepsi_approx.subs({rw: vals["rw"], r0: vals["r0"], q: q_val}) for q_val in q_axis]
    g_pepsi_approx_2_y = [g_pepsi_approx_2.subs({rw: vals["rw"], r0: vals["r0"], q: q_val}) for q_val in q_axis]
    # g_pepsi_approx_3_y = [g_pepsi_approx_3.subs({rw: vals["rw"], r0: vals["r0"], q: q_val}) for q_val in q_axis]

    fig, ax = plt.subplots(2, 1, height_ratios=[2, 1], figsize=(10, 6))
    plt.sca(ax[0])
    plt.plot(q_axis, g_crysol_y, label='crysol')
    plt.plot(q_axis, g_pepsi_11_y, label='pepsi-11')
    plt.plot(q_axis, g_pepsi_12_y, label='pepsi-12')
    plt.plot(q_axis, g_pepsi_approx_y,   "--", lw=2, label='pepsi approx',   color='g')
    plt.plot(q_axis, g_pepsi_approx_2_y, ":",  lw=2, label='pepsi approx 2', color='b')
    # plt.plot(q_axis, g_pepsi_approx_3_y, "-.", lw=2, label='pepsi approx 3', color='y')
    plt.xlabel('q')
    plt.ylabel('Value')
    plt.title('Comparing crysol and pepsi functions')
    plt.legend()
    plt.grid(True)

    plt.sca(ax[1])
    plt.axhline(0, color='k', linestyle='--')
    plt.plot(q_axis, [100*(g_pepsi_12_y[i] - g_crysol_y[i]) / g_crysol_y[i] for i in range(len(q_axis))], "r")
    plt.plot(q_axis, [100*(g_pepsi_11_y[i] - g_crysol_y[i]) / g_crysol_y[i] for i in range(len(q_axis))], "b")
    plt.plot(q_axis, [100*(g_pepsi_approx_y[i] - g_crysol_y[i]) / g_crysol_y[i] for i in range(len(q_axis))],   "--", lw=2, color="g")
    # plt.plot(q_axis, [100*(g_pepsi_approx_2_y[i] - g_crysol_y[i]) / g_crysol_y[i] for i in range(len(q_axis))], ":" , lw=2, color="b")
    # plt.plot(q_axis, [100*(g_pepsi_approx_3_y[i] - g_crysol_y[i]) / g_crysol_y[i] for i in range(len(q_axis))], "-.", lw=2, color="y")
    plt.xlabel('q')
    plt.ylabel('Deviation [%]')
    plt.show()

pepsi_comparison()