import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

# definitions
dr, A, q = sp.symbols('dr A q')
Vi = A * (1 - dr)**3
g = Vi * sp.exp(-sp.pi * q**2 * sp.cbrt(Vi)**2)

rm = 1.62

# series expansion
g2 = g.series(dr, x0=0, n=2).removeO()
g3 = g.series(dr, x0=0, n=3).removeO()
g_pepsi = Vi*sp.exp(-sp.pi*q**2*sp.cbrt(Vi)**2)*(1 + dr*(3 - sp.cbrt(4*sp.pi/3)**2 * sp.pi*q*q*rm*rm))

# prepare the values for plotting
vals = {"A": 1, "dr": 0.1}
q_axis = np.linspace(0, 1, 100)
g_y =       [ g.subs({A: vals["A"], dr: vals["dr"], q: q_val}) for q_val in q_axis]
g2_y =      [g2.subs({A: vals["A"], dr: vals["dr"], q: q_val}) for q_val in q_axis]
g3_y =      [g3.subs({A: vals["A"], dr: vals["dr"], q: q_val}) for q_val in q_axis]
g_pepsi_y = [g_pepsi.subs({A: vals["A"], dr: vals["dr"], q: q_val}) for q_val in q_axis]

# Plot the original function g and its series expansion g3
fig, ax = plt.subplots(2, 1, height_ratios=[2, 1], figsize=(10, 6))

plt.sca(ax[0])
plt.plot(q_axis, g_y, "k", label='exact')
plt.plot(q_axis, g2_y, "--", label='first order')
plt.plot(q_axis, g3_y, "--", label='second order')
plt.plot(q_axis, g_pepsi_y, "--", label='pepsi')
plt.xlabel('q')
plt.ylabel('Value')
plt.title('Original function g and its series expansion')
plt.legend()
plt.grid(True)

# percentage deviation from the exact value
plt.sca(ax[1])
plt.axhline(0, color='k', linestyle='--')
plt.plot(q_axis, [100 * (g2_y[i] - g_y[i]) / g_y[i] for i in range(len(q_axis))])
plt.plot(q_axis, [100 * (g3_y[i] - g_y[i]) / g_y[i] for i in range(len(q_axis))])
plt.plot(q_axis, [100 * (g_pepsi_y[i] - g_y[i]) / g_y[i] for i in range(len(q_axis))])
plt.xlabel('q')
plt.ylabel('Deviation [%]')

plt.show()