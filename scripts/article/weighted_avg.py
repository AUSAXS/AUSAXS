import numpy as np

import matplotlib.pyplot as plt

file_names = [
    'output/datcombine/1.dat', 
    'output/datcombine/2.dat'
]
x_values = []
y_values = []
yerr_values = []

for file_name in file_names:
    data = np.loadtxt(file_name, skiprows=1)
    x_values.append(data[:, 0])
    y_values.append(data[:, 1])
    yerr_values.append(data[:, 2])

x_values = np.array(x_values)
y_values = np.array(y_values)
yerr_values = np.array(yerr_values)
weights = 1 / yerr_values**2
weighted_avg = np.sum(weights * y_values, axis=0) / np.sum(weights, axis=0)
weighted_avg_err = 1 / np.sqrt(np.sum(weights, axis=0))

combined_data = np.loadtxt('output/datcombine/combined.dat', skiprows=2)

fig, ax = plt.subplots(2, 1, sharex=True, figsize=(16, 10))
plt.sca(ax[0])
plt.errorbar(x_values[0], yerr_values[0], label='Weighted average')
plt.errorbar(combined_data[:, 0], combined_data[:, 2], label='datcombine')
plt.legend()
plt.loglog()
plt.ylabel("calculated error")

weighted_avg_err /= yerr_values[0]
combined_data[:, 2] /= yerr_values[0]
plt.sca(ax[1])
plt.axhline(1/np.sqrt(2), linestyle='--', color='red', label='1/sqrt(N)')
plt.plot(x_values[0], weighted_avg_err, label='Weighted average')
plt.plot(combined_data[:, 0], combined_data[:, 2], label='datcombine')
plt.xlabel('q')
plt.ylabel('(calculated error) / (input error)')
plt.loglog()
plt.legend()
plt.show()
