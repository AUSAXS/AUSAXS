import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("temp/debye_waller_fitting.txt", skiprows=1)
exv = data[:, 0]
atomic = data[:, 1]
exv_single = data[:, 2]

plt.figure()
plt.hist(exv, bins=25, label='exv factors')
plt.show()

plt.figure()
plt.hist(atomic, bins=25, label='atomic factors')
plt.show()

plt.figure()
plt.hist(exv_single, bins=25, label='exv factors (no atomic)')
plt.show()