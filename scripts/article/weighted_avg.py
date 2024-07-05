import numpy as np
import subprocess

import matplotlib.pyplot as plt

# file = "data/SASDJG5/SASDJG5.dat"
file = "data/rigidbody/lysozyme/2epe.dat"
repeats = 20

data = np.loadtxt(file, skiprows=4)
with open(f'output/datcombine/1.dat', 'w') as f:
    for row in data:
        f.write(f"{row[0]:.8e}\t {row[1]:.8e}\t {row[2]:.8e}\n")

factors = {"wavg": [], "datcombine": []}
for i in range(1, repeats):
    x_values = []
    y_values = []
    yerr_values = []

    data = np.loadtxt(file, skiprows=4)
    x_values.append(data[:, 0])
    y_values.append(data[:, 1])
    yerr_values.append(data[:, 2])
    with open(f'output/datcombine/{i+1}.dat', 'w') as f:
        for row in data:
            f.write(f"{row[0]:.8e}\t {row[1] + i*1e-9:.8e}\t {row[2]:.8e}\n")

    for j in range(i):
        x_values.append(data[:, 0])
        y_values.append(data[:, 1] + j*1e-9)
        yerr_values.append(data[:, 2])

    x_values = np.array(x_values)
    y_values = np.array(y_values)
    yerr_values = np.array(yerr_values)
    weights = 1 / yerr_values**2
    weighted_avg = np.sum(weights * y_values, axis=0) / np.sum(weights, axis=0)
    weighted_avg_err = 1 / np.sqrt(np.sum(weights, axis=0))
    factors["wavg"].append(np.mean((weighted_avg_err / yerr_values[0])))

    input_files = [f'{j}.dat' for j in range(1, i+2)]
    subprocess.run([
        "datcombine", 
        "1.dat", 
        *input_files,
        "-o",
        "combined.dat",
        "--error-filter=N",
        "--outlier-filter=N"
    ], cwd='output/datcombine')
    combined_data = []
    for line in open('output/datcombine/combined.dat'):
        row = line.split()
        if len(row) != 3:
            continue
        try:
            row[0] = float(row[0])
            row[1] = float(row[1])
            row[2] = float(row[2])
            combined_data.append([row[0], row[1], row[2]])
        except ValueError:
            pass
    combined_data = np.array(combined_data)
    factors["datcombine"].append(np.mean((combined_data[:, 2] / yerr_values[0])))

    # fig, ax = plt.subplots(2, 1, sharex=True, figsize=(16, 10))
    # plt.sca(ax[0])
    # plt.errorbar(x_values[0], yerr_values[0], label='Weighted average')
    # plt.errorbar(combined_data[:, 0], combined_data[:, 2], label='datcombine')
    # plt.legend()
    # plt.loglog()
    # plt.ylabel("calculated error")

    # weighted_avg_err /= yerr_values[0]
    # combined_data[:, 2] /= yerr_values[0]
    # plt.sca(ax[1])
    # plt.axhline(1/np.sqrt(2), linestyle='--', color='red', label='1/sqrt(N)')
    # plt.plot(x_values[0], weighted_avg_err, label='Weighted average')
    # plt.plot(combined_data[:, 0], combined_data[:, 2], label='datcombine')
    # plt.xlabel('q')
    # plt.ylabel('(calculated error) / (input error)')
    # plt.loglog()
    # plt.legend()
    # plt.show()

print(factors)
N = np.arange(2, repeats+1)
plt.figure(figsize=(16, 10))
plt.plot(N, factors["wavg"], "*", label='Weighted average')
plt.plot(N, factors["datcombine"], "*", label='datcombine')
plt.plot(N, 1/np.sqrt(N), "k-", label='1/sqrt(N)')
plt.plot(N, 1/np.sqrt(N)*factors["datcombine"][0]/factors["wavg"][0], "k--")
plt.xlabel('N')
plt.ylabel('average scaling')
plt.legend()
plt.savefig("output/datcombine/comparison.png")
plt.show()