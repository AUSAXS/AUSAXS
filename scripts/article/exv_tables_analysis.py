import matplotlib.pyplot as plt
import numpy as np
from enum import Enum
import pandas as pd
import sys

class Data(Enum):
    chi2 = 0
    H = 1
    C = 2
    CH = 3
    CH2 = 4
    CH3 = 5
    N = 6
    NH = 7
    NH2 = 8
    NH3 = 9
    O = 10
    OH = 11
    S = 12
    SH = 13

file = sys.argv[1]
target = sys.argv[2]
with open(file, 'r') as f:
    lines = f.readlines()

data_voronoi = []
data_mf = []
section = -1
for line in lines:
    line = line.strip()

    if not line:
        continue
    if line.startswith("#"):
        section += 1
        continue
    if line.startswith("chi2"):
        continue

    tokens = line.split()
    if len(tokens) != 14:
        print(f"{line}\nCRITICAL ERROR")
        exit(1)
    tokens = [float(t) for t in tokens]
    match (section):
        case 0: 
            data_voronoi.append(tokens)
        case 1:
            data_mf.append(tokens)

data_voronoi = np.array(data_voronoi)
data_mf = np.array(data_mf)
data = data_voronoi if target == 'voronoi' else data_mf

##########################
### CORRELATION MATRIX ###
##########################
if (True):
    def get_correlation_matrix(data):
        df = pd.DataFrame(data, columns=[e.name for e in Data])
        correlation_matrix = df[['CH', 'CH2', 'CH3', 'NH', 'NH2', 'NH3', 'O', 'OH', 'S', 'SH', 'chi2']].corr()
        correlation_matrix[correlation_matrix[:] == 1] = np.NaN
        for i in range(len(correlation_matrix.columns)):
            for j in range(i, len(correlation_matrix.columns)):
                correlation_matrix.iloc[i, j] = np.NaN
        return correlation_matrix

    # df = pd.DataFrame(data, columns=[e.name for e in Data])
    # correlation_matrix = df[['CH', 'CH2', 'CH3', 'NH', 'NH2', 'NH3', 'O', 'OH', 'S', 'SH', 'chi2']].corr()
    # correlation_matrix[correlation_matrix[:] == 1] = np.NaN
    # for i in range(len(correlation_matrix.columns)):
    #     for j in range(i, len(correlation_matrix.columns)):
    #         correlation_matrix.iloc[i, j] = np.NaN
    # plt.figure(figsize=(8, 6))
    # cmap = plt.get_cmap('coolwarm', 20)
    # cmap.set_bad(color='white')
    # plt.imshow(correlation_matrix, cmap=cmap, vmin=-0.1, vmax=0.5, origin='upper')
    # plt.colorbar(label='Correlation Coefficient')
    # plt.xticks(ticks=np.arange(len(correlation_matrix.columns)), labels=correlation_matrix.columns, rotation=45)
    # plt.yticks(ticks=np.arange(len(correlation_matrix.columns)), labels=correlation_matrix.columns)
    # plt.title('Correlation Matrix')
    # plt.show()

    plt.figure(figsize=(12, 12))
    plt.subplots_adjust(hspace=0.35, wspace=0)
    plt.subplot(2, 1, 1)
    bins = np.linspace(0, 200, 200)
    plt.hist(data_voronoi[:, Data.chi2.value], bins=bins, color='tab:blue', alpha=0.5)
    plt.hist(data_mf[:, Data.chi2.value], bins=bins, color='tab:orange', alpha=0.5)

    plt.axvline(55.86, color='tab:blue', linestyle='--', label="Voronoi")
    plt.axvline(5.75, color='tab:orange', linestyle='--', label="Minimum Fluctuation")
    plt.axvline(10.18, color='tab:red', linestyle='--', label="Traube")

    plt.yscale('log')
    plt.xlabel('$\chi^2_r$')
    plt.ylabel('Density')
    plt.legend()

    corr = get_correlation_matrix(data_voronoi)
    labels = corr.columns.tolist()
    labels[-1] = '$\chi^2_r$'
    ax1 = plt.subplot(2, 2, 3)
    im1 = plt.imshow(corr, cmap='coolwarm', vmin=-0.1, vmax=0.5, origin='upper')
    plt.xticks(ticks=np.arange(len(labels)), labels=labels, rotation=90)
    plt.yticks(ticks=np.arange(len(labels)), labels=labels)
    plt.title('Voronoi')
    ax1.xaxis.set_ticks_position('none')
    ax1.yaxis.set_ticks_position('none')

    corr = get_correlation_matrix(data_mf)
    ax2 = plt.subplot(2, 2, 4, sharey=ax1)
    im2 = plt.imshow(corr, cmap='coolwarm', vmin=-0.1, vmax=0.5, origin='upper')
    plt.xticks(ticks=np.arange(len(labels)), labels=labels, rotation=90)
    plt.yticks(ticks=np.arange(len(labels)), labels=labels)
    plt.title('Minimum Fluctuation')
    ax2.xaxis.set_ticks_position('none')
    ax2.yaxis.set_visible(False)
#    plt.colorbar(label='Correlation coefficient')

    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im2, cax=cax, label='Correlation coefficient')

    plt.savefig("exv_table_analysis.png", dpi=600)

    print(
        "Largest outliers: " \
        f"\n\tVoronoi: {data_voronoi[:, Data.chi2.value].max()}" \
        f"\n\tMinimum Fluctuation: {data_mf[:, Data.chi2.value].max()}"
    )

    print(
        "Percentage higher than 100: " \
        f"\n\tVoronoi: {100*np.sum(100 < data_voronoi[:, Data.chi2.value])/len(data_voronoi[:, Data.chi2.value]):.1f}%" \
        f"\n\tMinimum Fluctuation: {100*np.sum(100 < data_mf[:, Data.chi2.value])/len(data_mf[:, Data.chi2.value]):.1f}%" \
    )
    plt.show()

#######################
### 1D CORRELATIONS ###
#######################
means = np.mean(data, axis=0)
stds = np.std(data, axis=0)
if (False):
    alpha = 0.1
    fig, ax = plt.subplots(4, 3, figsize=(14, 10), sharey='row')
    plt.sca(ax[0, 0])
    plt.plot(data[:, Data.CH.value], data[:, Data.chi2.value], '.k', alpha=alpha)
    plt.axvline(means[Data.CH.value] - stds[Data.CH.value], color='r', linestyle='--')
    plt.axvline(means[Data.CH.value] + stds[Data.CH.value], color='r', linestyle='--')

    plt.ylabel('$\chi^2_r$')
    plt.sca(ax[0, 1])
    plt.plot(data[:, Data.CH2.value], data[:, Data.chi2.value], '.k', alpha=alpha)
    plt.axvline(means[Data.CH2.value] - stds[Data.CH2.value], color='r', linestyle='--')
    plt.axvline(means[Data.CH2.value] + stds[Data.CH2.value], color='r', linestyle='--')
    plt.sca(ax[0, 2])
    plt.plot(data[:, Data.CH3.value], data[:, Data.chi2.value], '.k', alpha=alpha)    
    plt.axvline(means[Data.CH3.value] - stds[Data.CH3.value], color='r', linestyle='--')
    plt.axvline(means[Data.CH3.value] + stds[Data.CH3.value], color='r', linestyle='--')

    plt.sca(ax[1, 0])
    plt.plot(data[:, Data.NH.value], data[:, Data.chi2.value], '.k', alpha=alpha)
    plt.axvline(means[Data.NH.value] - stds[Data.NH.value], color='r', linestyle='--')
    plt.axvline(means[Data.NH.value] + stds[Data.NH.value], color='r', linestyle='--')
    plt.ylabel('$\chi^2_r$')
    plt.sca(ax[1, 1])
    plt.plot(data[:, Data.NH2.value], data[:, Data.chi2.value], '.k', alpha=alpha)
    plt.axvline(means[Data.NH2.value] - stds[Data.NH2.value], color='r', linestyle='--')
    plt.axvline(means[Data.NH2.value] + stds[Data.NH2.value], color='r', linestyle='--')
    plt.sca(ax[1, 2])
    plt.plot(data[:, Data.NH3.value], data[:, Data.chi2.value], '.k', alpha=alpha)
    plt.axvline(means[Data.NH3.value] - stds[Data.NH3.value], color='r', linestyle='--')
    plt.axvline(means[Data.NH3.value] + stds[Data.NH3.value], color='r', linestyle='--')
    plt.xlabel('Volume [A$^3$]')

    plt.sca(ax[2, 0])
    plt.plot(data[:, Data.O.value], data[:, Data.chi2.value], '.k', alpha=alpha)
    plt.axvline(means[Data.O.value] - stds[Data.O.value], color='r', linestyle='--')
    plt.axvline(means[Data.O.value] + stds[Data.O.value], color='r', linestyle='--')
    plt.ylabel('$\chi^2_r$')
    plt.sca(ax[2, 1])
    plt.plot(data[:, Data.OH.value], data[:, Data.chi2.value], '.k', alpha=alpha)
    plt.axvline(means[Data.OH.value] - stds[Data.OH.value], color='r', linestyle='--')
    plt.axvline(means[Data.OH.value] + stds[Data.OH.value], color='r', linestyle='--')
    fig.delaxes(ax[2, 2])

    plt.sca(ax[3, 0])
    plt.plot(data[:, Data.S.value], data[:, Data.chi2.value], '.k', alpha=alpha)
    plt.axvline(means[Data.S.value] - stds[Data.S.value], color='r', linestyle='--')
    plt.axvline(means[Data.S.value] + stds[Data.S.value], color='r', linestyle='--')
    plt.xlabel('Volume [A$^3$]')
    plt.ylabel('$\chi^2_r$')
    plt.sca(ax[3, 1])
    plt.plot(data[:, Data.SH.value], data[:, Data.chi2.value], '.k', alpha=alpha)
    plt.axvline(means[Data.SH.value] - stds[Data.SH.value], color='r', linestyle='--')
    plt.axvline(means[Data.SH.value] + stds[Data.SH.value], color='r', linestyle='--')
    plt.xlabel('Volume [A$^3$]')
    fig.delaxes(ax[3, 2])

    plt.text(0.5, 1.1, '1H', ha='center', va='center', transform=ax[0, 0].transAxes)
    plt.text(0.5, 1.1, '2H', ha='center', va='center', transform=ax[0, 1].transAxes)
    plt.text(0.5, 1.1, '3H', ha='center', va='center', transform=ax[0, 2].transAxes)

    plt.text(-0.3, 0.5, 'C', ha='center', va='center', transform=ax[0, 0].transAxes, rotation=90)
    plt.text(-0.3, 0.5, 'N', ha='center', va='center', transform=ax[1, 0].transAxes, rotation=90)
    plt.text(-0.3, 0.5, 'O', ha='center', va='center', transform=ax[2, 0].transAxes, rotation=90)
    plt.text(-0.3, 0.5, 'S', ha='center', va='center', transform=ax[3, 0].transAxes, rotation=90)

    plt.tight_layout()
    plt.show()

#######################
### 2D CORRELATIONS ###
#######################
if (False):
    def perform_plot(datax, datay):
        num_bins = 30
        H_sum, xedges, yedges = np.histogram2d(datax, datay, bins=num_bins, weights=data[:, Data.chi2.value])
        H_count, _, _ = np.histogram2d(datax, datay, bins=num_bins)

        H_mean = np.divide(H_sum, H_count, out=np.zeros_like(H_sum), where=H_count != 0)
        H_mean[H_count == 0] = np.nan

        plt.imshow(H_mean.T, origin='lower', aspect='auto', 
                extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                cmap='viridis')

    use_data = [Data.CH.value, Data.CH2.value, Data.CH3.value, Data.NH.value, Data.NH2.value, Data.NH3.value, Data.O.value, Data.OH.value, Data.S.value, Data.SH.value]
    fig, ax = plt.subplots(10, 10, figsize=(20, 20), sharex='col', sharey='row')
    for i in range(10):
        for j in range(i):
            plt.sca(ax[i, j])
            perform_plot(data[:, use_data[j]], data[:, use_data[i]])
        
        plt.text(0.5, -0.4, str(Data(use_data[i]).name), ha='center', va='center', transform=ax[-1, i].transAxes)
        plt.text(-0.3, 0.5, str(Data(use_data[i]).name), ha='center', va='center', transform=ax[i, 0].transAxes, rotation=90)

    plt.tight_layout()
    plt.show()
