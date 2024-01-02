import matplotlib.pyplot as plt
import numpy as np
import os
import json

data = {
    "pepsi": [],
    "foxs": [],
    "crysol": [],
    "ausaxs": [],
    "ausaxs-st": [],
    "ausaxs-serial": []
}

atoms = {
    "SASDPT4": 2470,
    "SASDPP4": 1860,
    "SASDPS4": 3459,
    "SASDE35": 2917,
    "SASDQ59": 3540,
    "SASDJG5": 4734,
    "SASDJG5_serial": 4734,
    "SASDME4": 5354,
    "SASDPB9": 7715,
    "urateox": 9436,
    "SASDDD3": 12043,
    "SASDPR4": 12248,
    "SASDEL9": 16640,
    "SASDPQ4": 24180,
    "A2M_native": 43652
}

for _, dirs, _ in os.walk("test/benchmarks"):
    for dir in dirs:
        for _, _, files in os.walk("test/benchmarks/" + dir):
            files.sort()
            for file in files:
#                print(f"Processing {dir}/{file}")
                jfile = json.load(open("test/benchmarks/" + dir + "/" + file))["results"][0]
                data[jfile["command"]].append([atoms[dir], jfile["mean"]*1e3, jfile["stddev"]*1e3])


# atoms = [1001, 1004, 1481, 2955, 4734, 7729, 9488, 11745, 12340, 16640, 43655]
# pepsi_saxs = [17, 18, 17, 45, 69, 115, 228, 98, 425, 181, 2247]
# pepsi_saxs_unc = [1, 1, 1, 3, 3, 9, 14, 4, 29, 14, 135]
# ausaxs_mt = [61, 63, 72, 96, 142, 206, 234, 300, 302, 442, 2171]
# ausaxs_mt_unc = [8, 8, 8, 10, 10, 14, 14, 12, 13, 18, 54]
# ausaxs_st = [68, 70, 86, 161, 284, 600, 850, 1380, 1406, 2517, 17512]
# ausaxs_st_unc = [13, 12, 16, 20, 37, 93, 110, 113, 134, 124, 386]
# foxs = [145, 140, 203, 479, 891, 1775, 2419, 3606, 3758, 6420, 40083]
# foxs_unc = [12, 10, 11, 16, 16, 57, 64, 82, 94, 152, 655]
# crysol = [1182, 1195, 1218, 18425, 1453, 18339, 18353, 18681, 18792, 2468, 4836]
# crysol_unc = [21, 24, 29, 714, 39, 229, 329, 232, 457, 55, 179]

atom_mean = [(192155+148408)/2, (148040+118570)/2, (118310+95332)/2, (95106+76358)/2, (76182+60166)/2, (60006+46584)/2, (46460+35298)/2, (35198+26322)/2, (26245+19144)/2, (19072+13542)/2, (13502+9486)/2, (9446+6368)/2, (6338+4182)/2, (4162+2558)/2, (2542+1440)/2, (1428+736)/2]
# ausaxs_serial = [1787, 1163, 884, 720, 585, 493, 423, 359, 310, 267, 232, 201, 132, 110, 94, 84]
# ausaxs_serial_unc = [18, 11, 6, 4, 1, 3, 5, 8, 4, 6, 3, 3, 1, 1, 1, 1]

pepsi = np.array(data["pepsi"])
foxs = np.array(data["foxs"])
crysol = np.array(data["crysol"])
ausaxs = np.array(data["ausaxs"])
ausaxs_st = np.array(data["ausaxs-st"])
ausaxs_ser = np.array(data["ausaxs-serial"])

pepsi.sort(axis=0)
foxs.sort(axis=0)
crysol.sort(axis=0)
ausaxs.sort(axis=0)
ausaxs_st.sort(axis=0)
ausaxs_ser[:, 1] = ausaxs_ser[:, 1] / 100
ausaxs_ser[:, 2] = ausaxs_ser[:, 2] / 100

plt.figure()
if data["pepsi"]:           plt.errorbar(pepsi[:, 0],       pepsi[:, 1],        yerr=pepsi[:, 2],     fmt='.', color='blue',      label="Pepsi-SAXS")
if data["foxs"]:            plt.errorbar(foxs[:, 0],        foxs[:, 1],         yerr=foxs[:, 2],      fmt='.', color='orange',    label="FoXS")
if data["crysol"]:          plt.errorbar(crysol[:, 0],      crysol[:, 1],       yerr=crysol[:, 2],    fmt='.', color='cyan',      label="CRYSOL")
if data["ausaxs-st"]:       plt.errorbar(ausaxs_st[:, 0],   ausaxs_st[:, 1],    yerr=ausaxs_st[:, 2], fmt='.', color='pink',      label="$AUSAXS_{singlethreaded}$")
if data["ausaxs"]:          plt.errorbar(ausaxs[:, 0],      ausaxs[:, 1],       yerr=ausaxs[:, 2],    fmt='.', color='purple',    label="$AUSAXS_{multithreaded}$")
if data["ausaxs-serial"]:   plt.errorbar(atom_mean,         ausaxs_ser[:,1],    yerr=ausaxs_ser[:, 2],fmt='.', color='brown',     label="$AUSAXS_{serial}$")
if data["pepsi"]:           plt.plot(pepsi[:, 0],       pepsi[:, 1],     '--', color='blue')
if data["foxs"]:            plt.plot(foxs[:, 0],        foxs[:, 1],      '--', color='orange')
if data["crysol"]:          plt.plot(crysol[:, 0],      crysol[:, 1],    '--', color='cyan')
if data["ausaxs-st"]:       plt.plot(ausaxs_st[:, 0],   ausaxs_st[:, 1], '--', color='pink')
if data["ausaxs"]:          plt.plot(ausaxs[:, 0],      ausaxs[:, 1],    '--', color='purple')
if data["ausaxs-serial"]:   plt.plot(atom_mean,         ausaxs_ser[:, 1],'--', color='brown')
plt.legend()
plt.xlabel("Number of atoms")
plt.ylabel("Execution time (ms)")
plt.semilogy()
plt.axis([0, 50000, 1, 100000])
plt.savefig("benchmark.png", dpi=300)
