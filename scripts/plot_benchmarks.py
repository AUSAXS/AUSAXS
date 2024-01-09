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
    "SASDPT4": 1001,
    "SASDPP4": 1004,
    "SASDPS4": 1481,
    "SASDE35": 2917,
    "SASDQ59": 3540,
    "SASDJG5": 4734,
    "SASDJG5_serial": 4734,
    "SASDME4": 5380,
    "SASDPB9": 6210,
    "urateox": 9436,
    "SASDPQ4": 9488,
    "SASDDD3": 10775,
    "SASDPR4": 12332,
    "SASDEL9": 16640,
    "A2M_native": 43652
}

for _, dirs, _ in os.walk("test/benchmarks"):
    for dir in dirs:
        print("")
        for _, _, files in os.walk("test/benchmarks/" + dir):
            files.sort()
            for file in files:
#                print(f"Processing {dir}/{file}")
                jfile = json.load(open("test/benchmarks/" + dir + "/" + file))["results"][0]
                data[jfile["command"]].append([atoms[dir], jfile["mean"]*1e3, jfile["stddev"]*1e3])
                print(f"{dir} {jfile['command']} {jfile['mean']*1e3} {jfile['stddev']*1e3}")

atom_mean = [(192155+148408)/2, (148040+118570)/2, (118310+95332)/2, (95106+76358)/2, (76182+60166)/2, (60006+46584)/2, (46460+35298)/2, (35198+26322)/2, (26245+19144)/2, (19072+13542)/2, (13502+9486)/2, (9446+6368)/2, (6338+4182)/2, (4162+2558)/2, (2542+1440)/2, (1428+736)/2]

pepsi = np.array(data["pepsi"])
foxs = np.array(data["foxs"])
crysol = np.array(data["crysol"])
ausaxs = np.array(data["ausaxs"])
ausaxs_st = np.array(data["ausaxs-st"])
ausaxs_ser = np.array(data["ausaxs-serial"])

pepsi = pepsi[pepsi[:, 0].argsort()]
foxs = foxs[foxs[:, 0].argsort()]
crysol = crysol[crysol[:, 0].argsort()]
ausaxs = ausaxs[ausaxs[:, 0].argsort()]
ausaxs_st = ausaxs_st[ausaxs_st[:, 0].argsort()]
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
plt.savefig("benchmark.png", dpi=600)
