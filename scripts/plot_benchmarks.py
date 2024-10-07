import matplotlib.pyplot as plt
import numpy as np
import os
import json

data = {
    "pepsi": [],
    "foxs": [],
    "crysol": [],
    "ausaxs_simple": [],
    "ausaxs_fraser": [],
    "ausaxs_grid": [],
    "ausaxs_simple_st": [],
    "ausaxs_serial": []
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
    "SASDA92": 16068,
    "SASDEL9": 16640,
    "A2M_native": 43652
}

for _, dirs, _ in os.walk("tests/benchmarks"):
    for dir in dirs:
        for _, _, files in os.walk("tests/benchmarks/" + dir):
            files.sort()
            for file in files:
                print(f"Processing {dir}/{file}")
                jfile = json.load(open("tests/benchmarks/" + dir + "/" + file))["results"][0]
                key = file.split(".")[0]
                print(f"key is {key}")
                data[key].append([atoms[dir], jfile["mean"]*1e3, jfile["stddev"]*1e3])
                print(f"{dir} {jfile['command']} {jfile['mean']*1e3} {jfile['stddev']*1e3}")

atom_mean = [(192155+148408)/2, (148040+118570)/2, (118310+95332)/2, (95106+76358)/2, (76182+60166)/2, (60006+46584)/2, (46460+35298)/2, (35198+26322)/2, (26245+19144)/2, (19072+13542)/2, (13502+9486)/2, (9446+6368)/2, (6338+4182)/2, (4162+2558)/2, (2542+1440)/2, (1428+736)/2]

pepsi               = np.array(data["pepsi"])
foxs                = np.array(data["foxs"])
crysol              = np.array(data["crysol"])
ausaxs_simple       = np.array(data["ausaxs_simple"])
ausaxs_fraser       = np.array(data["ausaxs_fraser"])
ausaxs_grid         = np.array(data["ausaxs_grid"])
ausaxs_simple_st    = np.array(data["ausaxs_simple_st"])
ausaxs_ser          = np.array(data["ausaxs_serial"])

if data["pepsi"]:            pepsi               = pepsi[pepsi[:, 0].argsort()]
if data["foxs"]:             foxs                = foxs[foxs[:, 0].argsort()]
if data["crysol"]:           crysol              = crysol[crysol[:, 0].argsort()]
if data["ausaxs_simple"]:    ausaxs_simple       = ausaxs_simple[ausaxs_simple[:, 0].argsort()]
if data["ausaxs_fraser"]:    ausaxs_fraser       = ausaxs_fraser[ausaxs_fraser[:, 0].argsort()]
if data["ausaxs_grid"]:      ausaxs_grid         = ausaxs_grid[ausaxs_grid[:, 0].argsort()]
if data["ausaxs_simple_st"]: ausaxs_simple_st    = ausaxs_simple_st[ausaxs_simple_st[:, 0].argsort()]
if data["ausaxs_serial"]:    ausaxs_ser[:, 1]    = ausaxs_ser[:, 1] / 100
if data["ausaxs_serial"]:    ausaxs_ser[:, 2]    = ausaxs_ser[:, 2] / 100

plt.figure()
if data["pepsi"]:           plt.errorbar(pepsi[:, 0],               pepsi[:, 1],                yerr=pepsi[:, 2],            fmt='.', color='blue',      label="Pepsi-SAXS")
if data["foxs"]:            plt.errorbar(foxs[:, 0],                foxs[:, 1],                 yerr=foxs[:, 2],             fmt='.', color='orange',    label="FoXS")
if data["crysol"]:          plt.errorbar(crysol[:, 0],              crysol[:, 1],               yerr=crysol[:, 2],           fmt='.', color='cyan',      label="CRYSOL")
if data["ausaxs_simple"]:   plt.errorbar(ausaxs_simple[:, 0],       ausaxs_simple[:, 1],        yerr=ausaxs_simple[:, 2],    fmt='.', color='purple',    label="Simple")
if data["ausaxs_simple_st"]:plt.errorbar(ausaxs_simple_st[:, 0],    ausaxs_simple_st[:, 1],     yerr=ausaxs_simple_st[:, 2], fmt='.', color='pink',      label="Simple$_{st}$")
if data["ausaxs_fraser"]:   plt.errorbar(ausaxs_fraser[:, 0],       ausaxs_fraser[:, 1],        yerr=ausaxs_fraser[:, 2],    fmt='.', color='green',     label="Fraser")
if data["ausaxs_grid"]:     plt.errorbar(ausaxs_grid[:, 0],         ausaxs_grid[:, 1],          yerr=ausaxs_grid[:, 2],      fmt='.', color='red',       label="Grid-based")
if data["ausaxs_serial"]:   plt.errorbar(atom_mean,                 ausaxs_ser[:,1],            yerr=ausaxs_ser[:, 2],       fmt='.', color='brown',     label="Serial")

if data["pepsi"]:           plt.plot(pepsi[:, 0],               pepsi[:, 1],            '--', color='blue')
if data["foxs"]:            plt.plot(foxs[:, 0],                foxs[:, 1],             '--', color='orange')
if data["crysol"]:          plt.plot(crysol[:, 0],              crysol[:, 1],           '--', color='cyan')
if data["ausaxs_simple"]:   plt.plot(ausaxs_simple[:, 0],       ausaxs_simple[:, 1],    '--', color='purple')
if data["ausaxs_simple_st"]:plt.plot(ausaxs_simple_st[:, 0],    ausaxs_simple_st[:, 1], '--', color='pink')
if data["ausaxs_fraser"]:   plt.plot(ausaxs_fraser[:, 0],       ausaxs_fraser[:, 1],    '--', color='green')
if data["ausaxs_grid"]:     plt.plot(ausaxs_grid[:, 0],         ausaxs_grid[:, 1],      '--', color='red')
if data["ausaxs_serial"]:   plt.plot(atom_mean,                 ausaxs_ser[:, 1],       '--', color='brown')
plt.legend()
plt.xlabel("Number of atoms")
plt.ylabel("Execution time (ms)")
plt.semilogy()
plt.axis([0, 20000, 1, 100000])
plt.savefig("benchmark.png", dpi=600)
