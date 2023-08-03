import matplotlib.pyplot as plt
import numpy as np

atoms = [1001, 1004, 1481, 2955, 4734, 7729, 9488, 11745, 12340, 16640, 43655]
pepsi_saxs = [17, 18, 17, 45, 69, 115, 228, 98, 425, 181, 2247]
pepsi_saxs_unc = [1, 1, 1, 3, 3, 9, 14, 4, 29, 14, 135]
ausaxs_mt = [61, 63, 72, 96, 142, 206, 234, 300, 302, 442, 2171]
ausaxs_mt_unc = [8, 8, 8, 10, 10, 14, 14, 12, 13, 18, 54]
ausaxs_st = [68, 70, 86, 161, 284, 600, 850, 1380, 1406, 2517, 17512]
ausaxs_st_unc = [13, 12, 16, 20, 37, 93, 110, 113, 134, 124, 386]
foxs = [145, 140, 203, 479, 891, 1775, 2419, 3606, 3758, 6420, 40083]
foxs_unc = [12, 10, 11, 16, 16, 57, 64, 82, 94, 152, 655]
crysol = [1182, 1195, 1218, 18425, 1453, 18339, 18353, 18681, 18792, 2468, 4836]
crysol_unc = [21, 24, 29, 714, 39, 229, 329, 232, 457, 55, 179]

atom_mean = [(192155+148408)/2, (148040+118570)/2, (118310+95332)/2, (95106+76358)/2, (76182+60166)/2, (60006+46584)/2, (46460+35298)/2, (35198+26322)/2, (26245+19144)/2, (19072+13542)/2, (13502+9486)/2, (9446+6368)/2, (6338+4182)/2, (4162+2558)/2, (2542+1440)/2, (1428+736)/2]
ausaxs_serial = [1787, 1163, 884, 720, 585, 493, 423, 359, 310, 267, 232, 201, 132, 110, 94, 84]
ausaxs_serial_unc = [18, 11, 6, 4, 1, 3, 5, 8, 4, 6, 3, 3, 1, 1, 1, 1]

plt.figure()
plt.errorbar(atoms, pepsi_saxs, yerr=pepsi_saxs_unc, fmt='.', color='blue', label="Pepsi-SAXS")
plt.errorbar(atoms, foxs, yerr=foxs_unc, fmt='.', color='orange', label="FoXS")
plt.errorbar(atoms, crysol, yerr=crysol_unc, fmt='.', color='cyan', label="CRYSOL")
plt.errorbar(atoms, ausaxs_st, yerr=ausaxs_st_unc, fmt='.', color='pink', label="$AUSAXS_{singlethreaded}$")
plt.errorbar(atoms, ausaxs_mt, yerr=ausaxs_mt_unc, fmt='.', color='purple', label="$AUSAXS_{multithreaded}$")
plt.errorbar(atom_mean, ausaxs_serial, yerr=ausaxs_serial_unc, fmt='.', color='brown', label="$AUSAXS_{serial}$")
plt.plot(atoms, pepsi_saxs, '--', color='blue')
plt.plot(atoms, foxs, '--', color='orange')
plt.plot(atoms, crysol, '--', color='cyan')
plt.plot(atoms, ausaxs_st, '--', color='pink')
plt.plot(atoms, ausaxs_mt, '--', color='purple')
plt.plot(atom_mean, ausaxs_serial, '--', color='brown')
plt.legend()
plt.xlabel("Number of atoms")
plt.ylabel("Execution time (ms)")
plt.semilogy()
plt.axis([0, 50000, 1, 100000])
plt.savefig("benchmark.png", dpi=300)