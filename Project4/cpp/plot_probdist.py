"""
This program reads in the number of times each energy has been observed,
and makes one histogram for each temperature
"""

import numpy as np
import matplotlib.pyplot as plt


#Arrays to contain all observed energies:
energies = {"1_ordered":[], "1_random":[], "24_ordered":[], "24_random":[]}

for key in energies:
    infile = open("results/probdist_T%s_10000000.txt" % key)
    for i in range(3):
        infile.readline()
    for line in infile:
        data = line.split()
        new_energies = int(data[1])*[float(data[0])]
        energies[key].extend(new_energies)
    infile.close()

plt.figure()
plt.hist([energies["1_ordered"], energies["1_random"]], \
label=["Ordered", "Random"])
plt.legend(fontsize="14")
plt.title("T = 1")
plt.xlabel("Energy per spin",fontsize="14")
plt.ylabel("Counts",fontsize="14")
plt.tight_layout()
plt.savefig("results/plot_probdist_T=1.pdf")

plt.figure()
plt.hist([energies["24_ordered"], energies["24_random"]], \
bins=20,label=["Ordered", "Random"])
plt.legend(fontsize="14")
plt.title("T = 2.4")
plt.xlabel("Energy per spin",fontsize="14")
plt.ylabel("Counts",fontsize="14")
plt.tight_layout()
plt.savefig("results/plot_probdist_T=2,4.pdf")
