import numpy as np
import matplotlib.pyplot as plt


for n in [10,50,100,500]:
    fname = "results/beam_res_%d_1e-10" % n
    r = np.linspace(0,1,n+2)
    j = np.arange(1,n+1)
    h = 1./(n+1)
    e_vals_analytic = 2/h**2-2*np.cos(j*np.pi/(n+1))/h**2
    infile = open(fname,"r")
    infile.readline()
    infile.readline()            #Skip first two lines
    e_vals = []
    v0 = [0.]
    v1 = [0.]
    v2 = [0.]
    for line in infile:
        data = line.split()
        e_vals.append(float(data[0]))
        v0.append(float(data[2]))
        v1.append(float(data[3]))
        v2.append(float(data[4]))
    v0.append(0)
    v1.append(0)
    v2.append(0)
    plt.figure()
    plt.plot(j,e_vals-e_vals_analytic,"x")
    plt.legend()
    plt.xlabel("Eigenvalue number")
    plt.ylabel("Numeric-analytic eigenvalue")
    plt.savefig("results/beam_evals_%d.pdf" % n)

    plt.figure()
    plt.plot(r,v0,label="v0")
    plt.plot(r,v1,label="v1")
    plt.plot(r,v2,label="v2")
    plt.legend()
    plt.xlabel(r"$\rho = x/L$")
    plt.ylabel(r"u($\rho$)")
    plt.savefig("results/beam_evecs_%d.pdf" % n)
