import numpy as np
import matplotlib.pyplot as plt

wrlist = [0.01, 0.5, 1, 5]
n = 300
tol = 1e-10
max_r = [50.,10.,10.,10.]
for i in range(4):
    fname = "results/twoparticles_%d_%g_%g_%d" % (n,tol,max_r[i],i)
    r = np.linspace(0,max_r[i],n+2)
    h = max_r[i]/(n+1)
    v0 = np.zeros(n+2)
    e_vals, v0[1:-1] = np.loadtxt(fname,skiprows=2,unpack=True)
    plt.figure()
    plt.plot(r,v0**2)
    plt.title(r"$\omega_r$ = %g" % (wrlist[i]))
    plt.xlabel(r"r/$\alpha$")
    plt.ylabel(r"$\propto r^2R^2$")
    plt.axis([0,max_r[i],0,1.1*max(v0**2)])
    plt.savefig("results/twoparticles_gs_%d.pdf" % i)
