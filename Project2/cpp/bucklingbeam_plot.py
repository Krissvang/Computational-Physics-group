import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

nvalues = np.array([5,10,20,50,100,200,500])
times =  []
iterations = []

for n in nvalues:
    fname = "results/beam_res_%d_1e-10" % n
    r = np.linspace(0,1,n+2)
    j = np.arange(1,n+1)
    h = 1./(n+1)
    e_vals_analytic = 2/h**2-2*np.cos(j*np.pi/(n+1))/h**2
    infile = open(fname,"r")
    line = infile.readline()
    data = line.split()
    #Read time and number of iterations from file, remove "("
    times.append(float(data[1]))
    iterations.append(int((data[2])[1:]))
    infile.readline()                     #Skip second line
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
    plt.title("n = %d" %n)
    plt.xlabel("Eigenvalue number")
    plt.ylabel("Numeric-analytic eigenvalue")
    plt.savefig("results/beam_evals_%d.pdf" % n)

    plt.figure()
    plt.plot(r,v0,label="v0")
    plt.plot(r,v1,label="v1")
    plt.plot(r,v2,label="v2")
    plt.title("n = %d" %n)
    plt.legend()
    plt.xlabel(r"$\rho = x/L$")
    plt.ylabel(r"u($\rho$)")
    plt.savefig("results/beam_evecs_%d.pdf" % n)

times = np.array(times)
iterations =np.array(iterations)
x = np.log10(nvalues)
y1 = np.log10(times)
y2 = np.log10(iterations)
plt.figure()
plt.plot(x,y1,"x",label="times")
plt.plot(x,y2,"x",label="iterations")
#Fit data to straight lines:
m1, c1, r1, p1, s1 = linregress(x,y1)
m2, c2, r2, p2, s2 = linregress(x,y2)
plt.plot(x,m1*x+c1,label="%gx+%g"%(m1,c1))
plt.plot(x,m2*x+c2,label="%gx+%g"%(m2,c2))
plt.legend()
plt.xlabel("log10(n)")
plt.ylabel("log10(times or iterations)")
plt.savefig("results/beam_times-and-iterations.pdf")
