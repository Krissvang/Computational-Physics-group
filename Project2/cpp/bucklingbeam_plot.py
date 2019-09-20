import numpy as np
import matplotlib.pyplot as plt

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
    infile.close()
    v0 =np.zeros(n+2)
    v1 =np.zeros(n+2)
    v2 =np.zeros(n+2)
    e_vals, e_vals_Jacobi, v0[1:-1], v1[1:-1], v2[1:-1] = \
    np.loadtxt(fname,skiprows=2,unpack=True)
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
p1 = np.polyfit(x,y1,1)
p2 = np.polyfit(x,y2,1)
plt.plot(x,p1[0]*x+p1[1],label="%gx+%g"%(p1[0],p1[1]))
plt.plot(x,p2[0]*x+p2[1],label="%gx+%g"%(p2[0],p2[1]))
plt.legend()
plt.xlabel("log10(n)")
plt.ylabel("log10(times or iterations)")
plt.savefig("results/beam_times-and-iterations.pdf")
