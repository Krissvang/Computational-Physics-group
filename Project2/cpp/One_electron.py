import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


#Definition of variable
n=910
tol = 1e-10
max_r = 5.03
#Definition of the interval
h = max_r/(n+1)
#Definition of the position vector
r = np.linspace(0,max_r,n+2)
#Here we read the finel results frome the file
fname = "Results_one_electron/final_result"
v0 = np.zeros(n+2)
v1 = np.zeros(n+2)
v2 = np.zeros(n+2)
e_vals, v0[1:-1], v1[1:-1], v2[1:-1] = np.loadtxt(fname,skiprows=4,unpack=True)
#Here we open an image
plt.figure()
#Plot of the wave functions for the fist three eigenvalues
eigenvalue_0,=plt.plot(r,v0**2, label='eigenvalue 0')
eigenvalue_1,=plt.plot(r,v1**2, label='eigenvalue 1')
eigenvalue_2,=plt.plot(r,v2**2, label='eigenvalue 2')
#Settings of the plot property
plt.legend([eigenvalue_0, eigenvalue_1, eigenvalue_2], ['eigenvalue 0', 'eigenvalue 1', 'eigenvalue 2'])
plt.legend()
plt.title(r"Wave function for one electron in a harmonic oscillator potential ")
plt.xlabel(r"r/$\alpha$")
plt.ylabel(r"$\phi$")
plt.axis([0,max_r,0,1.1*max(v0**2)])
plt.ylim((-0.0004,0.005))
plt.xlim((-0.3,5.3))
plt.savefig("Results_one_electron/Plot.pdf")


#Definition of variables for the other plots
file_values = np.array([0,1,2,3,4,5,6,7,8,9])
r_values = np.array([4,6,8,10,12,14,16,18,20,22])
eigen_0 = np.zeros(10)
eigen_1 = np.zeros(10)
eigen_2= np.zeros(10)

r_values1 = np.array([4,4.2,4.4,4.6,4.8,5,5.2,5.4,5.6,5.8])
eigen1_0 = np.zeros(10)
eigen1_1 = np.zeros(10)
eigen1_2= np.zeros(10)

n_values = np.array([50,100,150,200,250,300,400,500,600,800])
eigen2_0 = np.zeros(10)
eigen2_1 = np.zeros(10)
eigen2_2= np.zeros(10)

#Here we read the finel results frome the file
for file_val in file_values:
    fname = "Results_one_electron/n_constant_%d" % file_val
    fname1 = "Results_one_electron/n_constant1_%d" % file_val
    fname2 = "Results_one_electron/final_result_%d" % file_val
    infile = open(fname,"r")
    lines=infile.readlines()
    eigen_0[file_val]=float(lines[4].split()[0])
    eigen_1[file_val]=float(lines[5].split()[0])
    eigen_2[file_val]=float(lines[6].split()[0])
    
    infile1 = open(fname1,"r")
    lines1=infile1.readlines()
    eigen1_0[file_val]=float(lines1[4].split()[0])
    eigen1_1[file_val]=float(lines1[5].split()[0])
    eigen1_2[file_val]=float(lines1[6].split()[0])

    infile2 = open(fname1,"r")
    lines2=infile2.readlines()
    eigen2_0[file_val]=float(lines2[4].split()[0])
    eigen2_1[file_val]=float(lines2[5].split()[0])
    eigen2_2[file_val]=float(lines2[6].split()[0])

#Here we open an image   
plt.figure()
#Here we plot the data
eigenvalue_0,=plt.plot(r_values,eigen_0-3, "b.", label='eigenvalue 0')
eigenvalue_1,=plt.plot(r_values,eigen_1-7, "r.", label='eigenvalue 1')
eigenvalue_2,=plt.plot(r_values,eigen_2-11, "g.",label='eigenvalue 2')
plt.plot(r_values,eigen_0-3, 'b')
plt.plot(r_values,eigen_1-7, 'r')
plt.plot(r_values,eigen_2-11, 'g')
#Settings of the plot property
plt.legend([eigenvalue_0, eigenvalue_1, eigenvalue_2], ['eigenvalue 0', 'eigenvalue 1', 'eigenvalue 2'])
plt.legend()
plt.title(r"Difference between the computed and the theoretical eigenvalues")
plt.ylabel(r"a.u.")
plt.xlabel(r"$\rho_{max}$")
plt.savefig("Results_one_electron/Plot_1.pdf")

#Here we open an image
plt.figure()
#Here we plot the data
eigenvalue_0,=plt.plot(r_values1,eigen1_0-3, "b.", label='eigenvalue 0')
eigenvalue_1,=plt.plot(r_values1,eigen1_1-7, "r.", label='eigenvalue 1')
eigenvalue_2,=plt.plot(r_values1,eigen1_2-11, "g.", label='eigenvalue 2')
plt.plot(r_values1,eigen1_0-3, 'b')
plt.plot(r_values1,eigen1_1-7, 'r')
plt.plot(r_values1,eigen1_2-11, 'g')
#Settings of the plot property
plt.legend([eigenvalue_0, eigenvalue_1, eigenvalue_2], ['eigenvalue 0', 'eigenvalue 1', 'eigenvalue 2'])
plt.legend()
plt.title(r"Difference between the computed and the theoretical eigenvalues")
plt.ylabel(r"a.u.")
plt.xlabel(r"$\rho_{max}$")
plt.savefig("Results_one_electron/Plot_2.pdf")

#Here we open an image
plt.figure()
#Here we plot the data
eigenvalue_0,=plt.plot(n_values,eigen2_0-3, "b.", label='eigenvalue 0')
eigenvalue_1,=plt.plot(n_values,eigen2_1-7, "r.", label='eigenvalue 1')
eigenvalue_2,=plt.plot(n_values,eigen2_2-11, "g.", label='eigenvalue 2')
plt.plot(n_values,eigen1_0-3, 'b')
plt.plot(n_values,eigen1_1-7, 'r')
plt.plot(n_values,eigen1_2-11, 'g')
#Settings of the plot property
plt.legend([eigenvalue_0, eigenvalue_1, eigenvalue_2], ['eigenvalue 0', 'eigenvalue 1', 'eigenvalue 2'])
plt.legend()
plt.title(r"Difference between the computed and the theoretical eigenvalues")
plt.ylabel(r"a.u.")
plt.xlabel(r"n")
plt.savefig("Results_one_electron/Plot_3.pdf")


