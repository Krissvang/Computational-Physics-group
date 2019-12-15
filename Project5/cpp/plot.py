import matplotlib
matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt

#Have to use this for unix subsystem
#This program plots all the plots used in the project


w= "1"

def plot_energy_vs_alpha(omega):
  alpha, pot_e, var_pot_e, kin_e, var_kin_e, r12, Accepted =\
  np.loadtxt("results/energy_vs_alpha_w=%s"%(omega), skiprows=2,unpack=True) 
  plt.figure("Plot")
  plt.grid(1)
  plt.plot(alpha,pot_e,'#E57375',marker='o',linestyle='dashed', markersize='3')
  plt.xlabel("$\\alpha$",fontsize="15")
  plt.ylabel("$E_0$",fontsize="15")
  plt.savefig("results/energy_vs_alpha_w=%s.pdf"%(omega))
  plt.close()

def plot_energy_vs_alpha_var(omega):
  alpha, pot_e, var_pot_e, kin_e, var_kin_e, r12, Accepted =\
  np.loadtxt("results/energy_vs_alpha_w=%s"%(omega), skiprows=2,unpack=True) 
  plt.figure("Plot")
  plt.grid(1)
  plt.plot(alpha,var_pot_e,'#E57375',marker='o',
  linestyle='dashed', markersize='3')
  plt.xlabel("$\\alpha$",fontsize="15")
  plt.ylabel("$\sigma^2_{E_L}$",fontsize="15")
  plt.savefig("results/energy_vs_alpha_var_w=%s.pdf"%(omega))
  plt.close()


def plot_r12_vs_omega(filename,filename2,identifier,identifier2):
  var, pot_e, var_pot_e, kin_e, var_kin_e,pe_wo,\
    var_pe_wo,pe_w,var_pe_w, r12, Accepted =\
  np.loadtxt("results/%s"%(filename), skiprows=2,unpack=True) 
  var1, pot_e1, var_pot_e1, kin_e1, var_kin_e1, pe_wo,\
    var_pe_wo1, pe_w,var_pe_w1, r12_1, Accepted1 =\
  np.loadtxt("results/%s"%(filename2), skiprows=2,unpack=True) 
  plt.figure("Plot")
  plt.grid(1)
  plt.plot(var,r12,'#E57375',marker='o',linestyle='dashed',
   markersize='3',label="With Jastrow")
  plt.plot(var,r12_1,'#4d80e4',marker='o',linestyle='dashed',
   markersize='3', label="Without Jastrow")

  plt.legend(fontsize="12")

  plt.xlabel("$\omega$",fontsize="15")
  plt.ylabel("$\langle r_{12}\\rangle$",fontsize="15")

  plt.savefig("results/r12_vs_omega.pdf")
  plt.close()

def plot_r12_vs_omega_diff(filename,filename2,identifier,identifier2):
  var, pot_e, var_pot_e, kin_e, var_kin_e,pe_wo,\
    var_pe_wo,pe_w,var_pe_w, r12, Accepted =\
  np.loadtxt("results/%s"%(filename), skiprows=2,unpack=True) 
  var1, pot_e1, var_pot_e1, kin_e1, var_kin_e1, pe_wo,\
    var_pe_wo1, pe_w,var_pe_w1, r12_1, Accepted1 =\
  np.loadtxt("results/%s"%(filename2), skiprows=2,unpack=True) 
  plt.figure("Plot")
  plt.grid(1)
  plt.plot(var,r12-r12_1,'#E57375',marker='o',linestyle='dashed',
   markersize='3',label="Difference between $\langle r_{12}\
   \\rangle$ with Jastrow and without")


  plt.legend(fontsize = "11")

  plt.xlabel("$\omega$",fontsize="15")
  plt.ylabel("$\Delta\langle r_{12}\\rangle$",fontsize="15")

  plt.savefig("results/r12_vs_omega_diff.pdf")
  plt.close()

def plot_energy_vs_omega(filename,filename2):
  var, e, var_pot_e, kin_e, var_kin_e,pe_wo,var_pe_wo,pe_w,\
    var_pe_w, r12, Accepted =\
  np.loadtxt("results/%s"%(filename), skiprows=2,unpack=True) 

  var, e1, var_pot_e, kin_e, var_kin_e,pe_wo,var_pe_wo,\
    pe_w,var_pe_w, r12, Accepted =\
  np.loadtxt("results/%s"%(filename2), skiprows=2,unpack=True) 

  plt.figure("Plot")
  plt.grid(1)

  plt.plot(var,e,'#E57375',marker='o',linestyle='dashed',
   markersize='3',label="With Jastrow")
  plt.plot(var,e1,'#4d80e4',marker='o',linestyle='dashed',
   markersize='3',label="Without Jastrow")

  plt.xlabel("$\omega$",fontsize="15")
  plt.ylabel("$\langle E_{min} \\rangle$",fontsize="15")
  plt.legend(fontsize="12")
  plt.savefig("results/energy_vs_omega.pdf")
  plt.close()

def plot_var_vs_omega(filename,filename2):
  var, e, var_e, kin_e, var_kin_e,pe_wo,var_pe_wo,pe_w,var_pe_w, r12, Accepted =\
  np.loadtxt("results/%s"%(filename), skiprows=2,unpack=True) 

  var, e, var_e1, kin_e, var_kin_e,pe_wo,var_pe_wo,pe_w,var_pe_w, r12, Accepted =\
  np.loadtxt("results/%s"%(filename2), skiprows=2,unpack=True) 

  plt.figure("Plot")
  plt.grid(1)
  plt.plot(var,var_e,'#E57375',marker='o',linestyle='dashed',
   markersize='3',label="With Jastrow")
  plt.plot(var,var_e1,'#4d80e4',marker='o',linestyle='dashed',
   markersize='3',label="Without Jastrow")

  plt.xlabel("$\omega$",fontsize="15")
  plt.ylabel("$\sigma^2$",fontsize="15")
  plt.legend(fontsize="12")
  plt.savefig("results/variance_vs_omega.pdf")
  plt.close()


def plot_kin_vs_pot_wo(filename,identifier):
  var, e, var_pot_e, kin_e, var_kin_e,pe_wo,var_pe_wo,\
    pe_w,var_pe_w, r12, Accepted =\
  np.loadtxt("results/%s"%(filename), skiprows=2,unpack=True) 
  plt.figure("Plot")
  plt.grid(1)
  plt.plot(var,kin_e/pe_wo,'#E57375',marker='o',linestyle='dashed', markersize='3')
  plt.xlabel("$\omega$",fontsize="14")
  plt.ylabel("$\\frac{\langle T \\rangle}{\langle V \\rangle}$",fontsize="14")
  plt.savefig("results/ke_vs_pe_unperturbed%s.pdf"%(identifier))
  plt.close()

def plot_kin_vs_pot(filename, filename3):
  var, e, var_pot_e, kin_e, var_kin_e, pe_wo, var_pe_wo,\
     pe_w,var_pe_w, r12, Accepted =\
  np.loadtxt("results/%s"%(filename), skiprows=2,unpack=True) 

  var, e, var_pot_e, kin_e3, var_kin_e3, pe_wo3,var_pe_wo3,\
     pe_w3, var_pe_w3, r12, Accepted =\
  np.loadtxt("results/%s"%(filename3), skiprows=2,unpack=True) 

  plt.figure("Plot")
  plt.grid(1)

  plt.plot(var,kin_e/pe_w,'#E57375',marker='o',linestyle='dashed',
   markersize='3',label="With repulsion")


  plt.plot(var,kin_e3/pe_wo3,'#4d80e4',marker='o',linestyle='dashed',
   markersize='3', label = "Without repultion")

  plt.xlabel("$\omega$",fontsize="14")
  plt.ylabel("$\\frac{\langle T \\rangle}{\langle V \\rangle}$",fontsize="14")
  plt.legend(fontsize="12")
  plt.savefig("results/ke_vs_pe.pdf")
  plt.close()

def plot_mc_stab():
  alpha, pot_e1, var_pot_e1, kin_e, var_kin_e, r12, Accepted =\
  np.loadtxt("results/mc_stab_1", skiprows=2,unpack=True) 
  alpha, pot_e2, var_pot_e2, kin_e, var_kin_e, r12, Accepted =\
  np.loadtxt("results/mc_stab_2", skiprows=2,unpack=True) 
  alpha, pot_e3, var_pot_e3, kin_e, var_kin_e, r12, Accepted =\
  np.loadtxt("results/mc_stab_3", skiprows=2,unpack=True) 
  alpha, pot_e4, var_pot_e4, kin_e, var_kin_e, r12, Accepted =\
  np.loadtxt("results/mc_stab_4", skiprows=2,unpack=True) 
  alpha, pot_e5, var_pot_e5, kin_e, var_kin_e, r12, Accepted =\
  np.loadtxt("results/mc_stab_5", skiprows=2,unpack=True) 
  alpha, pot_e6, var_pot_e6, kin_e, var_kin_e, r12, Accepted =\
  np.loadtxt("results/mc_stab_6", skiprows=2,unpack=True)

  plt.figure("Plot")
  plt.grid(1)

  plt.plot(alpha,pot_e1,'#76c9bf',marker='o',linestyle='dashed',
   markersize='3',label="$N_{mc}=10^1$", alpha=0.8)
  plt.plot(alpha,pot_e2,'#d9d942',marker='o',linestyle='dashed',
   markersize='3',label="$N_{mc}=10^2$", alpha=0.8)
  plt.plot(alpha,pot_e3,'#f07790',marker='o',linestyle='dashed',
   markersize='3',label="$N_{mc}=10^3$", alpha=0.8)
  plt.plot(alpha,pot_e4,'#d976a4',marker='o',linestyle='dashed',
   markersize='3',label="$N_{mc}=10^4$", alpha=0.8)
  plt.plot(alpha,pot_e5,'#f8ad49',marker='o',linestyle='dashed',
   markersize='3',label="$N_{mc}=10^5$", alpha=0.8)
  plt.plot(alpha,pot_e6,'#236964',marker='o',linestyle='dashed',
   markersize='3',label="$N_{mc}=10^6$", alpha=0.6)

  plt.xlabel("$\\alpha$",fontsize="14")
  plt.ylabel("$\langle E_L \\rangle$",fontsize="14")
  plt.legend()
  plt.savefig("results/mc_stab.pdf")
  plt.close()

def plot_mc_stab2():
  alpha, pot_e1, var_pot_e1, kin_e, var_kin_e, r12, Accepted =\
  np.loadtxt("results/mc_stab_1", skiprows=2,unpack=True) 
  alpha, pot_e2, var_pot_e2, kin_e, var_kin_e, r12, Accepted =\
  np.loadtxt("results/mc_stab_2", skiprows=2,unpack=True) 
  alpha, pot_e3, var_pot_e3, kin_e, var_kin_e, r12, Accepted =\
  np.loadtxt("results/mc_stab_3", skiprows=2,unpack=True) 
  alpha, pot_e4, var_pot_e4, kin_e, var_kin_e, r12, Accepted =\
  np.loadtxt("results/mc_stab_4", skiprows=2,unpack=True) 
  alpha, pot_e5, var_pot_e5, kin_e, var_kin_e, r12, Accepted =\
  np.loadtxt("results/mc_stab_5", skiprows=2,unpack=True) 
  alpha, pot_e6, var_pot_e6, kin_e, var_kin_e, r12, Accepted =\
  np.loadtxt("results/mc_stab_6", skiprows=2,unpack=True)

  plt.figure("Plot")
  plt.grid(1)

  plt.plot(alpha,var_pot_e1,'#76c9bf',marker='o',linestyle='dashed',
   markersize='3',label="$N_{mc}=10^1$", alpha=0.8)
  plt.plot(alpha,var_pot_e2,'#d9d942',marker='o',linestyle='dashed',
   markersize='3',label="$N_{mc}=10^2$", alpha=0.8)
  plt.plot(alpha,var_pot_e3,'#f07790',marker='o',linestyle='dashed',
   markersize='3',label="$N_{mc}=10^3$", alpha=0.8)
  plt.plot(alpha,var_pot_e4,'#d976a4',marker='o',linestyle='dashed',
   markersize='3',label="$N_{mc}=10^4$", alpha=0.8)
  plt.plot(alpha,var_pot_e5,'#f8ad49',marker='o',linestyle='dashed',
   markersize='3',label="$N_{mc}=10^5$", alpha=0.8)
  plt.plot(alpha,var_pot_e6,'#236964',marker='o',linestyle='dashed',
   markersize='3',label="$N_{mc}=10^6$", alpha=0.6)

  plt.xlabel("$\\alpha$",fontsize="14")
  plt.ylabel("$\langle E_L \\rangle$",fontsize="14")
  plt.legend()
  plt.savefig("results/mc_stab_var.pdf")
  plt.close()


def plot_test(filename):
  var, e, var_pot_e, kin_e, var_kin_e, r12, Accepted =\
  np.loadtxt("results/%s"%(filename), skiprows=2,unpack=True) 
  plt.figure("Plot")
  plt.grid(1)
  plt.plot(var,e,'#E57375',marker='o',linestyle='dashed', markersize='3')
  plt.xlabel("Variating",fontsize="14")
  plt.ylabel("$E$",fontsize="14")
  plt.savefig("results/test.pdf")
  plt.close()

file= "running_omega_w_jastrow"
file2= "running_omega_wo_jastrow"
file3="running_omega_wo_repuls"
identifier_1="_w_jastrow"
identifier_2="_wo_jastrow"
identifier_3="_wo_repuls"

test = "5d_run_1"

plot_test(test)

plot_r12_vs_omega(file,file2,identifier_1,identifier_2)

plot_r12_vs_omega_diff(file,file2,identifier_1,identifier_2)

plot_var_vs_omega(file,file2)

plot_energy_vs_omega(file,file2)


plot_kin_vs_pot(file,file3)


plot_energy_vs_alpha("1")
plot_energy_vs_alpha_var("1")

plot_mc_stab()
plot_mc_stab2()