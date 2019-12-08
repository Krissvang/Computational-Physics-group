import matplotlib
matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt

#Have to use this for unix subsystem



w= "1"

def plot_energy_vs_alpha(omega):
  alpha, pot_e, var_pot_e, kin_e, var_kin_e, r12, Accepted =\
  np.loadtxt("results/energy_vs_alpha_w=%s"%(omega), skiprows=2,unpack=True) 
  plt.figure("Plot")
  plt.grid(1)
  plt.plot(alpha,pot_e,'#E57375',marker='o',linestyle='dashed', markersize='3')
  plt.xlabel("Alpha",fontsize="14")
  plt.ylabel("$E_0$",fontsize="14")
  plt.savefig("results/energy_vs_alpha_w=%s.png"%(omega))
  plt.close()

def plot_r12_vs_omega(filename,identifier):
  var, pot_e, var_pot_e, kin_e, var_kin_e,pe_wo,var_pe_wo,pe_w,var_pe_w, r12, Accepted =\
  np.loadtxt("results/%s"%(filename), skiprows=2,unpack=True) 
  plt.figure("Plot")
  plt.grid(1)
  plt.plot(var,r12,'#E57375',marker='o',linestyle='dashed', markersize='3')
  plt.xlabel("$\omega$",fontsize="14")
  plt.ylabel("$r_{12}$",fontsize="14")
  plt.savefig("results/r12_vs_omega%s.png"%(identifier))
  plt.close()

def plot_energy_vs_omega(filename,identifier):
  var, e, var_pot_e, kin_e, var_kin_e,pe_wo,var_pe_wo,pe_w,var_pe_w, r12, Accepted =\
  np.loadtxt("results/%s"%(filename), skiprows=2,unpack=True) 
  plt.figure("Plot")
  plt.grid(1)
  plt.plot(var,e,'#E57375',marker='o',linestyle='dashed', markersize='3')
  plt.xlabel("$\omega$",fontsize="14")
  plt.ylabel("$\langle E_L \\rangle$",fontsize="14")
  plt.savefig("results/energy_vs_omega%s.png"%(identifier))
  plt.close()

def plot_var_vs_omega(filename,identifier):
  var, e, var_e, kin_e, var_kin_e,pe_wo,var_pe_wo,pe_w,var_pe_w, r12, Accepted =\
  np.loadtxt("results/%s"%(filename), skiprows=2,unpack=True) 
  plt.figure("Plot")
  plt.grid(1)
  plt.plot(var,var_e,'#E57375',marker='o',linestyle='dashed', markersize='3')
  plt.xlabel("$\omega$",fontsize="14")
  plt.ylabel("$\sigma^2$",fontsize="14")
  plt.savefig("results/variance_vs_omega%s.png"%(identifier))
  plt.close()


def plot_kin_vs_pot_wo(filename,identifier):
  var, e, var_pot_e, kin_e, var_kin_e,pe_wo,var_pe_wo,pe_w,var_pe_w, r12, Accepted =\
  np.loadtxt("results/%s"%(filename), skiprows=2,unpack=True) 
  plt.figure("Plot")
  plt.grid(1)
  plt.plot(var,kin_e/pe_wo,'#E57375',marker='o',linestyle='dashed', markersize='3')
  plt.xlabel("$\omega$",fontsize="14")
  plt.ylabel("$\\frac{\langle T \\rangle}{\langle V_{wo} \\rangle}$",fontsize="14")
  plt.savefig("results/ke_vs_pe_unperturbed%s.png"%(identifier))
  plt.close()

def plot_kin_vs_pot_w(filename,identifier):
  var, e, var_pot_e, kin_e, var_kin_e,pe_wo,var_pe_wo,pe_w,var_pe_w, r12, Accepted =\
  np.loadtxt("results/%s"%(filename), skiprows=2,unpack=True) 
  plt.figure("Plot")
  plt.grid(1)
  plt.plot(var,kin_e/pe_w,'#E57375',marker='o',linestyle='dashed', markersize='3')
  plt.xlabel("$\omega$",fontsize="14")
  plt.ylabel("$\\frac{\langle T \\rangle}{\langle V_{w} \\rangle}$",fontsize="14")
  plt.savefig("results/ke_vs_pe_perturbed%s.png"%(identifier))
  plt.close()



def plot_test(filename):
  var, e, var_pot_e, kin_e, var_kin_e,pe_wo,var_pe_wo,pe_w,var_pe_w, r12, Accepted =\
  np.loadtxt("results/%s"%(filename), skiprows=2,unpack=True) 
  plt.figure("Plot")
  plt.grid(1)
  plt.plot(var,r12,'#E57375',marker='o',linestyle='dashed', markersize='3')
  plt.xlabel("Variating",fontsize="14")
  plt.ylabel("$E$",fontsize="14")
  plt.savefig("results/test.png")
  plt.close()








file= "running_omega_w_jastrow"
file2= "running_omega_wo_jastrow"
file3="running_omega_wo_repuls"
identifier_1="_w_jastrow"
identifier_2="_wo_jastrow"
identifier_3="_wo_repuls"
plot_test(file)
plot_r12_vs_omega(file,identifier_1)
plot_var_vs_omega(file,identifier_1)
plot_energy_vs_omega(file,identifier_1)

plot_kin_vs_pot_w(file,identifier_1)

plot_kin_vs_pot_wo(file3,identifier_3)

plot_r12_vs_omega(file2,identifier_2)
plot_var_vs_omega(file2,identifier_2)
plot_energy_vs_omega(file2,identifier_2)
#plot_energy_vs_alpha(w)