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

def plot_r12_vs_omega(filename):
  var, pot_e, var_pot_e, kin_e, var_kin_e,pe_wo,var_pe_wo,pe_w,var_pe_w, r12, Accepted =\
  np.loadtxt("results/%s"%(filename), skiprows=2,unpack=True) 
  plt.figure("Plot")
  plt.grid(1)
  plt.plot(var,r12,'#E57375',marker='o',linestyle='dashed', markersize='3')
  plt.xlabel("$\omega$",fontsize="14")
  plt.ylabel("$r_{12}$",fontsize="14")
  plt.savefig("results/r12_vs_omega.png")
  plt.close()

def plot_energy_vs_omega(filename):
  var, e, var_pot_e, kin_e, var_kin_e,pe_wo,var_pe_wo,pe_w,var_pe_w, r12, Accepted =\
  np.loadtxt("results/%s"%(filename), skiprows=2,unpack=True) 
  plt.figure("Plot")
  plt.grid(1)
  plt.plot(var,e,'#E57375',marker='o',linestyle='dashed', markersize='3')
  plt.xlabel("$\omega$",fontsize="14")
  plt.ylabel("$\langle E_L \\rangle$",fontsize="14")
  plt.savefig("results/energy_vs_omega.png")
  plt.close()

def plot_var_vs_omega(filename):
  var, e, var_e, kin_e, var_kin_e,pe_wo,var_pe_wo,pe_w,var_pe_w, r12, Accepted =\
  np.loadtxt("results/%s"%(filename), skiprows=2,unpack=True) 
  plt.figure("Plot")
  plt.grid(1)
  plt.plot(var,var_e,'#E57375',marker='o',linestyle='dashed', markersize='3')
  plt.xlabel("$\omega$",fontsize="14")
  plt.ylabel("$\sigma^2$",fontsize="14")
  plt.savefig("results/variance_vs_omega.png")
  plt.close()

def plot_test(filename):
  var, pot_e, var_pot_e, kin_e, var_kin_e,pe_wo,var_pe_wo,pe_w,var_pe_w, r12, Accepted =\
  np.loadtxt("results/%s"%(filename), skiprows=2,unpack=True) 
  plt.figure("Plot")
  plt.grid(1)
  plt.plot(var,r12,'#E57375',marker='o',linestyle='dashed', markersize='3')
  plt.xlabel("Variating",fontsize="14")
  plt.ylabel("$E$",fontsize="14")
  plt.savefig("results/test.png")
  plt.close()

file= "running_omega"
plot_test(file)
plot_r12_vs_omega(file)
plot_var_vs_omega(file)
plot_energy_vs_omega(file)
#plot_energy_vs_alpha(w)