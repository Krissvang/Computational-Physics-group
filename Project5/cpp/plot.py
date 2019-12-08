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

def plot_test(filename):
  alpha, pot_e, var_pot_e, kin_e, var_kin_e, r12, Accepted =\
  np.loadtxt("results/%s"%(filename), skiprows=2,unpack=True) 
  plt.figure("Plot")
  plt.grid(1)
  plt.plot(alpha,var_pot_e,'#E57375',marker='o',linestyle='dashed', markersize='3')
  plt.xlabel("Variating",fontsize="14")
  plt.ylabel("$Variance$",fontsize="14")
  plt.savefig("results/test.png")
  plt.close()

file= "5d_run_12"
plot_test(file)
#plot_energy_vs_alpha(w)