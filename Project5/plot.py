import matplotlib
matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt

#Have to use this for unix subsystem



omega= "1"

def plot_energy_vs_alpha(w):
  alpha, E_wo_c, V_wo_c, E_w_c, V_w_c, Accepted =\
  np.loadtxt("results/energy_vs_alpha_w=%s"%(w), skiprows=1,unpack=True) 
  plt.figure("Plot")
  plt.plot(alpha,E_wo_c,'#E57375',marker='o',linestyle='dashed', markersize='3')
  plt.xlabel("Alpha",fontsize="14")
  plt.ylabel("$E_0$",fontsize="14")
  plt.savefig("results/energy_vs_alpha_w=%s.png"%(w))
  plt.close()

plot_energy_vs_alpha(omega)