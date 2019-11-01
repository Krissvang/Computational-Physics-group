import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def plot_mcs_vs_energy(temperature):
  mcs, energy, heatcap, mag, absmag, sucept, time = np.loadtxt("results/most_likely_state_T=%s.txt"%(temperature), skiprows=2,unpack=True) 
  plt.figure("Plot")
  plt.plot(mcs,energy,'#E57375')
  plt.xlabel("Number of Monte Carlo Cycles")
  plt.ylabel("Average energy")
  plt.savefig("results/mcs_vs_energy_T=%s.png"%(temperature))
  plt.close()

def plot_mcs_vs_mag(temperature):
  mcs, energy, heatcap, mag, absmag, sucept, time = np.loadtxt("results/most_likely_state_T=%s.txt"%(temperature), skiprows=2,unpack=True) 
  plt.figure("Plot")
  plt.plot(mcs,absmag, '#E57375')
  plt.xlabel("Number of Monte Carlo Cycles")
  plt.ylabel("Average absolute value of magnetization")
  plt.savefig("results/mcs_vs_magnetization_T=%s.png"%(temperature))
  plt.close()


plot_mcs_vs_mag("1")
plot_mcs_vs_energy("1")

