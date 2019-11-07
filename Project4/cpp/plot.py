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

def plot_mcs_vs_configs(temperature):
  mcs, energy, heatcap, mag, absmag, sucept, time, configs = np.loadtxt("results/most_likely_state_T=%s.txt"%(temperature), skiprows=2,unpack=True) 
  plt.figure("Plot")
  plt.plot(mcs,configs, '#E57375')
  plt.xlabel("Number of Monte Carlo Cycles")
  plt.ylabel("Accepted configurations")
  plt.savefig("results/mcs_vs_configs_T=%s.png"%(temperature))
  plt.close()


def plot_energy(L,DT,MC):
  temp, energy, heatcap, mag, absmag, sucept, time = np.loadtxt("results/L=%s.txt"%(L), skiprows=2,unpack=True) 
  plt.figure("Plot")
  plt.plot(temp,energy,'ro')
  plt.xlabel("Temperature")
  plt.ylabel("Energy")
  plt.savefig("results/energy_L=%s_dT=%s_mc=%s.png"%(L,DT,MC))
  plt.close()

def plot_heatcap(L,DT,MC):
  temp, energy, heatcap, mag, absmag, sucept, time = np.loadtxt("results/L=%s.txt"%(L), skiprows=2,unpack=True) 
  plt.figure("Plot")
  plt.plot(temp,heatcap,'ro')
  plt.xlabel("Temperature")
  plt.ylabel("Heatcapacity")
  plt.savefig("results/heatcap_L=%s_dT=%s_mc=%s.png"%(L,DT,MC))
  plt.close()


def plot_absmag(L,DT,MC):
  temp, energy, heatcap, mag, absmag, sucept, time = np.loadtxt("results/L=%s.txt"%(L), skiprows=2,unpack=True) 
  plt.figure("Plot")
  plt.plot(temp,heatcap,'ro')
  plt.xlabel("Temperature")
  plt.ylabel("$\langle|M|\\rangle$")
  plt.savefig("results/absmag_L=%s_dT=%s_mc=%s.png"%(L,DT,MC))
  plt.close()

def plot_sucept(L,DT,MC):
  temp, energy, heatcap, mag, absmag, sucept, time = np.loadtxt("results/L=%s.txt"%(L), skiprows=2,unpack=True) 
  plt.figure("Plot")
  plt.plot(temp,sucept,'ro')
  plt.xlabel("Temperature")
  plt.ylabel("Heatcapacity ($\chi$)")
  plt.savefig("results/sucept_L=%s_dT=%s_mc=%s.png"%(L,DT,MC))
  plt.close()

#plot_mcs_vs_mag("1")
#plot_mcs_vs_energy("1")
#plot_mcs_vs_configs("")
plot_energy(60,0.005,7)
plot_heatcap(60,0.005,7)
plot_absmag(60,0.005,7)
plot_sucept(60,0.005,7)