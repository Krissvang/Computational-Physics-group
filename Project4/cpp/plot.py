import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def plot_mcs_vs_energy(temperature):
  mcs, energy, heatcap, mag, absmag, sucept, time, configs = np.loadtxt("results/most_likely_state_T=%s.txt"%(temperature), skiprows=2,unpack=True) 
  plt.figure("Plot")
  plt.plot(mcs,energy/20**2,'#E57375')
  plt.xlabel("Number of Monte Carlo Cycles")
  plt.ylabel("Average energy per spin")
  plt.savefig("results/mcs_vs_energy_T=%s.pdf"%(temperature))
  plt.close()

def plot_mcs_vs_mag(temperature):
  mcs, energy, heatcap, mag, absmag, sucept, time, configs = np.loadtxt("results/most_likely_state_T=%s.txt"%(temperature), skiprows=2,unpack=True) 
  plt.figure("Plot")
  plt.plot(mcs,absmag/20**2, '#E57375')
  plt.xlabel("Number of Monte Carlo Cycles")
  plt.ylabel("Average absolute value of magnetization per spin")
  plt.savefig("results/mcs_vs_magnetization_T=%s.pdf"%(temperature))
  plt.close()

def plot_mcs_vs_configs(temperature):
  mcs, energy, heatcap, mag, absmag, sucept, time, configs = np.loadtxt("results/most_likely_state_configs_T=%s.txt"%(temperature), skiprows=2,unpack=True) 
  plt.figure("Plot")
  plt.plot(mcs,configs/20**2, '#E57375')
  plt.xlabel("Number of Monte Carlo Cycles")
  plt.ylabel("Accepted configurations per spin")
  plt.savefig("results/mcs_vs_configs_T=%s.pdf"%(temperature))
  plt.close()


def plot_energy():
  temp, energy40, heatcap, mag, absmag, sucept, time = np.loadtxt("results/L=40.txt", skiprows=2,unpack=True) 
  temp, energy60, heatcap, mag, absmag, sucept, time = np.loadtxt("results/L=60.txt", skiprows=2,unpack=True)
  temp, energy80, heatcap, mag, absmag, sucept, time = np.loadtxt("results/L=80.txt", skiprows=2,unpack=True)
  temp, energy100, heatcap, mag, absmag, sucept, time = np.loadtxt("results/L=100.txt", skiprows=2,unpack=True)
  plt.figure("Plot")
  plt.plot(temp,energy40/40**2,c='#4D4D4D',marker='o',linestyle="", alpha=.7, markersize=4, label="L=40")
  plt.plot(temp,energy60/60**2,c='#8CDCDA',marker='o', linestyle="", alpha=.7, markersize=4, label="L=60")
  plt.plot(temp,energy80/80**2,c='#B1D877',marker='o', linestyle="", alpha=.7, markersize=4, label="L=80")
  plt.plot(temp,energy100/100**2,c='#F16A70',marker='o', linestyle="", alpha=.7, markersize=4, label="L=100")
  plt.legend()
  plt.xlabel("Temperature")
  plt.ylabel("Energy per spin")
  plt.savefig("results/energy.pdf")
  plt.close()

def plot_heatcap():
  temp, energy40, heatcap40, mag, absmag, sucept, time = np.loadtxt("results/L=40.txt", skiprows=2,unpack=True) 
  temp, energy60, heatcap60, mag, absmag, sucept, time = np.loadtxt("results/L=60.txt", skiprows=2,unpack=True)
  temp, energy80, heatcap80, mag, absmag, sucept, time = np.loadtxt("results/L=80.txt", skiprows=2,unpack=True)
  temp, energy100, heatcap100, mag, absmag, sucept, time = np.loadtxt("results/L=100.txt", skiprows=2,unpack=True)
  plt.figure("Plot")
  plt.plot(temp,heatcap40/40**2,c='#4D4D4D',marker='o',linestyle="", alpha=.7, markersize=4, label="L=40")
  plt.plot(temp,heatcap60/60**2,c='#8CDCDA',marker='o', linestyle="", alpha=.7, markersize=4, label="L=60")
  plt.plot(temp,heatcap80/80**2,c='#B1D877',marker='o', linestyle="", alpha=.7, markersize=4, label="L=80")
  plt.plot(temp,heatcap100/100**2,c='#F16A70',marker='o', linestyle="", alpha=.7, markersize=4, label="L=100")
  plt.legend()
  plt.xlabel("Temperature")
  plt.ylabel("Heatcapacity per spin")
  plt.savefig("results/heatcap.pdf")
  plt.close()


def plot_absmag():
  temp, energy40, heatcap, mag, absmag40, sucept, time = np.loadtxt("results/L=40.txt", skiprows=2,unpack=True) 
  temp, energy60, heatcap, mag, absmag60, sucept, time = np.loadtxt("results/L=60.txt", skiprows=2,unpack=True)
  temp, energy80, heatcap, mag, absmag80, sucept, time = np.loadtxt("results/L=80.txt", skiprows=2,unpack=True)
  temp, energy100, heatcap, mag, absmag100, sucept, time = np.loadtxt("results/L=100.txt", skiprows=2,unpack=True)
  plt.figure("Plot")
  plt.plot(temp,absmag40/40**2,c='#4D4D4D',marker='o',linestyle="", alpha=.7, markersize=4, label="L=40")
  plt.plot(temp,absmag60/60**2,c='#8CDCDA',marker='o', linestyle="", alpha=.7, markersize=4, label="L=60")
  plt.plot(temp,absmag80/80**2,c='#B1D877',marker='o', linestyle="", alpha=.7, markersize=4, label="L=80")
  plt.plot(temp,absmag100/100**2,c='#F16A70',marker='o', linestyle="", alpha=.7, markersize=4, label="L=100")
  plt.legend()
  plt.xlabel("Temperature")
  plt.ylabel("$\langle|M|\\rangle$ per spin")
  plt.savefig("results/absmag.pdf")
  plt.close()

def plot_sucept():
  temp, energy40, heatcap, mag, absmag, sucept40, time = np.loadtxt("results/L=40.txt", skiprows=2,unpack=True) 
  temp, energy60, heatcap, mag, absmag, sucept60, time = np.loadtxt("results/L=60.txt", skiprows=2,unpack=True)
  temp, energy80, heatcap, mag, absmag, sucept80, time = np.loadtxt("results/L=80.txt", skiprows=2,unpack=True)
  temp, energy100, heatcap, mag, absmag, sucept100, time = np.loadtxt("results/L=100.txt", skiprows=2,unpack=True)
  plt.figure("Plot")
  plt.plot(temp,sucept40/40**2,c='#4D4D4D',marker='o',linestyle="", alpha=.7, markersize=4, label="L=40")
  plt.plot(temp,sucept60/60**2,c='#8CDCDA',marker='o', linestyle="", alpha=.7, markersize=4, label="L=60")
  plt.plot(temp,sucept80/80**2,c='#B1D877',marker='o', linestyle="", alpha=.7, markersize=4, label="L=80")
  plt.plot(temp,sucept100/100**2,c='#F16A70',marker='o', linestyle="", alpha=.7, markersize=4, label="L=100")
  plt.legend()
  plt.xlabel("Temperature")
  plt.ylabel("Heatcapacity ($\chi$) per spin")
  plt.savefig("results/sucept.pdf")
  plt.close()


T="1"
T1="2,4"
plot_mcs_vs_mag(T)
plot_mcs_vs_energy(T)
plot_mcs_vs_configs(T)
plot_mcs_vs_mag(T1)
plot_mcs_vs_energy(T1)
plot_mcs_vs_configs(T1)
plot_energy()
plot_heatcap()
plot_absmag()
plot_sucept()