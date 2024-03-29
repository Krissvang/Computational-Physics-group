import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def crit_temp():
  temp, energy40, heatcap40, mag, absmag, sucept, time = \
    np.loadtxt("results/L=40.txt", skiprows=2,unpack=True)
  temp, energy60, heatcap60, mag, absmag, sucept, time = \
     np.loadtxt("results/L=60.txt", skiprows=2,unpack=True)
  temp, energy80, heatcap80, mag, absmag, sucept, time = \
    np.loadtxt("results/L=80.txt", skiprows=2,unpack=True)
  temp, energy100, heatcap100, mag, absmag, sucept, time = \
    np.loadtxt("results/L=100.txt", skiprows=2,unpack=True)
  x = np.array([1/100,1/80,1/60,1/40])
  index_max40=np.argmax(heatcap40)
  index_max60=np.argmax(heatcap60)
  index_max80=np.argmax(heatcap80)
  index_max100=np.argmax(heatcap100)
  y=np.array([temp[index_max100],temp[index_max80],
  temp[index_max60],temp[index_max40]])
  #slope, crit_temp, r_value, p_value, std_err = stats.linregress(x, y)

  n = len(x)
  D = sum(x**2)-(1./n)*(sum(x))**2
  E = sum(x*y)-(1./n)*sum(x)*sum(y)
  F = sum(y**2)-(1./n)*(sum(y))**2
  xm = np.mean(x)
  ym = np.mean(y)
  m = E/D
  c = ym-m*xm
  dm = np.sqrt(1./(n-2)*(D*F-E**2)/D**2)
  dc = np.sqrt(1./(n-2)*(D/n+xm**2)*(D*F-E**2)/D**2)
  crit_temp_exact=2/np.log(1+np.sqrt(2))
  slope = m; crit_temp = c; std_err = dc
  print("Calculated critical temperature = ", 
    crit_temp, "STD ERROR ", std_err)
  print("Difference between exact and calulated critical temp: ",
   abs(crit_temp_exact-crit_temp))
  plt.figure("Plot")
  plt.xlabel("$1/L$")
  plt.ylabel("T(1/L)")
  plt.plot(x,y, 'ro',label="Calculatet critical temp")
  plt.plot(x,crit_temp+slope*x,'b',label="Linreg of the points")
  plt.legend()
  plt.savefig("results/critical_temp_graph.pdf")
  plt.close()
crit_temp()
