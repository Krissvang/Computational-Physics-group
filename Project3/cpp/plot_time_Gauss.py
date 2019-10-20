import numpy as np
import matplotlib.pyplot as plt

x = np.log10(np.array([10,15,20,25,35,45]))
y = np.log10(np.array([0.0269326,0.224246,1.20471,4.27829,32.0729,143.811]))

plt.figure()
plt.plot(x,y,"x",label=r"$\sigma_N$")
p = np.polyfit(x,y,1)
plt.plot(x,p[0]*x+p[1],label="%gx+%g"%(p[0],p[1]))
plt.xlabel("log10(N)")
plt.ylabel("log10(t)")
plt.legend()
plt.savefig("img/plot_time_Gauss.pdf")
