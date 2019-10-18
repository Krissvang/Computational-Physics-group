import numpy as np
import matplotlib.pyplot as plt

x = np.arange(4,10)
y = np.log10(np.array([0.009,0.004,0.0011,0.0003,0.0001,0.00005]))

plt.figure()
plt.plot(x,y,"x",label=r"$\sigma_N$")
p = np.polyfit(x,y,1)
plt.plot(x,p[0]*x+p[1],label="%gx+%g"%(p[0],p[1]))
plt.xlabel("log10(N)")
plt.ylabel(r"log10($\sigma_N$)")
plt.legend()
plt.axis([3.8,9.2,1.1*min(y),0.9*max(y)])
plt.savefig("img/plot_std.pdf")
