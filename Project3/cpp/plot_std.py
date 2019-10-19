import numpy as np
import matplotlib.pyplot as plt

x = np.arange(4,10)
y = np.log10(np.array([0.00943,0.00359,0.00106,0.00033,0.00011,0.00003]))

plt.figure()
plt.plot(x,y,"x",label=r"$\sigma_N$")
p = np.polyfit(x,y,1)
plt.plot(x,p[0]*x+p[1],label="%gx+%g"%(p[0],p[1]))
plt.xlabel("log10(N)")
plt.ylabel(r"log10($\sigma_N$)")
plt.legend()
plt.axis([3.8,9.2,1.1*min(y),0.9*max(y)])
plt.savefig("img/plot_std.pdf")
