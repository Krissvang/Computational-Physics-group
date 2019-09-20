import numpy as np
import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.figure()
plt.plot([1, 2])
plt.ylabel('some numbers')
plt.savefig("test.png")