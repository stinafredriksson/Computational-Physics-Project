import math
import numpy as np
import matplotlib.pyplot as plt


rmin = 1
E = 3
U0 = 2
b = 0.5

def f(p):
    return 2*p/(rmin)/np.sqrt((1-U0/E)*(2*p**2*rmin+rmin**2)-b**2)


p = np.linspace(0,10,200)

# print(f(0))

plt.plot(p,f(p))
plt.show()