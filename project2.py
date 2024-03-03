#################################################################
## PROJECT 2
#################################################################
## IMPORTS

import math
import numpy as np
import matplotlib.pyplot as plt

#################################################################
## FUNCTIONS

a = 3
rmax = 3*a
V0 = 1

def abs_db_dtheta(b):
    return rmax/2*math.sqrt(1-b**2/rmax**2)

def cross_section(b,theta):
    return b/math.sin(theta)*abs_db_dtheta(b)

def potential(r):
    if r <= rmax:
        return V0*((a/r)**12-(a/r)**6)

def plot_potential():
    
    rs = np.linspace(0,rmax,200)

    plt.plot(rs,[potential(r) for r in rs])
    plt.show()
    

def main():
    plot_potential()


if __name__ == "__main__":
    main()