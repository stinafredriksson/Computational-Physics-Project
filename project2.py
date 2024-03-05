#################################################################
## PROJECT 2
#################################################################
## IMPORTS

import math
import numpy as np
import matplotlib.pyplot as plt

#################################################################
## FUNCTIONS

a = 5
rmax = 3*a
V0 = 1

def abs_db_dtheta(b):
    return rmax/2*math.sqrt(1-b**2/rmax**2)

def cross_section(b,theta):
    return b/math.sin(theta)*abs_db_dtheta(b)

def potential(r):
    if r <= rmax:
        return 4*V0*((a/r)**12-(a/r)**6)

def plot_potential():
    
    rs1 = np.linspace(0.2*a,a,200)
    rs2 = np.linspace(a,3*a,200)

    plt.plot(rs1,[potential(r) for r in rs1], c="C0",label=r"Positive")
    plt.plot(rs2,[potential(r) for r in rs2], c="C1",label=r"Negative")
    # plt.plot(rs,[4*V0*(a/r)**12 for r in rs], c="C1",linestyle="dotted",label=r"$4V_0\left(\frac{a}{r}\right)^{12}$")
    # plt.plot(rs,[a**6/(a**6-r**6) for r in rs])
    plt.legend()
    plt.grid(linestyle="--")
    plt.xlabel("r")
    plt.ylabel("V")
    plt.title("Lennard Jones potential")
    plt.yticks(color = "w")
    plt.xticks([0,a,3*a], [0,r"$a$",r"$r_\text{max}$"])
    plt.show()
    

def main():
    plot_potential()


if __name__ == "__main__":
    main()