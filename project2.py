#################################################################
## PROJECT 2
#################################################################
## IMPORTS

import math
import numpy as np
import matplotlib.pyplot as plt

from project1 import integral_1, boole
from functions import search

#################################################################
## GLOBALS

a = 5
rmax = 3*a
V0 = 2
E = 1

#################################################################
## FUNCTIONS

def r_min(b,V):

    def __integrand(r):
        return 1-b**2/r**2-4*V*((a/r)**12-(a/r)**6)

    return search(__integrand,a-0.4,1e-4)
    

def abs_db_dtheta(b):
    return rmax/2*math.sqrt(1-b**2/rmax**2)


def cross_section(b,theta):
    return b/math.sin(theta)*abs_db_dtheta(b)


def potential(r):
    if r <= rmax:
        return 4*V0*((a/r)**12-(a/r)**6)
    else:
        return 0
    

def integral_2(p: float, b: float, rmin: float) -> float:
    """numerical setup for the second integral"""
    if p==0:
        if b == 0:
            return 0
        return 2**0.5/b**1.5*(1-potential(p**2+rmin)/E)**0.25
    else:
        return 1/(p**2+rmin)**2/(1-b**2/(p**2+rmin)**2-potential(p**2+rmin)/E)**0.5*2*p


def integrate(b: float, N: int, V: float) -> float:
    """integrates the two integrals using Boole integration"""

    if b == 0:
        return math.pi

    rmin = r_min(b, V)

    theta_num_1=2*b*boole(integral_1,0,(rmax-b)**0.5,N,b,rmin)

    theta_num_2=2*b*boole(integral_2,0,(rmax-rmin)**0.5,N,b,rmin)
    
    return theta_num_1-np.real(theta_num_2)


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
    

def plot_rmin():

    bs = np.linspace(0,rmax,2000)
    Vs = [0.5,1,2,5,10]
    zeros = [[r_min(b,V) for b in bs] for V in Vs]


    for i,zero in enumerate(zeros):
        plt.plot(bs,zero,c=f"C{i}",label=f"V0={Vs[i]}E")
    plt.grid(linestyle="--")
    plt.xlabel("b")
    plt.ylabel(r"$r_{min}(b,V_0)$")
    plt.xticks([0,a,2*a,rmax],[4,r"$a$",r"$2a$",r"$r_{max}$"])
    plt.title(r"$r_{min}$ over b")
    plt.legend()
    plt.show()


def plot_integral():

    bs = np.linspace(0,rmax,200)
    Vs = [0.5,1,2,5,10]
    intes = [[],[],[],[],[]]

    for j,V in enumerate(Vs):
        for bi in bs:

            b = bi.item()

            intes[j].append(integrate(b,160*4,V))

    for i,inte in enumerate(intes):
        plt.plot(bs,inte,color=f"C{i}",label=f"V0 = {Vs[i]}E")
    plt.legend()
    plt.xlabel("b")
    plt.ylabel(r"$\Theta(b,V_0)$")
    plt.yticks([-2*math.pi,-math.pi,0,math.pi],[r"-2$\pi$",r"-$\pi$",0,r"$\pi$"])
    plt.xticks([0,a,2*a,3*a],[0,r"$a$",r"$2a$",r"$r_\text{max}$"])
    plt.grid(linestyle="--")
    plt.show()



def main():
    # plot_potential()
    # plot_rmin()
    plot_integral()

    # print(r_min(10,V0))



if __name__ == "__main__":
    main()