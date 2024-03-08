#################################################################
## PROJECT 2
#################################################################
## IMPORTS

import math
import numpy as np
import matplotlib.pyplot as plt

from project1 import integral_1
from functions import search

#################################################################
## GLOBALS

a = 1
rmax = 3*a
V0 = 1


#################################################################
## FUNCTIONS

def boole(f, a1: float, a2:float, N:int, b, rmin, E) -> float:
    """
        f: integrated function
        a: lower bound
        b: upper bound
        N: divisions
    """

    assert (not N % 4) and N>0, "N has to exist and be divisible by 4"

    h = (a2-a1)/N

    if E:
        return 2*h/45*(7*(f(a1,b,rmin,E)+f(a2,b,rmin,E)) + 32*sum([f(a1+k*h,b,rmin,E) for k in range(1,N,2)])+\
                   12*sum([f(a1+k*h,b,rmin,E) for k in range(2,N,4)])+\
                   14*sum([f(a1+k*h,b,rmin,E) for k in range(4,N,4)]))
    else:
        return 2*h/45*(7*(f(a1,b,rmin)+f(a2,b,rmin)) + 32*sum([f(a1+k*h,b,rmin) for k in range(1,N,2)])+\
                   12*sum([f(a1+k*h,b,rmin) for k in range(2,N,4)])+\
                   14*sum([f(a1+k*h,b,rmin) for k in range(4,N,4)]))

def r_min(b,E):

    def __integrand(r):
        return 1-b**2/r**2-4*V0*((a/r)**12-(a/r)**6)/E

    return search(__integrand,a-0.4,1e-10)



def cross_section(b,Theta,db_dTheta):
    if b>rmax:
        return 0
    else:
        return b/math.sin(Theta)*abs(db_dTheta)

def differentiate(b, N, E):
    h=0.01
    return (integrate(b+h,N,E)-integrate(b-h,N,E))/(2*h)

def potential_E(r, E):
    if r <= rmax:
        return 4*V0*((a/r)**12-(a/r)**6)/E
    else:
        return 0
    

def integral_2(p: float, b: float, rmin: float, E) -> float:
    """numerical setup for the second integral"""
    if p==0:
        if b == 0:
            return 0
        return 2**0.5/b**1.5*(1-potential_E(p**2+rmin,E))**0.25
    else:
        return 1/(p**2+rmin)**2/(1-b**2/(p**2+rmin)**2-potential_E(p**2+rmin,E))**0.5*2*p


def integrate(b: float, N: int, E) -> float:
    """integrates the two integrals using Boole integration"""

    if b == 0:
        return math.pi

    rmin = r_min(b, E)

    theta_num_1=2*b*boole(integral_1,0,(rmax-b)**0.5,N,b,rmin,None)

    theta_num_2=2*b*boole(integral_2,0,(rmax-rmin)**0.5,N,b,rmin,E)
    
    return theta_num_1-np.real(theta_num_2)

def plot_cross_section():

    bs = np.linspace(0.1,rmax-0.1,200)
    Vs = [0.1,0.5,1,5,10,100]
    cs = [[],[],[],[],[],[]]
    N=160

    for j,V in enumerate(Vs):
        for bi in bs:

            b = bi.item()
            E=V*V0

            dTheta_db = differentiate(b,N,E)
            Theta = integrate(b,N,E)
            cs[j].append(abs(cross_section(b,Theta,dTheta_db**-1)))
    
    for i,csi in enumerate(cs):
        plt.plot(bs,csi,color=f"C{i}",label=f"E={Vs[i]}V0")
    plt.yscale("log")
    plt.grid(linestyle='--')
    plt.xlabel(r"$b$")
    plt.legend()
    plt.ylabel(r"$d\sigma/d\Omega$")
    plt.title(r"Cross-Section vs $b$")
    plt.show()


def plot_potential():
    E=1
    rs1 = np.linspace(0.2*a,a,200)
    rs2 = np.linspace(a,3*a,200)

    plt.plot(rs1,[potential_E(r,E) for r in rs1], c="C0",label=r"Positive")
    plt.plot(rs2,[potential_E(r,E) for r in rs2], c="C1",label=r"Negative")
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

    bs = np.linspace(0,rmax,200)
    Vs = [0.1,0.5,1,5,10,100]
    zeros = [[r_min(b,V*V0) for b in bs] for V in Vs] # E=V*V0


    for i,zero in enumerate(zeros):
        plt.plot(bs,zero,c=f"C{i}",label=f"E={Vs[i]}V0")
    plt.grid(linestyle="--")
    plt.xlabel(r"$b$")
    plt.ylabel(r"$r_{min}(b,E)$")
    plt.xticks([0,a,2*a,rmax],[4,r"$a$",r"$2a$",r"$r_{max}$"])
    plt.title(r"$r_{min}$ over $b$")
    plt.legend()
    plt.show()


def plot_integral():

    bs = np.linspace(0,rmax,200)
    Vs = [0.1,0.5,1,5,10,100]
    intes = [[],[],[],[],[],[]]

    for j,V in enumerate(Vs):
        for bi in bs:

            b = bi.item()
            E=V0*V

            intes[j].append(integrate(b,160*4,E))

    for i,inte in enumerate(intes):
        plt.plot(bs,inte,color=f"C{i}",label=f"E={Vs[i]}V0")
    plt.legend()
    plt.xlabel("b")
    plt.ylabel(r"$\Theta(b,E)$")
    plt.yticks([-2*math.pi,-math.pi,0,math.pi],[r"-2$\pi$",r"-$\pi$",0,r"$\pi$"])
    plt.xticks([0,a,2*a,3*a],[0,r"$a$",r"$2a$",r"$r_\text{max}$"])
    plt.grid(linestyle="--")
    plt.show()



def main():
    # plot_potential()
    # plot_rmin()
    # plot_integral()
    # print(r_min(10,V0))
    plot_cross_section()




if __name__ == "__main__":
    main()