#################################################################
## PROJECT 1
#################################################################

import math
import time
import numpy as np
import matplotlib.pyplot as plt

#################################################################
## FUNCTIONS
    
def r_min(b: float, U0: float, E: float, rmax: float) -> float:
    """calculates rmin"""
    if E > U0:
        return b/math.sqrt(1-U0/E)
    else:
        return rmax
    

def boole(f, a1: float, a2:float, N:int, b, rmin) -> float:
    """
        f: integrated function
        a: lower bound
        b: upper bound
        N: divisions
    """

    assert (not N % 4) and N>0, "N has to exist and be divisible by 4"

    h = (a2-a1)/N

    return 2*h/45*(7*(f(a1,b,rmin)+f(a2,b,rmin)) + 32*sum([f(a1+k*h,b,rmin) for k in range(1,N,2)])+\
                   12*sum([f(a1+k*h,b,rmin) for k in range(2,N,4)])+\
                   14*sum([f(a1+k*h,b,rmin) for k in range(4,N,4)]))


def theta_b_less(b: float, rmax: float, U0: float, E: float) -> float:
    """analytical solution for E<=U0"""
    if U0 < 0:
        return theta_b_greater(b,rmax,U0,E)
    return math.pi - 2*math.asin(b/rmax)


def theta_b_greater(b: float, rmax: float, U0: float, E: float) -> float:
    """analytical solution for E>U0"""

    # particle enters the potential
    if 0 <= b < rmax*math.sqrt(1-U0/E):
        return 2*math.asin(b/rmax/(1-U0/E)**0.5) - 2*math.asin(b/rmax)
    # particle reflects from the potential
    elif b<rmax:
        return theta_b_less(b,rmax,U0,E)
    # particle misses the potential
    elif b>=rmax:
        return 0
    else:
        raise ValueError("b cannot be below zero")


def integral_1(p: float, b: float, rmin: float) -> float:
    """numerical setup for the first integral"""
    if p==0:
        if b == 0:
            return 0
        return 2**0.5*(1/b)**1.5
    else:
        return 1/(p**2+b)**2*(1-b**2/(p**2+b)**2)**-0.5*2*p
    

def integral_2(p: float, b: float, rmin: float) -> float:
    """numerical setup for the second integral"""
    if p==0:
        if b == 0:
            return 0
        return 2**0.5/b**1.5*(1-U0/E)**0.25
    else:
        return 1/(p**2+rmin)**2/(1-b**2/(p**2+rmin)**2-U0/E)**0.5*2*p
    

def integrate(b: float, N: float) -> float:
    """integrates the two integrals using Boole integration"""
    if E>U0:
        rmin = r_min(b,U0,E,rmax)
    else:
        rmin = rmax
    
    if b == 0 and E<=U0:
        theta_num_1 = math.pi
    else:
        theta_num_1=2*b*boole(integral_1,0,(rmax-b)**0.5,N,b,rmin)

    theta_num_2=2*b*boole(integral_2,0,(rmax-rmin)**0.5,N,b,rmin)
    
    return theta_num_1-np.real(theta_num_2)


def plot_analytical(U0: float, E:float, rmax:float) -> None:

    if 0 < E < U0:
        raise ValueError("Cannot use this function for $E<U_0$") 

    bs = np.linspace(0,rmax,200)

    plt.title("Deflection angle over impact parameter")
    plt.plot(bs,[theta_b_less(bi,rmax,U0,E) for bi in bs],c="C0",linestyle="--",label=r"$E\leq U_0$")
    plt.plot(bs,[theta_b_greater(bi,rmax,U0,E) for bi in bs],c="C1",linestyle="--",label=r"$E>U_0$")
    plt.grid(linestyle="--")
    plt.xlabel("b")

    if U0 > 0:
        plt.xticks([0,rmax*math.sqrt(1-U0/E),rmax],["0",r"$r_\text{max}\sqrt{1-\frac{U_0}{E}}$",r"$r_\text{max}$"])
        plt.yticks([0,math.pi/2,math.pi],[0,r"$\pi/2$",r"$\pi$"])
    else:
        plt.xticks([0,rmax/2,rmax],[0,r"$r_\text{max}/2$",r"$r_\text{max}$"])

    plt.ylabel(r"$\Theta(b)$")
    plt.legend()
    plt.subplots_adjust(bottom=0.145)
    plt.show()


def plot_square() -> None:

    def __square(r):
        if r<=rmax:
            return U0
        else:
            return 0

    rs = np.linspace(0,rmax+0.2*rmax,2000)

    plt.title("Square potential")
    plt.plot(rs,[__square(ri) for ri in rs])
    plt.xlabel("r")
    plt.ylabel(r"$V_\text{square}$")
    plt.yticks([0,U0/2,U0],["0",r"$U_0/2$",r"$U_0$"])
    plt.xticks([0,rmax/2,rmax],["0",r"$r_\text{max}/2$",r"$r_\text{max}$"])
    plt.grid(linestyle="--")
    plt.show()


def plot_difference() -> None:

    bs = np.linspace(0,rmax,2000)

    N = 160*4

    numerical = []
    theory = []
    error = []
    imag = []

    times = [time.time(),time.time()] # setup for time-left estimation

    for i, bnp in enumerate(bs):

        bi = bnp.item() # needed conversion numpy.float64 -> native float

        ## error bar
        ####
        times.pop(0)
        times.append(time.time())
        print(f"\r[{'#'*round(i/(len(bs)-1)*20):.<20}] {round(i/(len(bs)-1)*100):02}% |{round((times[1]-times[0])*(len(bs)-1-i)):02}s|",end="\r")
        ####

        theta_num = integrate(bi, N)

        if E>abs(U0) or U0<0:
            theta_analytic_values=theta_b_greater(bi,rmax,U0,E)
        else:
            theta_analytic_values=theta_b_less(bi,rmax,U0,E)

        error.append(abs(theta_num-theta_analytic_values))
        numerical.append(theta_num)
        theory.append(theta_analytic_values)
    
    plt.figure(1)
    plt.plot(bs,numerical,color="C0",linestyle="--",label="Numerical")
    plt.plot(bs,theory,color="C1",linestyle="dotted",label="Analytical")
    plt.xlabel("b")

    if E>abs(U0) and U0>0:
        title = r"$E>|U_0|$"
        plt.yticks([0,math.pi/4,math.pi/2],[0,r"$\pi/4$",r"$\pi/2$"])
        plt.xticks([0,rmax*math.sqrt(1-U0/E),rmax],[0,r"$r_\text{max}\sqrt{1-\frac{U_0}{E}}$",r"$r_\text{max}$"])
        plt.subplots_adjust(bottom=0.145)
    else:
        title = r"$E\leq |U_0|$"
        plt.xticks([0,rmax/2,rmax],[0,r"$r_\text{max}/2$",r"$r_\text{max}$"])
        if U0>0:
            plt.yticks([0,math.pi/2,math.pi],[0,r"$\pi/2$",r"$\pi$"])

    if U0 < 0:
        title2 = r"$U_0<0$"
    else:
        title2 = r"$U_0>0$"

    plt.title(f"{title}, {title2}")

    plt.grid(linestyle="--")

    # plt.plot(bs,imag)
    plt.legend()

    plt.figure(2)
    plt.title(f"{title}, {title2}")
    plt.plot(bs,error)
    plt.ylabel("|Error|")

    if E>abs(U0):
        title = r"$E>|U_0|$"
        if U0 > 0:
            plt.xticks([0,rmax*math.sqrt(1-U0/E),rmax],[0,r"$r_\text{max}\sqrt{1-\frac{U_0}{E}}$",r"$r_\text{max}$"])
            plt.subplots_adjust(bottom=0.145)
        else:
            plt.xticks([0,rmax/2,rmax],[0,r"$r_\text{max}/2$",r"$r_\text{max}$"])
    else:
        title = r"$E\leq |U_0|$"
        plt.xticks([0,rmax/2,rmax],[0,r"$r_\text{max}/2$",r"$r_\text{max}$"])
    
    plt.yscale("log")
    plt.grid(linestyle="--")
    plt.show()


def plot_ordo() -> None:

    b = 8

    Ns = [k*4 for k in range(10,101)]

    error = []

    conv = []

    for i, N in enumerate(Ns):

        print(f"{round(i/len(Ns)*100)}",end="\r")

        theta_num = integrate(b, N)

        if E>U0:
            theta_analytic_values=theta_b_greater(b,rmax,U0,E)
        else:
            theta_analytic_values=theta_b_less(b,rmax,U0,E)

        error.append(abs(theta_num-theta_analytic_values))

        if i != 0:
            conv.append(math.log(abs(error[i])/abs(error[i-1]))/math.log(Ns[i-1]/Ns[i]))
    
    print(f"q(N<80) = {round(np.mean(conv[:19]),2)}")

    # plt.plot(Ns[1:],conv)
    # plt.ylabel("order of convergence q")
            
    plt.plot(Ns,error)
    plt.xlabel("N")
    plt.yscale("log")
    plt.ylabel("log(|Error|)")
    plt.grid(linestyle="--")
    plt.show()

#################################################################
## GLOBALS

U0 = -1
E = 1
rmax = 10

#################################################################
## MAIN

def main():

    # plot_analytical(U0,E,rmax)
    plot_square()
    plot_difference()
    # plot_ordo()

    # b = 0
    # N = 40
    
    # theta_num=integrate(b,N)

#################################################################
## RUN CODE

if __name__ == '__main__':
    main()