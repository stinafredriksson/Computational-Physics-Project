#################################################################
## PROJECT 1
#################################################################

import math
import numpy as np
import matplotlib.pyplot as plt

#################################################################
## FUNCTIONS

def potential(U0, r, rmax):
    if r <= rmax:
        return U0
    else:
        return 0   
    
def r_min(b, U0, E, rmax):
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


def theta_b_less(b, rmax):
    return math.pi - 2*math.asin(b/rmax)

def theta_b_greater(b, rmax, U0, E):
    if 0 <= b < rmax*math.sqrt(1-U0/E):
        return 2*math.asin(b/rmax/(1-U0/E)**0.5) - 2*math.asin(b/rmax) 
    elif b<rmax:
        return theta_b_less(b,rmax)
    elif b>=rmax:
        return 0
    else:
        raise ValueError("b cannot be below zero")

def integral_1(p,b,rmin):
    if p==0:
        if b == 0:
            return 0
        return 2**0.5*(1/b)**1.5
    else:
        return 1/(p**2+b)**2*(1-b**2/(p**2+b)**2)**-0.5*2*p
    
def integral_2(p,b,rmin):
    if p==0:
        if b == 0:
            return 0
        return 2**0.5/b**1.5*(1-U0/E)**0.25
    else:
        return 1/(p**2+rmin)**2/(1-b**2/(p**2+rmin)**2-U0/E)**0.5*2*p
    

def plot_greater():

    bs = np.linspace(0,rmax,200)

    plt.title("Deflection angle over impact parameter")
    plt.plot(bs,[theta_b_less(bi,rmax) for bi in bs],c="C0",linestyle="--",label=r"$E\leq U_0$")
    plt.plot(bs, [theta_b_greater(bi,rmax, U0, E) for bi in bs],c="C1",linestyle="--",label=r"$E>U_0$")
    plt.grid(linestyle="--")
    plt.xlabel("b")
    plt.xticks([0,rmax*math.sqrt(1-U0/E),rmax],["0",r"$r_\text{max}\sqrt{1-\frac{U_0}{E}}$",r"$r_\text{max}$"])
    plt.ylabel(r"$\Theta(b)$")
    plt.legend()
    plt.yticks([0,math.pi/2,math.pi],[0,r"$\pi/2$",r"$\pi$"])
    plt.subplots_adjust(bottom=0.145)
    plt.show()


def plot_square():

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


def plot_difference():

    bs = np.linspace(0,rmax,100)

    N = 40

    numerical = []
    theory = []
    error = []

    for i, bi in enumerate(bs):

        print(f"{round(i/len(bs)*100)}",end="\r")

        if E>U0:
            rmin = r_min(bi, U0, E, rmax)
            theta_num_2=2*bi*boole(integral_2,0,(rmax-rmin)**0.5,N,bi,rmin)
        else:
            rmin = rmax
            theta_num_2 = 0
            

        theta_num_1=2*bi*boole(integral_1,0,(rmax-bi)**0.5,N,bi,rmin)
        
        theta_num=np.real(theta_num_1-theta_num_2)

        if E>U0:
            theta_analytic_values=theta_b_greater(bi,rmax,U0,E)
        else:
            theta_analytic_values=theta_b_less(bi,rmax)

        error.append(theta_num-theta_analytic_values)
        numerical.append(theta_num)
        theory.append(theta_analytic_values)

    plt.plot(bs,numerical,color="C0",linestyle="--",label="Numerical")
    plt.plot(bs,theory,color="C1",linestyle="dotted",label="Analytical")
    plt.legend()
    # plt.plot(bs,error)
    plt.show()

#################################################################
## GLOBALS

U0 = 2
E = 4
rmax = 10
# b = 0.5
# b = [rmax/14*i for i in range(15)]
# rmin = r_min(b, U0, E, rmax)
# b = [rmax/14*i for i in range(15)] # between 0 and rmax

#################################################################
## MAIN

def main():

    # plot_greater()
    # plot_square()
    plot_difference()

    b = 9
    
    # if E>U0:
    #     rmin = r_min(b, U0, E, rmax)
    # else:
    #     rmin = rmax

    # N = 40
    # theta_num_1=2*b*boole(integral_1,0,(rmax-b)**0.5,N,b,rmin)
    # theta_num_2=2*b*boole(integral_2,0,(rmax-rmin)**0.5,N,b,rmin)
    # theta_num=theta_num_1-theta_num_2
    # if E>U0:
    #     theta_analytic_values=theta_b_greater(b,rmax,U0,E)
    # else:
    #     theta_analytic_values=theta_b_less(b,rmax)

    
    # print(f"Numerical solution: {np.real(theta_num)}")
    # print(f"Analytical solution: {theta_analytic_values}")

#################################################################
## RUN CODE

if __name__ == '__main__':
    main()