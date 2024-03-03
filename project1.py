#################################################################
## PROJECT 1
#################################################################

import math
import numpy as np
from functions import boole
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


def f_E_less(p):
    if p ==0:
        return 0
    else:       
        return 1/(p**2+b)**2*(1-b**2/(p**2+b)**2)**-0.5*2*p

def f_E_greater_1(p):
    return 1/(p**2+rmin)**2*(1-b**2/(p**2+rmin)**2)**-0.5*2*p

def f_E_greater_2(p):
    if p ==0:
        return 0
    else:
        return 1/(p**2+rmin)**2*(1-b**2/(p**2+rmin)**2 - U0/E)**-0.5*2*p
    

def plot_greater():

    bs = np.linspace(0,rmax,200)


    plt.title("Deflection angle over impact parameter")
    plt.plot(bs,[theta_b_less(bi,rmax) for bi in bs],c="C0",linestyle="--",label=r"$E\leq V$")
    plt.plot(bs, [theta_b_greater(bi,rmax, U0, E) for bi in bs],c="C1",linestyle="--",label=r"E>V")
    plt.grid(linestyle="--")
    plt.xlabel("b")
    plt.xticks([0,rmax*math.sqrt(1-U0/E),rmax],["0",r"$r_\text{max}\sqrt{1-\frac{U_0}{E}}$",r"$r_\text{max}$"])
    plt.ylabel(r"$\Theta(b)$")
    plt.legend()
    plt.yticks([0,math.pi/2,math.pi],[0,r"$\pi/2$",r"$\pi$"])
    plt.subplots_adjust(bottom=0.145)
    plt.show()

#################################################################
## GLOBALS

U0 = 2
E = 4
rmax = 5
b = 0.5
# b = [rmax/14*i for i in range(15)]
rmin = r_min(b, U0, E, rmax)
# b = [rmax/14*i for i in range(15)] # between 0 and rmax

#################################################################
## MAIN

def main():

    plot_greater()

    # N = 20
    # if E < U0:
    #     theta_num = 2*b*boole(f_E_less, 0.00001, (rmax-b)**0.5, N)
    #     theta_analytic_values = theta_b_less(b, rmax)
    # else:
    #     theta_num = 2*b* boole(f_E_greater_2, 0, (rmax-rmin)**0.5, N)
    #     theta_analytic_values = theta_b_greater(b, rmax, U0, E)
    
    # print(f"Numerical solution: {theta_num}")
    # print(f"Analytical solution: {theta_analytic_values}")
    # print(rmin)


#################################################################
## RUN CODE

if __name__ == '__main__':
    main()