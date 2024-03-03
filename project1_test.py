#################################################################
## PROJECT 1
#################################################################

import math
import numpy as np
from functions import boole
import matplotlib.pyplot as plt
    
def r_min(b,U0,E,rmax):
    if E>U0:
        return b/math.sqrt(1-U0/E)
    else:
        return rmax

def theta_b_less(b,rmax):
    return math.pi-2*math.asin(b/rmax)

<<<<<<< HEAD
def theta_b_greater(b, rmax, U0, E):
    return 2*np.arcsin(b/rmax/(1-U0/E)**0.5) - 2*np.arcsin(b/rmax)
    
def integral_1(p):
    if p == 0:
=======
def theta_b_greater(b,rmax,U0,E):
    return 2*math.asin(b/rmax/(1-U0/E)**0.5) - 2*math.asin(b/rmax) 
    
def integral_1(p):
    if p==0:
>>>>>>> refs/remotes/origin/main
        return 2**0.5*(1/b)**1.5
    else:
        return 1/(p**2+b)**2*(1-b**2/(p**2+b)**2)**-0.5*2*p
    
def integral_2(p):
    if p==0:
        return 2**0.5/b**1.5*(1-U0/E)**0.25
    else:
        return 1/(p**2+rmin)**2*(1-b**2/(p**2+rmin)**2-U0/E)**-0.5*2*p
    
def check_large_b():

    def __arc(b):
        return 2*np.arcsin(b/rmax/(1-U0/E)**0.5)

    bs = np.linspace(0,rmax)
    # print(rmax*math.sqrt(1-U0/E))

    print(theta_b_greater(rmax*(1-U0/E)**0.5,rmax,U0,E)/math.pi)

    plt.plot(bs,theta_b_greater(bs,rmax,U0,E))
    plt.show()

## GLOBAL VARIABLES ##
U0 = 1
V0 = 1
E = 0.1
rmax = 5
# a
b = 0.5
rmin = r_min(b,U0,E,rmax)

def main():
<<<<<<< HEAD
    check_large_b()
    # N = 40
    # theta_num_1 = 2*b*boole(integral_1, 0, (rmax-b)**0.5, N)
    # theta_num_2 = 2*b*boole(integral_2, 0, (rmax-rmin)**0.5, N)
    # theta_num = theta_num_1 - theta_num_2
    # if E> U0:
    #     theta_analytic_values = theta_b_greater(b, rmax, U0, E)
    # else:
    #     theta_analytic_values = theta_b_less(b, rmax)
=======
    N = 40
    theta_num_1=2*b*boole(integral_1,0,(rmax-b)**0.5,N)
    theta_num_2=2*b*boole(integral_2,0,(rmax-rmin)**0.5,N)
    theta_num=theta_num_1-theta_num_2
    if E>U0:
        theta_analytic_values=theta_b_greater(b,rmax,U0,E)
    else:
        theta_analytic_values=theta_b_less(b,rmax)
>>>>>>> refs/remotes/origin/main
    
    # print(f"Numerical solution: {theta_num}")
    # print(f"Analytical solution: {theta_analytic_values}")



if __name__ == '__main__':
    main()