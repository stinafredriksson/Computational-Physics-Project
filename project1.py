#################################################################
## PROJECT 1
#################################################################

import math
import numpy as np
from functions import boole


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
    return 2*math.asin(b/rmax/(1-U0/E)**0.5) - 2*math.asin(b/rmax) 
    
def f_E_less(p):       
    return 1/(p**2+b)**2*(1-b**2/(p**2+b)**2)**-0.5*2*p

def f_E_greater_1(p):
    return 1/(p**2+rmin)**2*(1-b**2/(p**2+rmin)**2)**-0.5*2*p

def f_E_greater_2(p):
    # p = (r-rmin)**0.5
    return 1/(p**2+rmin)**2*(1-b**2/(p**2+rmin)**2 - U0/E)**-0.5*2*p


## GLOBAL VARIABLES ##

U0 = 1
E = 2
rmax = 5
b = 2
# b = [rmax/14*i for i in range(15)]
rmin = r_min(b, U0, E, rmax)
# b = [rmax/14*i for i in range(15)] # between 0 and rmax

def main():
    N = 20
    # print(rmin)
    if E < U0:
        # print('hello')
        theta_num = 2*b*boole(f_E_less, 0.00001, (rmax-b)**0.5, N)
        theta_analytic_values = theta_b_less(b, rmax)
    else:
        theta_num = 2*b*(boole(f_E_greater_1, (b-rmin)**0.5, (rmax-rmin)**0.5, N) - boole(f_E_greater_2, 0, (rmax-rmin)**0.5, N))
        theta_analytic_values = theta_b_greater(b, rmax, U0, E)
    
    print(f"Numerical solution: {theta_num}")
    print(f"Analytical solution: {theta_analytic_values}")

if __name__ == '__main__':
    main()