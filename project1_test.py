#################################################################
## PROJECT 1
#################################################################

import math
import numpy as np
from functions import boole
    
def r_min(b, U0, E, rmax):
    if E > U0:
        return b/math.sqrt(1-U0/E)
    else:
        return rmax

def theta_b_less(b, rmax):
    return math.pi - 2*math.asin(b/rmax)

def theta_b_greater(b, rmax, U0, E):
    return 2*math.asin(b/rmax/(1-U0/E)**0.5) - 2*math.asin(b/rmax) 
    
def integral_1(p):
    if p == 0:
        return 1.41421*(1/b)**1.5
    else:
        return 1/(p**2+b)**2*(1-b**2/(p**2+b)**2)**-0.5*2*p
    
def integral_2(p):
    if p == 0:
        return 0
    else:
        return 1/(p**2+rmin)**2*(1-b**2/(p**2+rmin)**2-U0/E)**-0.5*2*p

## GLOBAL VARIABLES ##

U0 = 1
E = 2
rmax = 5
b = 0.5
# b = [rmax/14*i for i in range(15)]
rmin = r_min(b, U0, E, rmax)
# b = [rmax/14*i for i in range(15)] # between 0 and rmax

def main():
    N = 40
    theta_num_1 = 2*b*boole(integral_1, 0, (rmax-b)**0.5, N)
    theta_num_2 = 2*b*boole(integral_2, 0, (rmax-rmin)**0.5, N)
    theta_num = theta_num_1 - theta_num_2
    if E> U0:
        theta_analytic_values = theta_b_greater(b, rmax, U0, E)
    else:
        theta_analytic_values = theta_b_less(b, rmax)
    
    print(f"Numerical solution: {theta_num}")
    print(f"Analytical solution: {theta_analytic_values}")



if __name__ == '__main__':
    main()