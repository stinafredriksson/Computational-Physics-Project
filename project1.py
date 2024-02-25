#################################################################
## PROJECT 1
#################################################################

import math
import numpy as np


def potential(U0, r, rmax):
    if r <= rmax:
        return U0
    else:
        return 0   
    
def rmin(b,V,E,rmax):
    if E>V:
        return b/math.sqrt(1-V/E)
    else:
        return rmax

def analytical(E, U0):
    def theta_b_greater(b, rmax):
        return math.pi - 2*math.asin(b/rmax)
    
    def theta_b_less(b, rmax, U0, E):
        return 2*math.asin(b/rmax/(1-U0/E)) - 2*math.asin(b/rmax) 
    
    if E < U0:
        return theta_b_less
    else:
        return theta_b_greater
    
def f_E_less(p, b):
    return 1/(p**2+b)**2*(1-b**2/(p**2+b)**2)**-0.5*2*p


def f_E_greater_1(p, rmin, b):
    return 1/(p**2+rmin)**2(1-b**2/(p**2+rmin)**2)**-0.5*2*p


def f_E_greater_2(p, rmin, b, U0, E):
    return 1/(p**2+rmin)**2*(1-b**2/(p**2+rmin)**2 - U0/E)**-0.5*2*p


def main():
    

    pass

if __name__ == '__main__':
    main()