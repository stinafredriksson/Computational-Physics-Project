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

def theta_b_less(b, rmax, U0, E):
    return 2*math.asin(b/rmax/(1-U0/E)) - 2*math.asin(b/rmax)    

def theta_b_greater(b, rmax):
    return math.pi - 2*math.asin(b/rmax)


def analytical(E, U0):

    if E < U0:
        return theta_b_less
    else:
        return theta_b_greater

    
