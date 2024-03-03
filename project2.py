#################################################################
## PROJECT 2
#################################################################
## IMPORTS

import math

#################################################################
## FUNCTIONS

rmax = 1

def abs_db_dtheta(b):
    return rmax/2*math.sqrt(1-b**2/rmax**2)

def cross_section(b,theta):
    return b/math.sin(theta)*abs_db_dtheta(b)

def Lennard_Jones(r, a):
    return 4*V0*((a/r)**12-(a/r)**6)


## GLOBAL VARIABLES ##
V0 = 1

def main():

    pass

if __name__ == "__main__":
    main()