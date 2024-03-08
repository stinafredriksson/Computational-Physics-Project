#################################################################
## IMPORTS

import math
from typing import Sequence

#################################################################
## FUNCTIONS

## DERIVATIVE
####

def euler(f:Sequence[float], h:float) -> float:
    """
        f: [f0, f1] forward
        f: [f-1, f0] backward
        h: step length
    """ 
    return (f[1]-f[0])/h

def three_point(f: Sequence[float], h:float) -> float:
    """
        f: [f-1, f1]
        h: step length
    """
    return (f[1]-f[0])/(2*h)

def five_point(f: Sequence[float], h:float) -> float:
    """
        f: [f-2, f-1, f1, f2]
        h: step length
    """
    return 1/(12*h)*(f[0]-8*f[1]+8*f[2]-f[3])
####


## INTEGRATION
####

def trapezoid(f, a: float, b: float, N:int) -> float:
    """
        f: integrated function
        a: lower bound
        b: upper bound
        N: divisions
    """
    h = (b-a)/N

    return h/2*(f(a) + 2*sum([f(a+k*h) for k in range(1,N)]) + f(b))


def simpson(f, a: float, b:float, N:int) -> float:
    """
        f: integrated function
        a: lower bound
        b: upper bound
        N: divisions
    """
    assert (not N % 2) and N>0, "N has to exist and be even"

    h = (b-a)/N

    return h/3*(f(a)+ 4*sum([f(a+k*h) for k in range(1,N,2)]) +\
                2*sum([f(a+k*h) for k in range(2,N,2)]) + f(b))


def boole(f, a: float, b:float, N:int) -> float:
    """
        f: integrated function
        a: lower bound
        b: upper bound
        N: divisions
    """
    assert (not N % 4) and N>0, "N has to exist and be divisible by 4"

    h = (b-a)/N

    return 2*h/45*(7*(f(a)+f(b)) + 32*sum([f(a+k*h) for k in range(1,N,2)])+\
                   12*sum([f(a+k*h) for k in range(2,N,4)])+\
                   14*sum([f(a+k*h) for k in range(4,N,4)]))
####


## ROOTS
####

def search(f, trial: float, tolerance: float, N: int = 100, positive: bool = True, debug: bool = False) -> float:
    """
        f: function to find root for
        trial: initial guess and output guess
        tolerance: approximate wanted closeness of answer
        N: maximum allowed iterations
        positive: True if you are looking for a positive root
        debug: True if you want information on number of iterations
    """

    if positive:
        sgn = 1
    else:
        sgn = -1


    def __sign(x):
        return int(x/abs(x)) if x else 0

    h = 0.1
    
    for i in range(N):
        
        if h < tolerance:
            if debug:
                print(f"{i+1} iterations")
            return trial
        j = 0
        while __sign(f(trial)) == __sign(f(trial+sgn*h)) and j < N:
            trial = trial + sgn*h
        trial = trial - sgn*h
        h = h/2
    else:
        if debug:
            print("max iterations")
        return trial


def newton_raphson(f, fprim, guess: float, tolerance: float, N: int = 100, debug: bool = False) -> float:
    """
        f: function to find root for
        fprim: derivative of f
        guess: initial guess
        tolerance: approximate wanted closeness of answer
        N: maximum allowed iterations
        debug: True if you want information on number of iterations
    """

    x = guess

    for i in range(N):
        
        x_last = x

        x = x - f(x)/fprim(x)

        if abs(x-x_last) < tolerance:
            if debug:
                print(f"{i+1} iterations")
            return x
    else:
        if debug:
            print("max iterations")
        return x


def secant(f, guess1: float, guess2: float, tolerance: float, N: int = 100, debug: bool = False) -> float:
    """
        f: function to find root for
        guess1: initial worse guess
        guess2: initial better guess
        tolerance: approximate wanted closeness of answer
        N: maximum allowed iterations
        debug: True if you want information on number of iterations
    """

    x = guess2
    x_last = guess1

    for i in range(N):

        x_next = x-f(x)*(x-x_last)/(f(x)-f(x_last))

        if abs(x_next-x) < tolerance:
            if debug:
                print(f"{i+1} iterations")
            return x_next

        x_last = x
        x = x_next
    else:
        if debug:
            print("max iterations")
        return x_next
####


## ODE SOLVERS
####

def euler_ODE(f, a: float, b: float, y_initial: float, N: int) -> float:
    """
        f: dy/dx ODE to solve
        a: initial x
        b: x to solve
        y_initial: IC, y(a)
        N: number of domain partitions
    """

    h = (b-a)/N

    y = y_initial

    for k in range(N):
        y_next = y + h*f(a+k*h, y)
        y = y_next
    
    return y_next


def taylor_ODE(f, fdx, fdy, a: float, b:float, y_initial: float, N: int) -> float:
    """
        f: dy/dx ODE to solve
        fdx: (dy/dx)dx
        fdy: (dy/dx)dy
        a: initial x
        b: x to solve
        y_initial: IC, y(a)
        N: number of domain partitions
    """

    h = (b-a)/N

    y = y_initial

    for k in range(N):

        y_next = y + h*f(a+k*h,y) + (h**2)/2*(fdx(a+k*h,y)+f(a+k*h,y)*fdy(a+k*h,y))

        y = y_next
    
    return y_next


def implicit_ODE(g, a: float, b: float, y_initial: float, N: int) -> float:
    """
        g: dy/dx = g(x)y ODE to solve
        a: initial x
        b: x to solve
        y_initial: IC, y(a)
        N: number of domain partitions
    """

    h = (b-a)/N

    y = y_initial

    for k in range(N):
        
        y_next = ((1+h/2*g(a+k*h))/(1-h/2*g(a+(k+1)*h)))*y

        y = y_next
    
    return y_next
####


## RUNGE-KUTTA
####

def __k1(f, x, y, h):
    return h*f(x,y)

def __k2(f, x, y, h, k1):
    return h*f(x+h/2,y+k1/2)

def RK2(f, h: float, a: float, b: float, y_initial: float, temporal: bool = False):

    x = a

    y = y_initial

    if temporal:
        x_lst = [x]
        y_lst = [y]

    while x <= b:
        
        k = __k1(f, x, y, h)

        y = y + __k2(f, x, y, h, k)

        x += h

        if temporal:
            x_lst.append(x)
            y_lst.append(y)
        
    if temporal:
        return y, [x_lst, y_lst]
    else:
        return y

def RK3(f, h: float, a: float, b: float, y_initial: float, temporal: bool = False):

    x = a

    y = y_initial

    if temporal:
        x_lst = [x]
        y_lst = [y]

    while x <= b:
        
        k1 = __k1(f, x, y, h)
        k2 = __k2(f, x, y, h, k1)
        k3 = h*f(x+h,y-k1+2*k2)
        k3_2 = h*f(x+h/2, y+k2/2)

        print(k3, k3_2)

        y = y + 1/6*(k1+4*k2+k3)

        x += h

        if temporal:
            x_lst.append(x)
            y_lst.append(y)
        
    if temporal:
        return y, [x_lst, y_lst]
    else:
        return y


def RK4(f, h: float, a: float, b: float, y_initial: float, temporal: bool = False):

    x = a

    y = y_initial

    if temporal:
        x_lst = [x]
        y_lst = [y]

    while x <= b:
        
        k1 = h*f(x,y)
        k2 = h*f(x+h/2, y+k1/2)
        k3 = h*f(x+h/2, y+k2/2)
        k4 = h*f(x+h, y+k3)

        y = y + 1/6*(k1+2*k2+2*k3+k4)

        x += h

        if temporal:
            x_lst.append(x)
            y_lst.append(y)
        
    if temporal:
        return y, [x_lst, y_lst]
    else:
        return y

####

#################################################################
## MAIN

def main():

    def df(x,y):
        return -x*y

    def dfdx(x,y):
        return -y
    
    def dfdy(x,y):
        return -x
    
    def g(x):
        return -x
    
    res = euler_ODE(df, 0, 3, 1, 100)

    res2 = taylor_ODE(df, dfdx, dfdy, 0, 3, 1, 100)

    res3 = implicit_ODE(g, 0, 3, 1, 100)

    print(res - math.exp(-9/2))
    print(res2 - math.exp(-9/2))
    print(res3 - math.exp(-9/2))

#################################################################
## RUN CODE

if __name__ == "__main__":
    main()