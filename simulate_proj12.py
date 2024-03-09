#################################################################
## RK4 SIMULATION OF PROJECTS 1 AND 2
#################################################################
## IMPORTS

import math
import time
import numpy as np
import matplotlib.pyplot as plt

from project2 import potential
from project2 import integrate as integrate_LJ

#################################################################
## CLASSES

class Vector():
    def __init__(self,x: float = 0, y: float = 0) -> None:
        self.x = x
        self.y = y

    def __abs__(self) -> float:
        return math.sqrt(self.x**2+self.y**2)
    
    def __add__(self, other):
        return Vector(self.x+other.x, self.y+other.y)
    
    def __sub__(self, other):
        return Vector(self.x-other.x, self.y-other.y)
    
    def __mul__(self, other):
        return Vector(self.x*other, self.y*other)
    
    def __truediv__(self, other):
        return Vector(self.x/other, self.y/other)
    
    def norm(self) -> float:
        return Vector(self.x/abs(self),self.y/abs(self))
    
    def __str__(self) -> str:
        return f"({round(self.x,2)},{round(self.y,2)})"
    
    def __rmul__(self,other):
        return self.__mul__(other)
    
    def dot(self,other):
        return self.x*other.x + self.y*other.y

class Particle():
    def __init__(self, r: Vector = Vector(0,0), v: Vector = Vector(0,0)) -> None:
        self.r = r
        self.v = v

#################################################################
## FUNCTIONS
        
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


def a(r: Vector, rmax: float, U0: float, a_const = 5, LJ_flag: bool = False) -> float:
    if LJ_flag:
        if r <= 3*a_const:
            return  potential(r, U0, a_const)
        else: 
            return 0
    if r <= rmax:
        return U0/rmax
    else:
        return 0
    

def rk4(particle: Particle, h: float, rmax: float, E: float, U0: float, LJ_flag: bool = False) -> None:

    r = particle.r
    v = particle.v

    k_v1 = a(abs(r),rmax, U0, LJ_flag = LJ_flag)*r.norm()
    k_r1 = v

    k_v2 = a(abs(r+k_r1*h/2), rmax, U0, LJ_flag = LJ_flag)*r.norm()
    k_r2 = v + k_v1*h/2

    k_v3 = a(abs(r+k_r2*h/2), rmax, U0, LJ_flag = LJ_flag)*r.norm()
    k_r3 = v + k_v2*h/2

    k_v4 = a(abs(r+k_r3*h), rmax, U0, LJ_flag = LJ_flag)*r.norm()
    k_r4 = v + k_v3*h

    particle.r = r + h/6*(k_r1 + 2*k_r2 + 2*k_r3 + k_r4)
    particle.v = v + h/6*(k_v1 + 2*k_v2 + 2*k_v3 + k_v4)

    # if E<U0 and abs(particle.r) < rmax:
    #     angle_loc = math.acos(-particle.r.x/abs(particle.r))

    #     angle = math.pi-2*angle_loc


    #     vel_dir = Vector(math.cos(angle),math.sin(angle)+1)

    #     particle.r = rmax*particle.r.norm()
    #     particle.v = abs(particle.v)*vel_dir


def plot_path(xs: list[float], ys: list[float], rmax: float, b: float, energy_fraction: float, a_const: float = 5, LJ_flag: bool = False) -> None:

    def __find_i(xs,target,r_flag=0):

        i =  0
        if r_flag:
            while xs[i]>target:
                i+=1
                if i == len(xs)-1:
                    return -1
        else:
            while xs[i]<target:
                i+=1
                if i == len(xs)-1:
                    return -1
        return i

    fig, ax = plt.subplots()
    if LJ_flag:
        potential_inner = plt.Circle((0,0),a_const,color="grey",linestyle="--", fill = False,label="Potential inner")
        rmax = 3*a_const
    potential_c = plt.Circle((0,0),rmax,color="k",linestyle="--", fill = False,label="Potential")
    source = plt.Circle((0,0),rmax/20,color="k")

    i1 = __find_i(xs,-rmax*0.7)
    i2 = __find_i(xs,rmax*0.8)
    if i1 == -1 or i2 == -1:
        i1 = __find_i(xs,-rmax*1.5)
        i2 = __find_i(xs[i1+1:],-rmax*1.5,1) + i1+1

    if LJ_flag:
        comp = "Lennard-Jones potential"
    elif abs(energy_fraction) > 1:
        comp = r"$E<|U_0|$"
    elif energy_fraction == 1:
        comp = r"$E=|U_0|$"
    else:
        comp = r"$E>|U_0|$"
    if LJ_flag:
        if energy_fraction > 1:
            unit = r"$E>V_0$"
        elif energy_fraction == 1:
            unit = r"$E=V_0$"
        else:
            unit = r"$E<V_0$"
    elif energy_fraction < 0:
        unit = r"$U_0<0$"
    else:
        unit = r"$U_0>0$"

    plt.title(comp + ", " + unit)

    ax.grid(linestyle="--")
    ax.hlines([b],-2*rmax,2*rmax,colors="C0",linestyles="--",label="y=b")
    ax.add_patch(potential_c)
    if LJ_flag:
        ax.add_patch(potential_inner)
    ax.add_patch(source)
    ax.hlines([0],-2*rmax,2*rmax,colors="k",linestyles="--")
    ax.plot(xs,ys,label="Path",color="C1")
    ax.arrow(xs[i1],ys[i1],xs[i1+1]-xs[i1-1],ys[i1+1]-ys[i1-1],shape='full', color="C1",lw=0, length_includes_head=True, head_width=1)
    ax.arrow(xs[i2],ys[i2],xs[i2+1]-xs[i2-1],ys[i2+1]-ys[i2-1],shape='full', color="C1",lw=0, length_includes_head=True, head_width=1)
    ax.set_aspect("equal")
    ax.set_ylim(-2*rmax,2*rmax)
    ax.set_xlim(-2*rmax,2*rmax)
    ax.legend()
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    plt.show()


def plot_r(ts: list[float], rs: list[float], i_in: int, i_out: int, b: float) -> None:

    plt.title("E > |U_0|")
    plt.plot(ts, rs, label="r")
    plt.vlines([ts[i_in], ts[i_out]],min(rs),max(rs), colors="C1",linestyles="--",label="Inside potential")
    plt.hlines([b],min(ts),max(ts),colors="C2",label="b",linestyles="--")
    plt.legend(loc="upper center")
    plt.grid(linestyle="--")
    plt.xlabel("time")
    plt.ylabel("r")
    plt.show()


def simulate(rmax: float, b: float, U0: float, E: float, a_const: float = 5, LJ_flag: bool = False):

    t = 0
    h = 0.01
    TMAX = 60

    if LJ_flag:
        rmax = 3*a_const

    p = Particle(Vector(-2*rmax,b), Vector(math.sqrt(2*E),0))
    rs = [abs(p.r)]
    xs = [p.r.x]
    ys = [p.r.y]
    ts = [0]

    i_in = 0
    i_out = 0
    i = 0

    in_flag = False

    while t < TMAX:

        if abs(p.r) < 10 and not in_flag:
            i_in = i
            in_flag = True

        if in_flag and abs(p.r)>=10:
            i_out = i
            in_flag = False

        rk4(p,h,rmax,E,U0, LJ_flag = LJ_flag)

        rs.append(abs(p.r))
        xs.append(p.r.x)
        ys.append(p.r.y)

        t += h
        ts.append(t)

        i += 1

        # if abs(p.r) < 0.01:
        #     break

    return ts,rs,[xs,ys],[i_in,i_out],p.v


def compare(U0: float, E: float, rmax: float, LJ_flag: bool = False):

    ta = []
    tc = []
    bs = np.linspace(0,rmax,100)
    times = [time.time(),time.time()]

    for i,b in enumerate(bs):

        times.pop(0)
        times.append(time.time())

        print(f"\r[{'#'*round(i/(len(bs)-1)*20):.<20}] {round(i/(len(bs)-1)*100):02}% |{round((times[1]-times[0])*(len(bs)-1-i)):02}s|",end="\r")

        ts,rs,[xs,ys],[i_in,i_out], vel = simulate(rmax,b, U0, E, LJ_flag)

        if abs(U0/E) > 1:
            theta_anal = theta_b_less(b,rmax,U0,E)
        else:
            theta_anal = theta_b_greater(b,rmax,U0,E)

        # ts,rs,[xs2,ys2],[i_in2,i_out2] = simulate(rmax,b, 4, 2)

        theta_comp = math.acos(vel.x/abs(vel))

        ta.append(theta_anal)
        tc.append(theta_comp)

        # print(f"Theta_anal = {round(180/math.pi*theta_anal,2)}")
        # print(f"Theta_comp = {round(180/math.pi*theta_comp,2)}")

    # print(rmax*math.sqrt(1-U0/E))

    plt.plot(bs,ta,label="Analytical")
    plt.plot(bs,tc,label="Numerical")
    plt.xlabel("b")
    plt.ylabel(r"$\Theta(b)$")
    plt.grid(linestyle="--")
    plt.legend()
    plt.show()


def run_single(U0: float, E: float, rmax: float, b: float, LJ_flag: bool = False) -> None:

    ts,rs,[xs,ys],[i_in,i_out], vel = simulate(rmax,b, U0, E, LJ_flag= LJ_flag)

    v0 = math.sqrt(2*E)


    theta_sim = np.sign(vel.y)*math.acos(vel.x/abs(vel))

    if LJ_flag:

        theta_anal = integrate_LJ(b,160,U0,E)

    else:
        if E>U0:
            theta_anal = theta_b_greater(b,rmax,U0,E)
        else:
            theta_anal = theta_b_less(b,rmax,U0,E)

    print(f"Simulated: {round(180/math.pi*theta_sim,2)}")
    print(f"Theoretical: {round(180/math.pi*theta_anal,2)}")

    plot_path(xs,ys,rmax,b,U0/E,LJ_flag=LJ_flag)

#################################################################
## MAIN

def main():

    U0 = 10
    E = 1

    rmax = 10
    b = 4

    # compare(U0,E,rmax)

    run_single(U0,E,rmax,b,LJ_flag=True)


#################################################################
## RUN CODE

if __name__ == "__main__":
    main()
