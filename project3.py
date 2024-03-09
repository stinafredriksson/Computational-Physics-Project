#################################################################
## PROJECT 3
#################################################################
## IMPORTS

import numpy as np
import time
import math
import matplotlib.pyplot as plt
from scipy.constants import Boltzmann as k

#################################################################
## CLASSES

class System():
    """This class handles the low-temperature behavior of a 
    lattice material"""

    def __init__(self,N: int, J: float = k, B: float = 0, T: float = 2) -> None:
        """
        parameters:
            N: lattice size [int]
            J: coupling strength [float] (usully on the order of kB)
            B: external magnetic field strngth [float]
            T: system temperature [float]
        """

        self.N = N
        self.J = J
        self.B = B
        self.T = T
        self.steps = 0
        self.accepted = 0
        # initializes the spin lattice with random spins
        self.lattice = np.random.choice([-1,1,1,1],size=(N,N))


    def magnetization(self) -> float:
        return np.sum(self.lattice)/self.N**2
    

    def hamiltonian(self) -> float:
        H = 0
        for i in range(self.N):
            for j in range(self.N):
                H += -self.J*self.lattice[i][j]*(self.lattice[(i+1)%self.N][j]+self.lattice[i-1][j]+self.lattice[i][(j+1)%self.N]+self.lattice[i][j-1])
                H += -self.B*np.sum(self.lattice)
        return H


    def calc_dE(self,i: int, j: int) -> float:

        return 2*self.J*self.lattice[i][j]*(self.lattice[(i+1)%self.N][j]+self.lattice[i-1][j]+self.lattice[i][(j+1)%self.N]+self.lattice[i][j-1])

    def step(self): ## METROPOLIS 
        i, j = np.random.randint(self.N), np.random.randint(self.N)
        dE = self.calc_dE(i,j)

        # input("")
        # print(np.exp(-dE/k/self.T))

        if dE < 0 or np.random.rand() < np.exp(-dE/k/self.T):
            self.accepted += 1
            self.lattice[i][j] = - self.lattice[i][j]
        
        self.steps += 1

    def equilibrialization(self, Nsweeps: int = 10) -> None:
        """sweeps for reaching thermal equilibrium
        the sweeps are not super intensive for N<=32"""

        for _ in range(Nsweeps):

            # do sweep
            for i in range(self.N):
                for j in range(self.N):
                    dE = self.calc_dE(i,j)
                    if dE < 0 or np.random.rand() < np.exp(-dE/k/self.T):
                        self.lattice[i][j] = - self.lattice[i][j]


    def show_state(self) -> None:
        plt.imshow(self.lattice,cmap="Greys")
        ax = plt.gca()
        # ax.set_xticks(np.arange(0, self.N, 1))
        # ax.set_yticks(np.arange(0, self.N, 1))

        ax.set_xticks([])
        ax.set_yticks([])

        # Labels for major ticks
        # ax.set_xticklabels(np.arange(1, self.N+1, 1))
        # ax.set_yticklabels(np.arange(1, self.N+1, 1))

        # Minor ticks
        ax.set_xticks(np.arange(-.5, self.N, 1), minor=True)
        ax.set_yticks(np.arange(-.5, self.N, 1), minor=True)

        # Gridlines based on minor ticks
        ax.grid(which='minor', color='grey', linestyle='dotted', linewidth=0.5)

        # Remove minor ticks
        ax.tick_params(which='minor', bottom=False, left=False)

        plt.xlabel("i")
        plt.ylabel("j")
        plt.title(f"Square Ising model with N={self.N}")
        plt.show()

#################################################################
## FUNCTIONS

def estimate_order(N,J,T,steps,number=10):

    order = []

    times = [time.time()]*2

    for j in range(number):
        system = System(N,J,T)
        for i in range(steps):

            times.pop(0)
            times.append(time.time())
            print(f"\r[{'#'*round((i+j*steps)/(steps*number -1)*20):.<20}] {round((i+j*steps)/(steps*number-1)*100):02}% |{round((times[1]-times[0])*(steps*number-1-i-j*steps)):04}s|",end="\r")
            system.step()

        order.append(abs(system.magnetization()))

    print("\nDone")

    return np.mean(order)    


def order_over_T(N,J,steps):

    Ts = [1,5,10,20,40,50,80,100]

    orders = [estimate_order(N,J,T,steps) for T in Ts]

    plt.plot(Ts,orders)
    plt.show()


def LLsizes(Magnetization_flag: None, Susceptibility_flag: None, Specific_flag:None):
    Ts = np.linspace(1,4,4)

    sizes = [4,8,16,32]
    magtot = [[] for _ in range(len(sizes))]
    avg = 100
    times = [time.time(), time.time()]
    suscept =  [[] for _ in range(len(sizes))]
    energy_tot = [[] for _ in range(len(sizes))]
    specific_tot = [[] for _ in range(len(sizes))]


    for t,Ti in enumerate(Ts):

        # averages over iterations
        if Magnetization_flag:
            mags = [[] for _ in range(len(sizes))]
        if Susceptibility_flag:
            mags_sq = [[] for _ in range(len(sizes))]
        if Specific_flag:
            energys = [[] for _ in range(len(sizes))]
            specifics = [[] for _ in range(len(sizes))]

        for j in range(avg):
            times.pop(0)
            times.append(time.time())
            print(f"\r[{'#'*round((j+t*avg)/(avg*len(Ts) -1)*20):.<20}] {round((j+t*avg)/(avg*len(Ts)-1)*100):03}% |{round((times[1]-times[0])*(avg*len(Ts)-1-j-t*avg)):04}s|",end="\r")

            # averages over steps
            if Magnetization_flag:
                temps_mag = [[] for _ in range(len(sizes))]
                if Susceptibility_flag:
                    temps_sus = [[] for _ in range(len(sizes))]
            if Specific_flag:
                temps_E = [[] for _ in range(len(sizes))]
                temp_spec = [[] for _ in range(len(sizes))]

            systems = [System(size,T=Ti) for size in sizes]

            for system in systems:
                system.equilibrialization()

            steps = 100
            for i in range(steps):
                # print(f"\r[{'#'*round(j/(len(Ts) -1)*20):.<20}] {round(j/(len(Ts)-1)*100):02}% |{round((times[1]-times[0])*(len(Ts)-1-i)):04}s|",end="\r")
                for i in range(len(sizes)): # temp,system in zip(temps,systems):
                    systems[i].step()
                    if Magnetization_flag:
                        M = abs(systems[i].magnetization())
                        temps_mag[i].append(M)
                        if Susceptibility_flag:
                            temps_sus[i].append(M**2)
                    if Specific_flag:
                        E = systems[i].hamiltonian()
                        temps_E[i].append(E)
                        temp_spec[i].append(E**2)
            
            if Magnetization_flag:
                for s in range(len(sizes)):
                    mags[s].append(np.mean(temps_mag[s]))
                    if Susceptibility_flag:
                        mags_sq[s].append(np.mean(temps_sus[s]))
            
            if Specific_flag:
                for s in range(len(sizes)):
                    energys[s].append(np.mean(temps_E[s]))
                    specifics[s].append(np.mean(temp_spec[s]))

        if Magnetization_flag:
            for s in range(len(sizes)):
                magtot[s].append(np.mean(mags[s]))
                if Susceptibility_flag:
                    suscept[s].append(np.mean(mags_sq[s]-np.array(mags[s])**2))
        if Specific_flag:
            for s in range(len(sizes)):
                energy_tot[s].append(np.mean(energys[s]))
                specific_tot[s].append(np.mean(specifics[s])-np.array(energy_tot[s])**2)

    print("\nDone")

    # print(system32.accepted/system32.steps)
    plotting_LLsize(magtot, suscept, specific_tot)

def plotting_LLsize(Magnet, Suscept, Specific):
    sizes = [4,8,16,32]
    Ts = np.linspace(1,4,20)
    labels = [f"N={size}"for size in sizes]
    if len(Magnet[0])>0:
        for i,magt in enumerate(Magnet):
            plt.plot(Ts,magt,color = f"C{i}",linestyle="--")
            plt.scatter(Ts,magt,facecolors = "none",edgecolors = f"C{i}",marker=["o","s"][i%2],label=labels[i])
        plt.grid(linestyle="--")
        plt.legend()
        plt.title("Magnetization over temperature for different grid sizes")
        plt.ylabel("Magnetization |M|")
        plt.xlabel(r"Temperature $Tk_B/J$")
        plt.show()
    if len(Suscept[0])>0:
        for i,suscept in enumerate(Suscept):
            plt.plot(Ts,suscept,color = f"C{i}",linestyle="--")
            plt.scatter(Ts,suscept,facecolors = "none",edgecolors = f"C{i}",marker=["o","s"][i%2],label=labels[i])
        plt.grid(linestyle="--")
        plt.legend()
        plt.title("Susceptibility over temperature for different grid sizes")
        plt.ylabel(r"Susceptibility $\chi$")
        plt.xlabel(r"Temperature $Tk_B/J$")
        plt.show()
    if len(Specific[0])>0:
        for i,specific in enumerate(Specific):
            plt.plot(Ts,specific,color = f"C{i}",linestyle="--")
            plt.scatter(Ts,specific,facecolors = "none",edgecolors = f"C{i}",marker=["o","s"][i%2],label=labels[i])
        plt.grid(linestyle="--")
        plt.legend()
        plt.title("Specific heat over temperature for different grid sizes")
        plt.ylabel(r"Specific heat")
        plt.xlabel(r"Temperature $Tk_B/J$")
        plt.show()

#################################################################
## MAIN

def main():

    # LLsizes(Magnetization_flag=None, Susceptibility_flag=None, Specific_flag=1)

    # print(estimate_order(16,1,100,10000))

    # order_over_T(16,1,10000)

    # Ts = np.linspace(1,4,20)

    # magsup = []

    # for j,Ti in enumerate(Ts):

    #     temp = []
    #     for _ in range(1):

    #         system = System(32, T=Ti)

    #         system.equilibrialization()

    #         times = [time.time(), time.time()]

    #         steps = 100

    #         tt = list(range(0,steps))
    #         mags = []

    #         for i in range(steps):

    #             times.pop(0)
    #             times.append(time.time())
    #             # print(f"\r[{'#'*round(i/(steps -1)*20):.<20}] {round(i/(steps-1)*100):02}% |{round((times[1]-times[0])*(steps-1-i)):04}s|",end="\r")
    #             print(f"\r[{'#'*round(j/(len(Ts) -1)*20):.<20}] {round(j/(len(Ts)-1)*100):02}% |{round((times[1]-times[0])*(len(Ts)-1-i)):04}s|",end="\r")
    #             system.step()

    #             mags.append(system.magnetization())

    #             # if i in [0,9,99,999]:

    #             #     system.show_state()
    #         temp.append(abs(np.mean(mags)))
            
    #     magsup.append(np.mean(temp))

    # mag = []

    # system = System(32, T=2.2)

    # # system.equilibrialization()

    # times = [time.time(), time.time()]

    # # steps = 100000

    # system.equilibrialization()

    # steps = 10000

    # tt = list(range(0,steps))

    # for i in range(steps):

    #     times.pop(0)
    #     times.append(time.time())
    #     print(f"\r[{'#'*round(i/(steps -1)*20):.<20}] {round(i/(steps-1)*100):02}% |{round((times[1]-times[0])*(steps-1-i)):04}s|",end="\r")
    #     # print(f"\r[{'#'*round(j/(len(Ts) -1)*20):.<20}] {round(j/(len(Ts)-1)*100):02}% |{round((times[1]-times[0])*(len(Ts)-1-i)):04}s|",end="\r")
    #     system.step()

    #     mag.append(system.magnetization())


    # plt.plot(tt,mag)
    # plt.show()


    # print(system.magnetization())
    # system.lattice[0][0] = 2
    # system.lattice[-1][-1] = 2

    # print(system.hamiltonian())
    # system.show_state()

    # print(system.calc_dE(2,2))

    # print(system.lattice)
    # print(lattice.magnetization())

#################################################################
## RUN CODE

if __name__ == "__main__":
    main()