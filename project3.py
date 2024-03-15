#################################################################
## PROJECT 3
#################################################################
## IMPORTS

import numpy as np
import time
import csv
import matplotlib.pyplot as plt
from scipy.constants import Boltzmann as k

#################################################################
## CLASSES

class System():
    """This class handles the low-temperature behavior of a 
    lattice material"""

    def __init__(self, N: int, J: float = -k, B: float = 0, T: float = 2) -> None:
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
        self.lattice = np.random.choice([-1,1],size=(N,N))


    def magnetization(self) -> float:
        """calculates the magnetization of the lattice"""

        return np.sum(self.lattice)
    

    def hamiltonian(self) -> float:
        """calculates the energy of the lattice"""

        H = 0
        for i in range(self.N):
            for j in range(self.N):
                H += -self.J*self.lattice[i][j]*(self.lattice[(i+1)%self.N][j]+self.lattice[i-1][j]+self.lattice[i][(j+1)%self.N]+self.lattice[i][j-1])
                
                # this if-statement saves a lot of computing time when B=0
                if self.B:
                    H += -self.B*np.sum(self.lattice)
        return H


    def calc_dE(self, i: int, j: int) -> float:
        """calculates the change in energy from a spin flip at i,j"""

        return 2*self.J*self.lattice[i][j]*(self.lattice[(i+1)%self.N][j]+self.lattice[i-1][j]+self.lattice[i][(j+1)%self.N]+self.lattice[i][j-1])


    def step(self) -> None:
        """a step in the Metropolis algorithm using a random next site"""

        # random next site
        i, j = np.random.randint(self.N), np.random.randint(self.N)
        
        dE = self.calc_dE(i,j)

        # accept the move or not
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
        """shows a figure of the state with binary spins"""
        
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


def LLsizes(sizes: list[int] = [4,8,16,32], n_avgs: int = 100,
            flags: list[bool] = [True, False, False, False]) -> None:
    
    """Observe that this function is also written in C++"""

    Magnetization_flag = flags[0]
    Susceptibility_flag = flags[1]
    Specific_flag = flags[2]
    Cumulant_flag = flags[3]

    Ts = np.linspace(1,4,20)

    times = [time.time(), time.time()]
    
    magtot = [[] for _ in range(len(sizes))]
    sus_tot =  [[] for _ in range(len(sizes))]
    energy_tot = [[] for _ in range(len(sizes))]
    specific_tot = [[] for _ in range(len(sizes))]
    cumul_tot = [[] for _ in range(len(sizes))]

    for t,Ti in enumerate(Ts):

        # averages over iterations
        if Magnetization_flag or Susceptibility_flag:
            mags = [[] for _ in range(len(sizes))]
        if Susceptibility_flag:
            sus = [[] for _ in range(len(sizes))]
        if Specific_flag:
            energys = [[] for _ in range(len(sizes))]
            specifics = [[] for _ in range(len(sizes))]
        if Cumulant_flag:
            cumulants = [[] for _ in range(len(sizes))]

        for j in range(n_avgs):
            times.pop(0)
            times.append(time.time())
            print(f"\r[{'#'*round((j+t*n_avgs)/(n_avgs*len(Ts) -1)*20):.<20}] {round((j+t*n_avgs)/(n_avgs*len(Ts)-1)*100):03}% |{round((times[1]-times[0])*(n_avgs*len(Ts)-1-j-t*n_avgs)):04}s|",end="\r")

            # averages over steps
            if Magnetization_flag or Susceptibility_flag:
                temps_mag = [[] for _ in range(len(sizes))]
                if Susceptibility_flag:
                    temps_sq = [[] for _ in range(len(sizes))]
            if Specific_flag:
                temps_E = [[] for _ in range(len(sizes))]
                temp_spec = [[] for _ in range(len(sizes))]
            if Cumulant_flag:
                temps_m2 = [[] for _ in range(len(sizes))]
                temps_m4 = [[] for _ in range(len(sizes))]

            systems = [System(size,T=Ti) for size in sizes]

            # lets the system approach thermal equilibrium
            for system in systems:
                system.equilibrialization()

            steps = 40

            ## runs steps in the simulation to sample data
            ########
            for step in range(steps):
                for s in range(len(sizes)): # temp,system in zip(temps,systems):
                    systems[s].equilibrialization(Nsweeps=4)
                    if Magnetization_flag or Susceptibility_flag:
                        M = abs(systems[s].magnetization())
                        temps_mag[s].append(M)
                        if Susceptibility_flag:
                            temps_sq[s].append(M**2)
                    if Specific_flag:
                        E = systems[s].hamiltonian()
                        temps_E[s].append(E)
                        temp_spec[s].append(E**2)
                    if Cumulant_flag:
                        M = abs(systems[s].magnetization())
                        temps_m2[s].append(M**2)
                        temps_m4[s].append(M**4)
            ########

            ## averages the samples taken during the steps
            ########
            if Magnetization_flag or Susceptibility_flag:
                for s in range(len(sizes)):
                    mags[s].append(np.mean(temps_mag[s]))
                    if Susceptibility_flag:
                        sus[s].append((np.mean(temps_sq[s])-np.mean(temps_mag[s])**2)/systems[s].N**2)
            if Specific_flag:
                for s in range(len(sizes)):
                    energys[s].append(np.mean(temps_E[s])/k)
                    specifics[s].append(np.mean((temp_spec[s])-np.mean(temps_E[s])**2)/k)
            if Cumulant_flag:
                for s in range(len(sizes)):
                    cumulants[s].append(1-np.mean(temps_m4[s])/(3*np.mean(temps_m2[s])**2))
            ########

        ## averages over a number of realizations
        ########
        if Magnetization_flag or Susceptibility_flag:
            for s in range(len(sizes)):
                magtot[s].append(np.mean(mags[s]))
                if Susceptibility_flag:
                    sus_tot[s].append(np.mean(sus[s]))
        if Specific_flag:
            for s in range(len(sizes)):
                energy_tot[s].append(np.mean(energys[s]))
                specific_tot[s].append(np.mean(specifics[s]))
        if Cumulant_flag:
            for s in range(len(sizes)):
                cumul_tot[s].append(np.mean(cumulants[s]))


    print("\nDone")

    # print(system32.accepted/system32.steps)
    plotting_LLsize(magtot, sus_tot, specific_tot,cumul_tot, sizes, Ts)


def plotting_LLsize(Magnet, Suscept, Specific, Cumul, sizes: list[int], Ts:list[float]):
    def __ferro(T):
        return 1/(T-Tc)
    
    Tc = 2.2691853
    x = np.linspace(Tc+0.1,Ts[-1],200)
    ferros = __ferro(x)

    labels = [f"N={size}"for size in sizes]
    ms = ["^","s","o","D"]
    if len(Magnet[0])>0:
        plt.figure(1)
        plt.vlines([Tc],0,1,linestyles="--", colors="grey", label=r"$T_c$")
        for i,magt in enumerate(Magnet):
            plt.plot(Ts,magt,color = f"C{i}",linestyle="--",marker=ms[i],markerfacecolor="none",label=labels[i])
            # plt.scatter(Ts,magt,facecolors = "none",edgecolors = f"C{i}",marker=["o","s"][i%2],label=labels[i])
        plt.grid(linestyle="--")
        plt.legend()
        plt.title("Magnetization over temperature for different grid sizes")
        plt.ylabel(r"Magnetization $|M|/N^2$")
        plt.xlabel(r"Temperature $Tk_B/J$")
    if len(Suscept[0])>0:
        plt.figure(2)
        plt.vlines([Tc],0,max(Suscept[-1]),linestyles="--", colors="grey", label=r"$T_c$")
        for i,suscept in enumerate(Suscept):
            plt.plot(Ts,suscept,color = f"C{i}",linestyle="--",marker=ms[i],markerfacecolor="none",label=labels[i])
            # plt.scatter(Ts,suscept,facecolors = "none",edgecolors = f"C{i}",marker=["o","s"][i%2],label=labels[i])
        plt.grid(linestyle="--")
        plt.legend()
        plt.title("Susceptibility over temperature for different grid sizes")
        plt.ylabel(r"Susceptibility $\chi/N^2$")
        plt.xlabel(r"Temperature $Tk_B/J$")
        # plt.plot(x,ferros)
    if len(Specific[0])>0:
        plt.figure(3)
        plt.vlines([Tc],0,max(Specific[-1]),linestyles="--", colors="grey", label=r"$T_c$")
        for i,specific in enumerate(Specific):
            plt.plot(Ts,specific,color = f"C{i}",linestyle="--",marker=ms[i],markerfacecolor="none",label=labels[i])
            # plt.scatter(Ts,specific,facecolors = "none",edgecolors = f"C{i}",marker=["o","s"][i%2],label=labels[i])
        plt.grid(linestyle="--")
        plt.legend()
        plt.title("Specific heat over temperature for different grid sizes")
        plt.ylabel(r"Specific heat $C_B/k_BN^2$")
        plt.xlabel(r"Temperature $Tk_B/J$")
    if len(Cumul[0])>0:
        plt.figure(4)
        plt.vlines([Tc],min(Cumul[-1]),1.2*max(Cumul[-1]),linestyles="--", colors="grey", label=r"$T_c$")
        for i,cum in enumerate(Cumul):
            plt.plot(Ts,cum,color = f"C{i}",linestyle="--",marker=ms[i],markerfacecolor="none",label=labels[i])
            # plt.scatter(Ts,cum,facecolors = "none",edgecolors = f"C{i}",marker=["o","s"][i%2],label=labels[i])
        plt.grid(linestyle="--")
        plt.legend()
        plt.title("Cumulant over temperature for different grid sizes")
        plt.ylabel(r"Cumulant $U^4_L$")
        plt.xlabel(r"Temperature $Tk_B/J$")
    plt.show()


def read_csv(filename):
    T = []
    M = []
    X = []
    E = []
    Cv = []
    U4 = []

    with open(filename, "r") as csvfile:
        reader = csv.reader(csvfile)

        next(reader)

        for row in reader:
            T.append(float(row[0]))
            M.append(float(row[1]))
            X.append(float(row[2]))
            E.append(float(row[3]))
            Cv.append(float(row[4]))
            U4.append(float(row[5]))
    return T,M,X,E,Cv,U4

def plot_CPP(sizes):

    Ms = []
    Xs = []
    Es = []
    Cvs = []
    U4s = []

    for size in sizes:
        filename = f"Ferro_80T/output_N{size}.csv"
        T,M,X,E,Cv,U4 = read_csv(filename)
        Ms.append(M)
        Xs.append(X)
        Es.append(E)
        Cvs.append(Cv)
        U4s.append(U4)
    
    plotting_LLsize(Ms,Xs,Cvs,U4s,sizes,T)


def show_evolution():

    sys = System(32,-k,0,1)

    n_steps = 200000

    for i in range(n_steps):
        sys.step()
        if i in [0,0.5*n_steps,n_steps-1]:
    # sys.equilibrialization(Nsweeps=100)
    
            sys.show_state()

#################################################################
## MAIN

def main():

    # LLsizes(flags=[False,True,False,False], n_avgs=2)

    sizes = [4,8,16,32]

    plot_CPP(sizes)


#################################################################
## RUN CODE

if __name__ == "__main__":
    main()