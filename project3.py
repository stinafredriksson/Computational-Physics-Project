#################################################################
## PROJECT 3
#################################################################

import numpy as np
import time
import matplotlib.pyplot as plt
from scipy.constants import Boltzmann as k

class System():

    def __init__(self,N, J, B=0, T=5) -> None:

        self.N = N
        self.J = J
        self.B = B
        self.T = T
        self.steps = 0
        # initializes the spin lattice with random spins
        self.lattice = np.random.choice([-1,1],size=(N,N))

    def magnetization(self):
        return np.sum(self.lattice)
    
    def hamiltonian(self):
        H = 0
        for i in range(self.N):
            for j in range(self.N):
                H += -self.J*self.lattice[i][j]*(self.lattice[(i+1)%self.N][j]+self.lattice[i-1][j]+self.lattice[i][(j+1)%self.N]+self.lattice[i][j-1])
                H += -self.B*np.sum(self.lattice)
        return H
    
    def calc_dE(self,i,j):

        # H_original = self.hamiltonian()
        # S_original = self.lattice[i][j]
        # self.lattice[i][j] = -S_original
        # H_flip = self.hamiltonian()
        # self.lattice[i][j] = -S_original
        # print()
        return 4*self.J*self.lattice[i][j]*(self.lattice[(i+1)%self.N][j]+self.lattice[i-1][j]+self.lattice[i][(j+1)%self.N]+self.lattice[i][j-1])

    def step(self):
        i, j = np.random.randint(self.N), np.random.randint(self.N)
        dE = self.calc_dE(i,j)

        if dE < 0 or np.random.rand() < np.exp(-dE/k/self.T):
            self.lattice[i][j] = - self.lattice[i][j]
        
        self.steps += 1

    def show_state(self):
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


def main():

    # print(estimate_order(16,1,100,10000))

    order_over_T(16,1,10000)

    # system = System(16, 1, T=100)

    # times = [time.time(), time.time()]

    # steps = 100000

    # for i in range(steps):

    #     times.pop(0)
    #     times.append(time.time())
    #     print(f"\r[{'#'*round(i/(steps -1)*20):.<20}] {round(i/(steps-1)*100):02}% |{round((times[1]-times[0])*(steps-1-i)):04}s|",end="\r")
    #     system.step()

    # print("\nDone")

    # system.show_state()


    # print(system.magnetization())
    # system.lattice[0][0] = 2
    # system.lattice[-1][-1] = 2

    # print(system.hamiltonian())
    # system.show_state()

    # print(system.calc_dE(2,2))

    # print(system.lattice)
    # print(lattice.magnetization())

if __name__ == "__main__":
    main()