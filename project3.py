#################################################################
## PROJECT 3
#################################################################

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import Boltzmann as k

class System():

    def __init__(self,N, J, B=0, T=5) -> None:

        self.N = N
        self.J = J
        self.B = B
        self.T = T
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

    def show_state(self):
        plt.imshow(self.lattice,cmap="Greys")
        plt.xlabel("i")
        plt.ylabel("j")
        plt.show()


def main():
    system = System(32, 1, T=100)

    for _ in range(10000):
        system.step()

    system.show_state()
    # system.lattice[0][0] = 2
    # system.lattice[-1][-1] = 2

    # print(system.hamiltonian())
    # system.show_state()

    print(system.calc_dE(2,2))

    # print(system.lattice)
    # print(lattice.magnetization())

if __name__ == "__main__":
    main()