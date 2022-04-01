import numpy as np
import matplotlib.pyplot as plt

eps = 1
ro = 1

def U(r):
    U = 4*eps*((ro/r)**12 - (ro/r)**6)
    return U


def main():
    r = np.linspace(0.9, 3, 100)
    plt.plot(r, U(r))
    plt.grid(True)
    plt.xlabel("Distance")
    plt.ylabel("Potential")
    plt.title("Plot of potential")

    plt.show()




if __name__ == "__main__":
    main()
