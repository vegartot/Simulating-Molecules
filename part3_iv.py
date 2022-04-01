import numpy as np, matplotlib.pyplot as plt


U = lambda r: 4*(pow(r, -12) - pow(r, -6))-4*(pow(3,-12)-pow(3,-6))

r = np.linspace(.9, 3.5, 200)
plt.plot(r, U(r))
plt.ylabel("Potential")
plt.xlabel("Distance [sigma]")
plt.title("Shifted potential")
plt.grid(True)
plt.show()

print(U(pow(2, 1/6)))
