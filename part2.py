import numpy as np, matplotlib.pyplot as plt


class Solver:
    def __init__(self, function):       # Provide a derivative of diff. eq. as function of t, u
        self.f = function

    def set_initial_condition(self, acc0, vel0, pos0):
        self.acc0 = float(acc0)
        self.vel0 = float(vel0)
        self.pos0 = float(pos0)

    def advance(self):
        pass

    def solve(self, timepoints):
        self.t = timepoints
        self.acc = np.zeros_like(self.t)
        self.vel = np.zeros_like(self.t)
        self.pos = np.zeros_like(self.t)
        self.acc[0], self.vel[0], self.pos[0] = self.acc0, self.vel0, self.pos0

        for k in range(len(self.t)-1):
            self.k = k
            self.acc[self.k+1], self.vel[self.k+1], self.pos[self.k+1] = self.advance()
        return self.t, self.acc, self.vel, self.pos

class Euler(Solver):
    def advance(self):
        f = self.f; k = self.k; t = self.t; acc = self.acc; vel = self.vel; pos = self.pos;
        dt = t[k+1] - t[k]
        acc_new = f(pos[k])
        vel_new = vel[k] + acc[k]*dt
        pos_new = pos[k] + vel[k]*dt
        return acc_new, vel_new, pos_new

class EulerCromer(Solver):
    def advance(self):
        f = self.f; k = self.k; t = self.t; acc = self.acc; vel = self.vel; pos = self.pos;
        dt = t[k+1] - t[k]
        acc_new = f(pos[k])
        vel_new = vel[k] + acc_new*dt
        pos_new = pos[k] + vel_new*dt
        return acc_new, vel_new, pos_new

class VelocityVerlet(Solver):
    def advance(self):
        f = self.f; k = self.k; t = self.t; acc = self.acc; vel = self.vel; pos = self.pos;
        dt = t[k+1] - t[k]
        pos_new = pos[k] + vel[k]*dt + 0.5*acc[k]*dt**2
        acc_new = f(pos_new)
        vel_new = vel[k] + 0.5*(acc[k] + acc_new)*dt
        return acc_new, vel_new, pos_new

def write_to_xyz(pos):
    # Provide array of distance between the atoms
    # Assume momentum only along x-axis.
    with open("part2_out.txt", "w") as infile:
        for i, p in enumerate(pos):
            infile.write("2\n")
            infile.write("type x y z\n")
            infile.write(f"Ar {-p/2} {0} {0}\n")
            infile.write(f"Ar {p/2} {0} {0}\n")





# Define:
f = lambda dr: 24*(2*abs(dr)**(-12) - abs(dr)**(-6)) * (dr)/(abs(dr)**2)
def KineticEnergy(velocity):
    return 1/2*velocity**2

def PotentialEnergy(distance):
    return 4*((distance)**(-12) - (distance)**(-6))

if __name__ == "__main__":
    dt = 0.01
    timepoints = int(5/dt + 1)
    t = np.linspace(0, 5, timepoints)

    # 2b)
    # Solve case 1:
    solver1 = EulerCromer(f)
    solver1.set_initial_condition(0, 0, 1.5)
    dt = int(5/0.01 + 1)
    t1, a1, v1, p1 = solver1.solve(t)

    # Solve case 2:
    solver2 = EulerCromer(f)
    solver2.set_initial_condition(0, 0, 0.95)
    t2, a2, v2, p2 = solver2.solve(t)

    # Plotting:
    plt.subplot(2,1, 1)
    plt.plot(t1, p1)
    plt.title("Case 1")
    plt.ylabel("Distance")
    plt.xlabel("Time")

    plt.subplot(2,1,2)
    plt.plot(t2, p2)
    plt.title("Case 2")
    plt.ylabel("Distance")
    plt.xlabel("Time")

    plt.tight_layout()
    plt.show()

    # 2c. i)
    # Case 1:
    plt.subplot(2,3,1)
    plt.plot(t1, KineticEnergy(v1))
    plt.ylabel("Kinetic Energy")
    plt.xlabel("Time")

    plt.subplot(2,3,2)
    plt.plot(t1, PotentialEnergy(p1))
    plt.ylabel("Potential Energy")
    plt.xlabel("Time")
    plt.title("Case 1")

    plt.subplot(2,3,3)
    plt.plot(t1, PotentialEnergy(p1) + KineticEnergy(v1))
    plt.ylabel("Total Energy")
    plt.xlabel("Time")

    # Case 2:
    plt.subplot(2,3,4)
    plt.plot(t2, KineticEnergy(v2))
    plt.ylabel("Kinetic Energy")
    plt.xlabel("Time")

    plt.subplot(2,3,5)
    plt.plot(t2, PotentialEnergy(p2))
    plt.ylabel("Potential Energy")
    plt.xlabel("Time")
    plt.title("Case 2")

    plt.subplot(2,3,6)
    plt.plot(t2, PotentialEnergy(p2) + KineticEnergy(v2))
    plt.ylabel("Total Energy")
    plt.xlabel("Time")

    plt.tight_layout()
    plt.show()

    # 2c. iv)
    # Case 1:
    solver3 = Euler(f)
    solver3.set_initial_condition(0, 0, 1.5)
    solver5 = VelocityVerlet(f)
    solver5.set_initial_condition(0, 0, 1.5)
    t3, a3, v3, p3 = solver3.solve(t)
    t5, a5, v5, p5 = solver5.solve(t)

    # Case 2:
    solver4 = Euler(f)
    solver4.set_initial_condition(0, 0, 0.95)
    solver6 = VelocityVerlet(f)
    solver6.set_initial_condition(0, 0, 0.95)
    t4, a4, v4, p4 = solver4.solve(t)
    t6, a6, v6, p6 = solver6.solve(t)

    # Plotting:
    # Case 1:
    plt.subplot(2, 3, 1)
    plt.plot(t1, PotentialEnergy(p1) + KineticEnergy(v1))
    plt.title("Euler-Cromer")
    plt.ylabel("Potential Energy")
    plt.xlabel("Time")

    plt.subplot(2, 3, 2)
    plt.plot(t3, PotentialEnergy(p3) + KineticEnergy(v3))
    plt.title("Case 1:\nEuler")
    plt.ylabel("Potential Energy")
    plt.xlabel("Time")

    plt.subplot(2, 3, 3)
    plt.plot(t5, PotentialEnergy(p5) + KineticEnergy(v5))
    plt.title("Velocity-Verlet")
    plt.ylabel("Potential Energy")
    plt.xlabel("Time")

    # Case 2:
    plt.subplot(2, 3, 4)
    plt.plot(t1, PotentialEnergy(p2) + KineticEnergy(v2))
    plt.title("Euler-Cromer")
    plt.ylabel("Potential Energy")
    plt.xlabel("Time")

    plt.subplot(2, 3, 5)
    plt.plot(t3, PotentialEnergy(p4) + KineticEnergy(v4))
    plt.title("Case 2:\nEuler")
    plt.ylabel("Potential Energy")
    plt.xlabel("Time")

    plt.subplot(2, 3, 6)
    plt.plot(t5, PotentialEnergy(p6) + KineticEnergy(v6))
    plt.title("Velocity-Verlet")
    plt.ylabel("Potential Energy")
    plt.xlabel("Time")

    plt.tight_layout()
    plt.show()

    # Re-define:
    dt = 0.00001
    timepoints = int(5/dt + 1)
    t = np.linspace(0, 1.5, timepoints)

    # Solve case 1:
    t1, a1, v1, p1 = solver1.solve(t)
    t2, a2, v2, p2 = solver2.solve(t)

    t3, a3, v3, p3 = solver3.solve(t)
    t4, a4, v4, p4 = solver4.solve(t)

    t5, a5, v5, p5 = solver5.solve(t)
    t6, a6, v6, p6 = solver6.solve(t)

    # Plotting:
    # Case 1:
    plt.subplot(2, 3, 1)
    plt.plot(t1, PotentialEnergy(p1) + KineticEnergy(v1))
    plt.title("Euler-Cromer")
    plt.ylabel("Potential Energy")
    plt.xlabel("Time")

    plt.subplot(2, 3, 2)
    plt.plot(t3, PotentialEnergy(p3) + KineticEnergy(v3))
    plt.title("Case 1:\nEuler")
    plt.ylabel("Potential Energy")
    plt.xlabel("Time")

    plt.subplot(2, 3, 3)
    plt.plot(t5, PotentialEnergy(p5) + KineticEnergy(v5))
    plt.title("Velocity-Verlet")
    plt.ylabel("Potential Energy")
    plt.xlabel("Time")

    # Case 2:
    plt.subplot(2, 3, 4)
    plt.plot(t1, PotentialEnergy(p2) + KineticEnergy(v2))
    plt.title("Euler-Cromer")
    plt.ylabel("Potential Energy")
    plt.xlabel("Time")

    plt.subplot(2, 3, 5)
    plt.plot(t3, PotentialEnergy(p4) + KineticEnergy(v4))
    plt.title("Case 2:\nEuler")
    plt.ylabel("Potential Energy")
    plt.xlabel("Time")

    plt.subplot(2, 3, 6)
    plt.plot(t5, PotentialEnergy(p6) + KineticEnergy(v6))
    plt.title("Velocity-Verlet")
    plt.ylabel("Potential Energy")
    plt.xlabel("Time")

    plt.tight_layout()
    plt.show()

    write_to_xyz(p1)
