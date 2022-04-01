import numpy as np, matplotlib.pyplot as plt, itertools


class Particle:
    def __init__(self, x=0, y=0, z=0, vx=0, vy=0, vz=0):
        self.x = x
        self.y = y
        self.z = z
        self.vx = vx
        self.vy = vy
        self.vz = vz
        self.coord = np.array([self.x, self.y, self.z], float)
        self.velocity = np.array([self.vx, self.vy, self.vz], float)
        self.pos = []
        self.vel = []
        self.acc = []


    def distance(self, other, t):
        dr = self.pos[t] - other.pos[t]
        return np.linalg.norm(dr)

class Solver:
    """
    def __init__(self, N):
        self.N = N
        self.particles = []
        for i in range(N):
            self.particles.append(Particle())
        self.particles = np.array(self.particles)
    """
    def __init__(self, N):
        self.N = N
        self.particles = np.ndarray(N, dtype=Particle)
        for i, p in enumerate(self.particles):
            self.particles[i] = Particle()

    def set_initial_positions(self, *args):
        if len(args) != self.N:
            raise IndexError("Initial conditions not matching number of particles.")
        for i, arg in enumerate(args):
            for j, coord in enumerate(arg):
                self.particles[i].coord[j] = coord

    def set_initial_velocities(self, *args):
        if len(args) != self.N:
            raise IndexError("Initial conditions not matching number of particles.")
        for i, arg in enumerate(args):
            for j, velocity in enumerate(arg):
                self.particles[i].velocity[j] = velocity

    def solve(self, timepoints):
        for i, p in enumerate(self.particles):
            p.pos = np.zeros((timepoints.size, 3))
            p.vel = np.zeros((timepoints.size, 3))
            p.acc = np.zeros((timepoints.size, 3))
            p.pos[0] = p.coord
            p.vel[0] = p.velocity
            p.acc[0] = np.zeros(3)

        for k in range(timepoints.size - 1):
            dt = timepoints[k+1] - timepoints[k]

            for i in self.particles:
                i.pos[k+1] = i.pos[k] + i.vel[k]*dt + 1/2*i.acc[k]*dt**2

            for i, j in itertools.combinations(self.particles, 2):
                dr = i.distance(j, k+1)
                if dr < 3:
                    acc_new = 24*(2*dr**(-12) -dr**(-6))*(i.pos[k+1] - j.pos[k+1])/(dr**2)
                    i.acc[k+1] += acc_new
                    j.acc[k+1] -= acc_new


            for i in self.particles:
                i.vel[k+1] = i.vel[k] + 1/2*(i.acc[k] + i.acc[k+1])*dt


def write_to_xyz(particles, nr):
    N = len(particles)
    timesteps = len(particles[0].pos)
    with open(f"part3_{nr}_out.txt", "w") as infile:
        for i in range(timesteps):
            infile.write("%d\n" %N)
            infile.write("type x y z\n")
            for j in particles:
                infile.write("Ar %f %f %f\n" %(j.pos[i][0], j.pos[i][1], j.pos[i][2]))

def kinetic_energy(particles):
    vel_new = 0
    for i in particles:
        s = np.sum(i.vel**2, axis=1)
        vel_new += s
    return 0.5 * vel_new

def potential_energy(particles):
    timesteps = len(particles[0].pos)
    position = np.zeros(timesteps)
    for k in range(timesteps):
        U_tot = 0
        distance = 0
        for i, j in itertools.combinations(particles, 2):
            distance = i.distance(j, k)
            U_tot += 4*(distance**(-12) - distance**(-6))
        position[k] = U_tot
    return position

if __name__ == "__main__":

    solve = Solver(2)
    t = np.linspace(0, 5, 500)
    solve.set_initial_positions([1.5,0,0], [0,0,0])
    solve.solve(t)
    write_to_xyz(solve.particles, 1)

    plt.plot(t, kinetic_energy(solve.particles), label="Kinetic Energy")
    plt.plot(t, potential_energy(solve.particles), label="Potential Energy")
    plt.plot(t, kinetic_energy(solve.particles) + potential_energy(solve.particles), label="Total Energy")

    plt.title("Initial position: (%.1f, %.1f, %.1f)\nN = %d atoms" %(solve.particles[0].pos[0][0],\
    solve.particles[0].pos[0][1], solve.particles[0].pos[0][2], 2))

    plt.xlabel("Time"); plt.ylabel("Energy")
    plt.legend()
    plt.tight_layout()
    plt.show()

    solve = Solver(2)
    t = np.linspace(0, 5, 500)
    solve.set_initial_positions([.95,0,0], [0,0,0])
    solve.solve(t)
    write_to_xyz(solve.particles, 1)

    plt.plot(t, kinetic_energy(solve.particles), label="Kinetic Energy")
    plt.plot(t, potential_energy(solve.particles), label="Potential Energy")
    plt.plot(t, kinetic_energy(solve.particles) + potential_energy(solve.particles), label="Total Energy")

    plt.title("Initial position: (%.1f, %.1f, %.1f)\nN = %d atoms" %(solve.particles[0].pos[0][0],\
    solve.particles[0].pos[0][1], solve.particles[0].pos[0][2], 2))

    plt.xlabel("Time"); plt.ylabel("Energy")
    plt.legend()
    plt.tight_layout()
    plt.show()


    # 3b)
    solve = Solver(4)
    solve.set_initial_positions([1,0,0], [0,1,0], [-1,0,0], [0,-1,0])
    solve.set_initial_velocities([0,0,0], [0,0,0], [0,0,0], [0,0,0])
    solve.solve(t)
    write_to_xyz(solve.particles, 2)

    plt.plot(t, kinetic_energy(solve.particles), label="Kinetic Energy")
    plt.plot(t, potential_energy(solve.particles), label="Potential Energy")
    plt.plot(t, kinetic_energy(solve.particles) + potential_energy(solve.particles), label="Total Energy")

    plt.title("Initial position: (%.1f, %.1f, %.1f)\nN = %d atoms" %(solve.particles[0].pos[0][0],\
    solve.particles[0].pos[0][1], solve.particles[0].pos[0][2], 4))

    plt.xlabel("Time"); plt.ylabel("Energy")
    plt.legend()
    plt.tight_layout()
    plt.show()

    # 3b)
    solve = Solver(4)
    solve.set_initial_positions([1,0.1,0], [0,1,0], [-1,0,0], [0,-1,0])
    solve.set_initial_velocities([0,0,0], [0,0,0], [0,0,0], [0,0,0])
    solve.solve(t)
    write_to_xyz(solve.particles, 2)

    plt.plot(t, kinetic_energy(solve.particles), label="Kinetic Energy")
    plt.plot(t, potential_energy(solve.particles), label="Potential Energy")
    plt.plot(t, kinetic_energy(solve.particles) + potential_energy(solve.particles), label="Total Energy")

    plt.title("Initial position: (%.1f, %.1f, %.1f)\nN = %d atoms" %(solve.particles[0].pos[0][0],\
    solve.particles[0].pos[0][1], solve.particles[0].pos[0][2], 4))

    plt.xlabel("Time"); plt.ylabel("Energy")
    plt.legend()
    plt.tight_layout()
    plt.show()
