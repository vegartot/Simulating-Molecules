import numpy as np, matplotlib.pyplot as plt

def getVel():
    with open("outputdata_vel.txt", 'r') as infile:
        N, endTime, timepoints = infile.readline()[1:].split()
        N, endTime, timepoints = int(N), float(endTime), int(timepoints)
        velocities = np.ndarray((N, timepoints), dtype=object)
        indexN = 0;
        indexT = 0;
        for i, line in enumerate(infile):
            if not line.startswith('*'):
                vx, vy, vz = line.split()
                vx, vy, vz = float(vx), float(vy), float(vz)
                velocities[indexN][indexT] = np.array([vx, vy, vz])
                indexN += 1
            else:
                indexN = 0;
                indexT += 1;
        t = np.linspace(0, endTime, timepoints)
        return t, velocities, N, endTime, timepoints
        
def getPot():
    with open("outputdata_pot.txt", 'r') as infile:
        endTime, timepoints = infile.readline().split()
        endTime, timepoints = float(endTime), int(timepoints)
        t = np.linspace(0, endTime, timepoints)
        U = np.zeros_like(t)
        for i, line in enumerate(infile):
            U[i] = float(line)
        return t, U

def getKinetic(velocities, N, timepoints):
    kinetic = np.zeros(int(timepoints))
    for i in range(timepoints):
        for p in velocities:
            kinetic[i] += np.linalg.norm(p[i])**2
    
    return 0.5 * kinetic

def getCorr(velocities, N, timepoints):
    correlation = np.zeros(timepoints)
    v0 = [np.linalg.norm(velocities[i][0]) for i in range(N)]
    for i in range(timepoints):
        for j, p in enumerate(velocities):
            vi = np.linalg.norm(p[i])
            vj0 = v0[j]
            correlation[i] += np.dot(vi, vj0) / vj0**2
    return 1/N * correlation

t, vel, N, endTime, timepoints = getVel()
t, U = getPot()
Ek = getKinetic(vel, N, timepoints)
C = getCorr(vel, N, timepoints)

plt.plot(t, U, label="Potential")
plt.plot(t, Ek, label="Kinetic")
plt.plot(t, U+Ek, label="Total")
plt.legend()
plt.show()

plt.plot(t, C)
plt.title("Velocity Correlation")
plt.show()

plt.plot(t, 2/N * Ek)
plt.title("Temperature")
plt.show()
