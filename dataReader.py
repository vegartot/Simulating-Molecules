import numpy as np, matplotlib.pyplot as plt

def getVel():
    with open("outputdata_vel.txt", 'r') as infile:
        N, endTime, timepoints = infile.readline()[1:].split()
        N, endTime, timepoints = int(N), float(endTime), int(timepoints)
        velocities = []
        vel = 0
        for i, line in enumerate(infile):
            if not line.startswith('*'):
                vx, vy, vz = line.split()
                vx, vy, vz = float(vx), float(vy), float(vz)
                vel += vx**2 + vy**2 + vz**2
            else:
                velocities.append(vel)
                vel = 0
        velocities.append(vel)
        velocities = np.array(velocities)
        t = np.linspace(0, endTime, timepoints)
        return t, velocities
        
def getPot():
    with open("outputdata_pot.txt", 'r') as infile:
        endTime, timepoints = infile.readline().split()
        endTime, timepoints = float(endTime), int(timepoints)
        t = np.linspace(0, endTime, timepoints)
        U = np.zeros_like(t)
        for i, line in enumerate(infile):
            U[i] = float(line)
        return t, U

t1, vel = getVel()
t2, U = getPot()

plt.plot(t1, vel + U)
plt.show()

plt.plot(t2, U)
plt.show()