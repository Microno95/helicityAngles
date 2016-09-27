import time as tm
import os, os.path
import shutil as shl
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
import rebound as rb
import sys as sys
import multiprocessing as mp
from numpy import arccos as acos
from math import atan2, fsum, sqrt

np.seterr(all='raise')
phaseSpaceParameters = ['a', 'e', 'omega', 'inc', 'M', 'Omega']
altphaseSpaceParameters = ['x', 'y', 'z', 'vx', 'vy', 'vz']
tSave = [0]

def ensure_dir(path): # Ensures that the path to a directory exists, and creates it if not
    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path): # Raises the OSError if and only if the specified path is a file and not a directory
            raise


def calcHelicity(sim1, sim2, currentDeviation, derivativeDeviation):
    calculatedOrbitsSim1 = sim1.particles
    calculatedOrbitsSim2 = sim2.particles
    currentDeviation.append([])
    derivativeDeviation.append([])
    for i in altphaseSpaceParameters:
        for j in range(1, 3):
            currentDeviation[-1].append(eval("calculatedOrbitsSim2[{}].{} - calculatedOrbitsSim1[{}].{}".format(j, i, j, i)))
            derivativeDeviation[-1].append((currentDeviation[-1][-1] - currentDeviation[-2][-1]) / ((sim1.dt + sim2.dt) / 2.0))


def helicityWrapper(reb_sim1, reb_sim2, reb_sim3, deviation1, deviation2, derivDeviation1, derivDeviation2, helicity1, twist1):
    tSave.append((reb_sim1.t + reb_sim2.t + reb_sim3.t) / 3.0)
    calcHelicity(reb_sim1, reb_sim2, deviation1, derivDeviation1)
    calcHelicity(reb_sim1, reb_sim3, deviation2, derivDeviation2)
    helicity1.append([])
    twist1.append([])
    heliMag1 = sqrt(fsum([i**2 for i in deviation1[-1]]))
    heliMag2 = sqrt(fsum([i**2 for i in deviation2[-1]]))
    twistMag1 = sqrt(fsum([i**2 for i in derivDeviation1[-1]]))
    twistMag2 = sqrt(fsum([i**2 for i in derivDeviation2[-1]]))
    dotHelicity = fsum([i * j for i, j in zip(deviation1[-1], deviation2[-1])])
    dotTwist = fsum([i * j for i, j in zip(derivDeviation1[-1], derivDeviation2[-1])])
    helicity1[-1].append(acos(dotHelicity / (heliMag1 * heliMag2)))
    twist1[-1].append(acos(dotTwist / (twistMag1 * twistMag2)) * 3 / (reb_sim1.dt + reb_sim2.dt + reb_sim3.dt))


def reboundHelicityWrapper(reb_sim):
    helicityWrapper(sim1, sim2, sim3, currentDeviation1, currentDeviation2, derivativeDeviation1, derivativeDeviation2, helicityAngles1, twistAngles1)

def magnitude(vari):
    length_vari = 0
    for i in vari.particles[:3]:
        length_vari += i.x**2 + i.y**2 + i.z**2
        length_vari += i.vx**2 + i.vy**2 + i.vz**2
    return np.sqrt(length_vari)


def main(a2=2.53, a1=0.861, endTime=6, part=1, sep=1.1, dt=0.01, recVar=False, avar1 = (0, 0), avar2 = (0, 0)):
    global sim1, sim2, sim3, currentDeviation1, currentDeviation2, derivativeDeviation1, derivativeDeviation2, helicityAngles1, twistAngles1
    currentDeviation1 = []
    currentDeviation2 = []
    helicityAngles1 = [[0]]
    twistAngles1 = [[0]]
    pi = np.pi
    dir_path = "./Variational Lyapunov/"
    ensure_dir(dir_path)
    t0 = tm.perf_counter()
    sim1 = rb.Simulation()
    sim1.units = ('Yr', 'AU', 'Msun')
    sim1.integrator = 'ias15'
    sim1.dt = dt
    sim1.add(m=1.31, x=0, y=0, z=0)
    sim1.add(primary=sim1.particles[0], a=a1 + avar1[0], m=14.57/1047.56, e=0.239, omega=(290*pi)/180, inc=(16.7*pi)/180, M=154.8*pi/180, Omega=(295.5*pi)/180)
    sim1.add(primary=sim1.particles[0], a=a2 + avar1[1], m=10.19/1047.56, e=0.274, omega=(240.8*pi)/180, inc=(13.5*pi)/180, M=82.5*pi/180, Omega=(115*pi)/180)
    sim1.move_to_com()
    sim2 = rb.Simulation()
    sim2.units = ('Yr', 'AU', 'Msun')
    sim2.integrator = 'ias15'
    sim2.dt = dt
    sim2.add(m=1.31, x=0, y=0, z=0)
    sim2.add(primary=sim1.particles[0], a=a1 + avar2[0], m=14.57/1047.56, e=0.239, omega=(290*pi)/180, inc=(16.7*pi)/180, M=154.8*pi/180, Omega=(295.5*pi)/180)
    sim2.add(primary=sim1.particles[0], a=a2 + avar2[1], m=10.19/1047.56, e=0.274, omega=(240.8*pi)/180, inc=(13.5*pi)/180, M=82.5*pi/180, Omega=(115*pi)/180)
    sim2.move_to_com()
    for i in range(1, len(sim2.particles)):
        currentDeviation1.append([getattr(sim2.particles[i], k) for k in altphaseSpaceParameters])
    sim3 = rb.Simulation()
    sim3.units = ('Yr', 'AU', 'Msun')
    sim3.integrator = 'ias15'
    sim3.dt = dt
    sim3.add(m=1.31, x=0, y=0, z=0)
    sim3.add(primary=sim1.particles[0], a=a1, m=14.57/1047.56, e=0.239, omega=(290*pi)/180, inc=(16.7*pi)/180, M=154.8*pi/180, Omega=(295.5*pi)/180)
    sim3.add(primary=sim1.particles[0], a=a2, m=10.19/1047.56, e=0.274, omega=(240.8*pi)/180, inc=(13.5*pi)/180, M=82.5*pi/180, Omega=(115*pi)/180)
    sim3.move_to_com()
    for i in range(1, len(sim3.particles)):
        currentDeviation2.append([getattr(sim2.particles[i], k) for k in altphaseSpaceParameters])
    derivativeDeviation1 = [[0] * len(currentDeviation1)]
    derivativeDeviation2 = [[0] * len(currentDeviation2)]
    print("All systems go for dt: {} and a2: {}!".format(dt, a2))
    m = 0
    t_x = np.linspace(0, 10**endTime, 100000)
    print("hot potato")
    for t in t_x:
        global tCurrent
        tCurrent = t
        try:
            sim1.integrate(t)
            sim2.integrate(t)
            sim3.integrate(t)
            reboundHelicityWrapper(sim1)
        except KeyboardInterrupt:
            print(" ")
            quit()
        else:
            if m%10 == 0:
                print("t: {:.2%} Completed - Sim_mag: {:.1e};{:.1e};{:.1e}".format(t / t_x[-1], magnitude(sim1), magnitude(sim2), magnitude(sim3)), end="\r")
                sys.stdout.flush()
            m += 1
    print("Cold potato")
    with open(dir_path + "HelicityTest - {}.txt".format(*initDeviationSimOne), 'w') as fle:
        toSave = (tSave, currentDeviation1, currentDeviation2, derivativeDeviation1, derivativeDeviation2, helicityAngles1, twistAngles1)
        for i in zip(*toSave):
            baseLine = '; '
            line = ["{}".format(j) for j in i]
            fle.write(baseLine.join(line) + '\n')
    tcomp = tm.perf_counter() - t0
    print("Completed in {:.2f} minutes".format(tcomp/60) + " for system a2: {}!".format(a2))

def worker(arguments):
    main(**arguments)


if __name__ == "__main__":
    worker({'a2': 2.53, 'avar1': tuple([0.1, 0]), 'avar2': tuple([0, 0.001]), 'endTime': 6})

