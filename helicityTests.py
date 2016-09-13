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
from math import atan2

np.seterr(all='raise')


def ensure_dir(path): # Ensures that the path to a directory exists, and creates it if not
    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path): # Raises the OSError if and only if the specified path is a file and not a directory
            raise


def calcHelicity(sim1, sim2, currentDeviation, derivativeDeviation):
    calculatedOrbitsSim1 = sim1.calculate_orbits()
    calculatedOrbitsSim2 = sim2.calculate_orbits()
    currentDeviation.append([])
    derivativeDeviation.append([])
    helicityAngles.append([])
    twistAngles.append([])
    tSave.append((sim1.t + sim2.t) / 2.0)
    for i in phaseSpaceParameters:
        for j in range(2):
            currentDeviation[-1].append(eval("calculatedOrbitsSim2[{}].{} - calculatedOrbitsSim1[{}].{}".format(j, i, j, i)))
            derivativeDeviation[-1].append(currentDeviation[-1][-1] / ((sim1.dt + sim2.dt) / 2.0))
    for i in range(len(currentDeviation[-1])):
        for j in range(len(currentDeviation[-1])):
            if j != i:
                helicityAngles[-1].append(atan2(currentDeviation[-1][i], currentDeviation[-1][j]))
                derivDevIJ = [derivativeDeviation[-1][i], derivativeDeviation[-1][j]]
                devIJ = [currentDeviation[-1][i], currentDeviation[-1][j]]
                if norm(devIJ) != 0.0:
                    devIJ = [k / norm(devIJ) for k in devIJ]
                else:
                    devIJ = [0, 0]
                if norm(derivDevIJ) != 0.0:
                    derivDevIJ = [k / norm(derivDevIJ) for k in derivDevIJ]
                else:
                    derivDevIJ = [0, 0]
                dotProdij = np.sqrt(np.dot(derivDevIJ, devIJ))
                try:
                    angProdij = acos(dotProdij) / ((sim1.dt + sim2.dt) / 2.0)
                except FloatingPointError:
                    angProdij = 0.0
                twistAngles[-1].append(angProdij)


def helicityWrapper(reb_sim):
    calcHelicity(sim1, sim2, currentDeviation, derivativeDeviation)


def magnitude(vari):
    length_vari = 0
    for i in vari.particles[:3]:
        length_vari += i.x**2 + i.y**2 + i.z**2
        length_vari += i.vx**2 + i.vy**2 + i.vz**2
    return np.sqrt(length_vari)


def main(a2=2.53, a1=0.861, endTime=6, part=1, sep=1.1, dt=0.01, recVar=False, initDeviation=(0, 0, )):
    global currentDeviation
    currentDeviation = [list(initDeviation[:])]
    global derivativeDeviation
    derivativeDeviation = [[0] * len(initDeviation)]
    global helicityAngles
    helicityAngles = [[1] + [0] * (len(initDeviation)**2 - len(initDeviation))]
    global twistAngles
    twistAngles = [[0] * (len(initDeviation)**2 - len(initDeviation))]
    pi = np.pi
    dir_path = "./Variational Lyapunov/"
    ensure_dir(dir_path)
    t0 = tm.perf_counter()
    global sim1
    sim1 = rb.Simulation()
    sim1.units = ('Yr', 'AU', 'Msun')
    sim1.integrator = 'ias15'
    sim1.dt = dt
    sim1.add(m=1.31, x=0, y=0, z=0)
    sim1.add(primary=sim1.particles[0], a=a1, m=14.57/1047.56, e=0.239, omega=(290*pi)/180, inc=(16.7*pi)/180, M=154.8*pi/180, Omega=(295.5*pi)/180)
    sim1.add(primary=sim1.particles[0], a=a2, m=10.19/1047.56, e=0.274, omega=(240.8*pi)/180, inc=(13.5*pi)/180, M=82.5*pi/180, Omega=(115*pi)/180)
    sim1.move_to_com()
    global sim2
    sim2 = rb.Simulation()
    sim2.units = ('Yr', 'AU', 'Msun')
    sim2.integrator = 'ias15'
    sim2.dt = dt
    sim2.add(m=1.31, x=0, y=0, z=0)
    sim2.add(primary=sim2.particles[0], a=a1 + initDeviation[0], m=14.57/1047.56, e=0.239 + initDeviation[2], omega=((290 + initDeviation[4])*pi)/180, inc=((16.7 +  + initDeviation[6])*pi)/180, M=(154.8 + initDeviation[8])*pi/180, Omega=((295.5 + initDeviation[10])*pi)/180)
    sim2.add(primary=sim2.particles[0], a=a2 + initDeviation[1], m=10.19/1047.56, e=0.274 + initDeviation[3], omega=((240.8  + initDeviation[5])*pi)/180, inc=((13.5 +  + initDeviation[7])*pi)/180, M=(82.5 + initDeviation[9])*pi/180, Omega=((115 + initDeviation[11])*pi)/180)
    sim2.move_to_com()
    #sim2.additional_forces = helicityWrapper
    print("All systems go for dt: {} and a2: {}!".format(dt, a2))
    m = 0
    global t_x
    t_x = np.linspace(0, 10**endTime, 5000)
    global tSave
    tSave = [0]
    print("hot potato")
    global phaseSpaceParameters
    phaseSpaceParameters = ['a', 'e', 'omega', 'inc', 'M', 'Omega']
    for t in t_x:
        global tCurrent
        tCurrent = t
        try:
            sim1.integrate(t)
            sim2.integrate(t)
            helicityWrapper(sim1)
        except KeyboardInterrupt:
            print(" ")
            quit()
        else:
            if m%10 == 0:
                print("t: {:.2e} Completed - Sim_mag: {:.1e};{:.1e}".format(t, magnitude(sim1), magnitude(sim2)), end="\r")
                sys.stdout.flush()
            m += 1
    print("Cold potato")
    with open(dir_path + "HelicityTest - {}.txt".format(*initDeviation), 'w') as fle:
        toSave = (tSave, currentDeviation, derivativeDeviation, helicityAngles, twistAngles)
        for i in zip(*toSave):
            baseLine = '; '
            line = ["{}".format(j) for j in i]
            fle.write(baseLine.join(line) + '\n')
    tcomp = tm.perf_counter() - t0
    print("Completed in {:.2f} minutes".format(tcomp/60) + " for system a2: {}!".format(a2))

def worker(arguments):
    main(**arguments)


if __name__ == "__main__":
    worker({'a2': 2.53, 'initDeviation': tuple([0.01] + [0] * 11), 'endTime': 6})

