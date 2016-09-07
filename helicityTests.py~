import time as tm
import os, os.path
import shutil as shl
import numpy as np
import matplotlib.pyplot as plt
import rebound as rb
import sys as sys
import multiprocessing as mp


def ensure_dir(path): # Ensures that the path to a directory exists, and creates it if not
    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path): # Raises the OSError if and only if the specified path is a file and not a directory
            raise


def magnitude(vari):
    length_vari = 0
    for i in vari.particles[:3]:
        length_vari += i.x**2 + i.y**2 + i.z**2
        length_vari += i.vx**2 + i.vy**2 + i.vz**2
    return np.sqrt(length_vari)


def main(a2=2.53, a1=0.861, endTime=6, part=1, sep=1.1, dt=0.01, recVar=False, initDeviation=(0, 0, )):
    currentDeviation = [list(initDeviation[:])]
    derivativeDeviation = [[0] * len(initDeviation)]
    pi = np.pi
    dir_path = "./Variational Lyapunov/"
    ensure_dir(dir_path)
    t0 = tm.perf_counter()
    sim1 = rb.Simulation()
    sim1.units = ('Yr', 'AU', 'Msun')
    sim1.integrator = 'ias15'
    sim1.dt = dt
    sim1.add(m=1.31, x=0, y=0, z=0)
    sim1.add(primary=sim1.particles[0], a=a1, m=14.57/1047.56, e=0.239, omega=(290*pi)/180, inc=(16.7*pi)/180, M=154.8*pi/180, Omega=(295.5*pi)/180)
    sim1.add(primary=sim1.particles[0], a=a2, m=10.19/1047.56, e=0.274, omega=(240.8*pi)/180, inc=(13.5*pi)/180, M=82.5*pi/180, Omega=(115*pi)/180)
    sim1.move_to_com()
    sim2 = rb.Simulation()
    sim2.units = ('Yr', 'AU', 'Msun')
    sim2.integrator = 'ias15'
    sim2.dt = dt
    sim2.add(m=1.31, x=0, y=0, z=0)
    sim2.add(primary=sim2.particles[0], a=a1 + initDeviation[0], m=14.57/1047.56, e=0.239 + initDeviation[2], omega=((290 + initDeviation[4])*pi)/180, inc=((16.7 +  + initDeviation[6])*pi)/180, M=(154.8 + initDeviation[8])*pi/180, Omega=((295.5 + initDeviation[10])*pi)/180)
    sim2.add(primary=sim2.particles[0], a=a2 + initDeviation[1], m=10.19/1047.56, e=0.274 + initDeviation[3], omega=((240.8  + initDeviation[5])*pi)/180, inc=((13.5 +  + initDeviation[7])*pi)/180, M=(82.5 + initDeviation[9])*pi/180, Omega=((115 + initDeviation[11])*pi)/180)
    sim2.move_to_com()
    print("All systems go for dt: {} and a2: {}!".format(dt, a2))
    m = 0
    t_x = np.linspace(0, 10**2, 100)
    print("hot potato")
    tPrev = 0
    phaseSpaceParameters = ['a', 'e', 'omega', 'inc', 'M', 'Omega']
    for t in t_x:
        sim1.integrate(t)
        sim2.integrate(t)
        try:
            if tPrev > 0:
                calculatedOrbitsSim1 = sim1.calculate_orbits()
                calculatedOrbitsSim2 = sim2.calculate_orbits()
                currentDeviation.append([])
                derivativeDeviation.append([])
                for i in phaseSpaceParameters:
                    for j in range(2):
                        currentDeviation[-1].append(eval("calculatedOrbitsSim2[{}].{} - calculatedOrbitsSim1[{}].{}".format(j, i, j, i)))
                        derivativeDeviation[-1].append(currentDeviation[-1][-1] / (t - tPrev))
        except OverflowError:
            break
        else:
            if m%100 == 0:
                print("t: {:.2e} Completed - Sim_mag: {:.1e}".format(t, magnitude(sim1)))
            tPrev = t
            m += 1
    print("Cold potato")
    with open(dir_path + "HelicityTest - {}.txt".format(*initDeviation), 'w') as fle:
        toSave = (t_x, derivativeDeviation, currentDeviation)
        for i in zip(*toSave):
            baseLine = '; '
            line = ["{}".format(j) for j in i]
            fle.write(baseLine.join(line) + '\n')
    tcomp = tm.perf_counter() - t0
    print("Completed in {:.2f} minutes".format(tcomp/60) + " for system a2: {}!".format(a2))

def worker(arguments):
    main(**arguments)


if __name__ == "__main__":
    worker({'a2': 2.53, 'initDeviation': tuple([1] + [0] * 11), 'endTime': 2})

