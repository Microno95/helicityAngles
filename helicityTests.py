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


def main(a2=2.53, end_time=6, part=1, sep=1.1, dt=0.01, a1=0.861, recVar=False):
    pi = np.pi
    dir_path = "./Variational Lyapunov/"
    ensure_dir(dir_path)
    t0 = tm.perf_counter()
    sim = rb.Simulation()
    sim.units = ('Yr', 'AU', 'Msun')
    sim.integrator = 'ias15'
    sim.dt = dt
    sim.add(m=1.31, x=0, y=0, z=0)
    sim.add(primary=sim.particles[0], a=a1, m=14.57/1047.56, e=0.239, omega=(290*pi)/180, inc=(16.7*pi)/180, M=154.8*pi/180, Omega=(295.5*pi)/180)
    sim.add(primary=sim.particles[0], a=a2, m=10.19/1047.56, e=0.274, omega=(240.8*pi)/180, inc=(13.5*pi)/180, M=82.5*pi/180, Omega=(115*pi)/180)
    if not any(fle.startswith("SemiSys_Lyap_a_{:.5e}_P_{:.0f}_{:.1e}_Yrs_dt".format(a2, part, 10**end_time)) for fle in os.listdir(dir_path)): 
        print("All systems go for dt: {} and a2: {}!".format(dt, a2))
        sim.move_to_com()
        m = 0
        t_x = np.linspace(0, 10**end_time, 10)
        print("hot potato")
        for t in t_x:
            sim.integrate(t)
            try:
                sim_mag.append(magnitude(sim))
            except OverflowError:
                break
            else:
                if m%100 == 0:
                    print("t: {:.2e} Completed - Sim_mag: {:.1e}".format(t, magnitude(sim)))
                m += 1
        print("Cold potato")
        with open(dir_path + "SemiSys_Lyap_a_{:.5e}_P_{:.0f}_{:.1e}_Yrs_dt_{:.2e}.txt".format(a2, part, 10**end_time, sim.dt), 'w') as fle:
            toSave = (t_x, sim_mag, e1, e2, a_1, a_2)
            for i in zip(*toSave)
                fle.write(str.join(["{};".format(k) for k in i]) + '\n'))
        tcomp = tm.perf_counter() - t0
        print("Completed in {:.2f} minutes".format(tcomp/60) + " for system a2: {}!".format(a2))
    else:
        print("File already exists, passing...")

def worker(arguments):
    main(**arguments)


if __name__ == "__main__":
    pool = mp.Pool(processes=1)
    argu = []
    try:    
        dir_path = "./Variational Lyapunov/"
        argu.append({'a2': 2.53, 'deviation': (1,1,)})
        pool.map_async(worker, argu)
        pool.close()
    except:
        pool.terminate()
        quit()
    finally:
        pool.join()

