import math
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from const import const_local_mesh, simulation_params as sim_params, const_jintegral
from utils import step2str
def Kv(v):
    return (5.99825e-4 * (v / 100) ** 10 + 0.0779521 * (v / 100) ** 2 + 1.6154) * 100

def getvJ(step):
    str_step = step2str.step2str(step)
    dirnametest = sim_params.DIR_NAME_TEST
    day = sim_params.DAY
    nodeldatans = []
    with open(f"Newton/{dirnametest}/{day}/step{str_step}/node.l.dat", "r") as f:
        lines = f.readlines()
        for line in lines[1:]:
            nodeldatans.append([float(i) for i in line.split()[1:]])
    with open(f"Newton/{dirnametest}/{day}/step{str_step}/Jlist.dat", "r") as f:
        lines = f.readlines()
        Jlist = []
        for line in lines:
            try:
                Jlist.append(float(line))
            except ValueError:
                pass
        Jlen = len(Jlist)
    with open(f"Newton/{dirnametest}/{day}/step{str_step}/input.dat", "r") as f:
        lines = f.readlines()
        dt = float(lines[3].split()[0])
        deltat = dt * const_jintegral.inc
        v = 2./deltat if sim_params.Local01 == 0 else (const_local_mesh.elesizeL*0.001/deltat if deltat>0 else 1.)

    elesizeL = const_local_mesh.elesizeL
    for i in range(len(nodeldatans)):
        if math.isclose(nodeldatans[i][0], (elesizeL*step-35.)*0.001):
            nnmlans = i
            break
    for i in range(len(nodeldatans)):
        if math.isclose(nodeldatans[i][0], (elesizeL*(step+1)-35.)*0.001):
            nnmlans2 = i
            break
    
    velfrontlist = [
        (((node[0] - node2[0])**2 + (node[2] - node2[2])**2)**0.5 / elesizeL * v * 1000)
        for node, node2 in zip(nodeldatans[nnmlans:nnmlans + Jlen], nodeldatans[nnmlans2:nnmlans2 + Jlen])
    ]

    def Kcal(Jd, v):
        ee = const_jintegral.ee
        nu = const_jintegral.Nu
        rho = const_jintegral.Rho

        v1 = ((1 - nu) * ee) / ((1 + nu) * (1 - 2 * nu) * rho) ** 0.5
        v2 = (ee / ((1 + nu) * 2 * rho)) ** 0.5
        beta1 = (1 - (v / v1) ** 2) ** 0.5
        beta2 = (1 - (v / v2) ** 2) ** 0.5
        AI = (beta1 * (1 - beta2 ** 2)) / (4 * beta1 * beta2 - (1 + beta2 ** 2) ** 2)
        out = ((ee * Jd) / ((1 + nu) * AI)) ** 0.5
        return out

    # Assuming zmax and velfrontlist are defined elsewhere
    K = [Kcal(J_val, v_val) / 10 ** 4 for J_val, v_val in zip(Jlist, velfrontlist)]
    vKcal = list(zip(velfrontlist, K))

    
    return [(v, Kv(v)) for v in velfrontlist]

def getdata():
    data1 = pd.read_excel("experiments_data/vel_J_1.xlsx", header=None, index_col=None).to_numpy()
    data2 = pd.read_excel("experiments_data/vel_J_2.xlsx", header=None, index_col=None).to_numpy()

    data1 = [(v, Kv(v)) for v in data1[:, 0]]
    data2 = [(v, Kv(v)) for v in data2[:, 0]]

    return data1, data2

if __name__ == "__main__":
    print(f"testname: {sim_params.DIR_NAME_TEST}")
    print(f"day: {sim_params.DAY}")
    for step in range(0, 300):
        try:
            l = getvJ(step)
        except FileNotFoundError:
            print(f"step{step} failed")
            continue
        if len(l) > 0:
            x, y = zip(*l)
        else:
            continue
        plt.scatter(x, y, c="blue", s=1)

    data1, data2 = getdata()
    plt.scatter([v[0] for v in data1], [v[1] for v in data1], c="red", s=1)
    plt.scatter([v[0] for v in data2], [v[1] for v in data2], c="yellow", s=1)

    plt.xlabel("v")
    plt.ylabel("K")
    plt.savefig("K.png")
