import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import numpy as np
import jintegral as j
import local_mesh

from utils.step2str import step2str
from const import simulation_params as sim_params
from const import const_local_mesh
import experiments_data

stepini = 3
steplast = 300
x, y = [], []
# str_step = step2str(step)
dirnametest = sim_params.DIR_NAME_TEST
day = sim_params.DAY
for step in range(stepini, steplast):
    str_step = step2str(step)
    path = f"/Newton/{dirnametest}/{day}/step{str_step}/Jlist.dat"
    with open(path, "r") as f:
        lines = f.readlines()
        for line in lines:
            if line == "None\n":
                continue
            j = float(line)
            postipx = const_local_mesh.elesizeL * step - 35
            velhis = experiments_data.velhis
            velcoint = experiments_data.velcoint
            REstart = sim_params.REstart
            Local01 = const_local_mesh.Local01
            if step == 0:
                REstart = 0
            # input data の作成
            v = velhis[0][1] if velhis[0][0] > postipx else velcoint(-35+const_local_mesh.elesizeL*(step-1))
            x.append(v)
            y.append(j)

plt.xlabel("step")
plt.ylabel("J")
plt.legend()
plt.scatter(x, y, s = 0.1)
plt.savefig("J.png")