import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import numpy as np
import jintegral as j
import local_mesh

from utils.step2str import step2str
from const import simulation_params as sim_params

stepini = 3
steplast = 300
x, y = [], []
# str_step = step2str(step)
dirnametest = sim_params.DIR_NAME_TEST
day = "20230622"
for step in range(stepini, steplast):
    str_step = step2str(step)
    path = f"/home/lab/sfem_linear/arrest/generation/Newton/R1_0__NO2/{day}/step{str_step}/Jlist.dat"
    with open(path, "r") as f:
        lines = f.readlines()
        for line in lines:
            if line == "None\n":
                continue
            j = float(line)
            # for j in Jlist:
            x.append(step)
            y.append(j)

plt.xlabel("step")
plt.ylabel("J")
plt.legend()
plt.scatter(x, y, s = 0.1)
plt.savefig("J.png")