import math
from const import simulation_params as sim_params

NOR_LIST = [0.8, 0.8, 1.0, 1.0, 1.2 * math.sqrt(0.99), 1.2 * math.sqrt(0.99)]
nor = NOR_LIST[sim_params.TEST_NUMBER-1]

mag = 1         # magnification
bx1 = 1.0
bx2 = 1.5
by1 = 1.5
r1 = 8.0 - nor  # [mm]
nGr1 = 2 * mag
nGtheta1 = 4 * mag
r1len = r1      # [mm]
r1Mesh = nGr1
theta1Mesh = nGtheta1

hGyMin = (r1len/r1Mesh + nor) * math.tan(math.pi / 2 / nGtheta1)

nGx1 = math.floor(8 / hGyMin)
bx1 = bx1
nGx2 = math.floor((22 - nor) / hGyMin)
bx2 = bx2
nGx3 = math.floor((150 - r1len - r1) / hGyMin)

nGy1 = 10 * mag
by1 = by1
hGMin = min(8.0 / nGx1, (22 - nor) / nGx2, r1len / r1Mesh / nGr1, (150 - r1len / r1Mesh) / nGx3)
hGMax = max(8.0 / nGx1, (22 - nor) / nGx2, r1len / r1Mesh / nGr1, (150 - r1len / r1Mesh) / nGx3)
nGz1 = int(15 / hGMin)

x1Mesh = nGx1
x1r = bx1
x2Mesh = nGx2
x2r = bx2
x4Mesh = nGx3
z1Mesh = nGz1

y1Mesh = nGy1
y1r = by1
r1 = r1len / r1Mesh

x9Mesh = 1
x10Mesh = 1
x11Mesh = 1
