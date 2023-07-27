import math
from const import const_global_mesh

hL = 0.3      # loal mesh size in the x direction[mm]
elesizeL = hL
rxzL = 2      # ratio of x - z
rxyL = 1        # ratio of x - y
hLz = rxzL * hL  # local mesh size in the z direction[mm]
hLy = rxyL * hL # local mesh size in the y direction[mm]

Localmag = 1
HL = math.ceil(1.2 * const_global_mesh.hGyMin / hLy * Localmag)         # num.of elements above the crack tip
aL = math.ceil(2.83 * const_global_mesh.hGMax / elesizeL * Localmag)    # num.of elements behind the crack tip
lL = 15                                                                 # (Localmag)

Wint = 30.      # width of interpolation
Hint = 30.      # height of interpolation
Inint = 3.0     # interval of interpolation

delz = 0.00001  # delta z
delx = 0.00001  # delta x
nlc = 1         # num.of iteration when coordinates of the local nodes are determined

elesizezL = hLz
elesizeyL = hLy
nback = aL
nfront = lL
ndyL = HL
verin = Wint
intin = Inint
horin = Hint
href = 3
