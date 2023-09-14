import numpy as np
from const import const_jintegral as j
from const import simulation_params as sim_params
from const import const_local_mesh
from utils.logger import logger
import experiments_data

def generate(step, REstart, INTERM) -> None:
    postipx = const_local_mesh.elesizeL * step - 35
    velhis = experiments_data.velhis
    velcoint = experiments_data.velcoint
    Local01 = const_local_mesh.Local01
    if step == 0:
        REstart = 0
    # input data の作成
    v = velhis[0][1] if velhis[0][0] > postipx else velcoint(-35+const_local_mesh.elesizeL*(step-1))
    deltat = 2./v if Local01 == 0 else (const_local_mesh.elesizeL*0.001/v if v>0 else 1.)
    dt = deltat/j.inc
    input_data = [       
            [2, "\t!>solutiontype(1:static 2:dynamic)"],
            [0, "\t!>isNLgeom(0:off 1:on)"],
            [5, "\t!> max NR step"],
            [dt, "\t!> dt"],
            [1, "\t!>max time step"],
            [j.inc, "\t!> nincrement"],
            [j.ee, "\t!> young 率"],
            [j.Nu, "\t!> Poisson 比"],
            [j.Rho, "\t!> density"],
            [j.Gamma, "\t!> gamma for newmark-beta"],
            [j.Beta, "\t!> beta for newmark-beta"],
            [j.Alpha_l, "\t!> Rm for Rayleigh damping for mass"],
            [j.Beta_l, "\t!> Rk for Rayleigh damping for mass"],
            [2, "\t!> local mesh integral point (not used)"],
            [sim_params.HREF, "\t!> local mesh h-refine"],
            [1, "\t!> is_Restart (0:off, 1:on)"],
            [2, "\t!> XFEM(1:global, 2:local)"],
            [0.015, "\t!> thickness"],
            [sim_params.OPENMP, "\t!> the number of OpenMP threads"]
        ]
    if step == (step if sim_params.DYNAMIC_01_LIST[0] == 0 else INTERM) and REstart == 0:
        logger.info("step: {} :: input.dat static".format(step))
        input_data[15][0] = 0 # is_Restart
        input_data[0][0] = 1 # solutiontype
    np.savetxt("input.dat", input_data, fmt="%s")
