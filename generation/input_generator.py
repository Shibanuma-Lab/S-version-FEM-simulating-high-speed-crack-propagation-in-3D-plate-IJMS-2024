import numpy as np
from const import const_jintegral as j
from const import simulation_params as sim_params
from const import const_local_mesh
from utils.logger import logger

def generate(step) -> None:
    REstart = sim_params.REstart
    if step == 0:
        REstart = 0
    # input data の作成
    input_data = [       
            [2, "\t!>solutiontype(1:static 2:dynamic)"],
            [0, "\t!>isNLgeom(0:off 1:on)"],
            [5, "\t!> max NR step"],
            [sim_params.STEP_TIME, "\t!> dt"],
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
            [const_local_mesh.href, "\t!> local mesh h-refine"],
            [REstart, "\t!> is_Restart (0:off, 1:on)"],
            [2, "\t!> XFEM(1:global, 2:local)"],
            [0.015, "\t!> thickness"],
            [sim_params.OPENMP, "\t!> the number of OpenMP threads"]
        ]
    if step == (step if sim_params.DYNAMIC_01_LIST[0] == 0 else sim_params.INTERM) and REstart == 0:
        logger.info("step: {} :: input.dat static".format(step))
        input_data[15][0] = 0
        input_data[0][0] = 1
    np.savetxt("input.dat", input_data, fmt="%s")
