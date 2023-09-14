#! /usr/bin/python3

import os
import numpy as np
from multiprocessing import Pool
import shutil
import argparse
import time
from const import const_local_mesh, const_global_mesh, simulation_params as sim_params
from initial import initial
import global_mesh
import local_mesh
import boundary
import input_generator
import load_generator
import linux_command
from jintegral import jintegral
from utils.logger import logger
from utils.step2str import step2str

def makemeshs(step, REstart, INTERM):
    logger.info(os.getcwd())
    str_step = step2str(step)
    try:
        os.mkdir(f"inputfiles/step{str_step}")
    except:
        pass
    os.chdir(f"inputfiles/step{str_step}")
    logger.info(f"step: {step} :: Generate Global Mesh")
    g = global_mesh.GlobalMesh()
    g.generate(step)
    logger.info(f"step: {step} :: Generate Local Mesh")
    l = local_mesh.LocalMesh(step)
    l.make_local_mesh()
    l.generate()
    logger.info(f"step: {step} :: Generate Boundary")
    b = boundary.Boundary(l)
    logger.info(f"step: {step} :: Define Global Boundary")
    b.define_global_boundary()
    logger.info(f"step: {step} :: Define Local Boundary")
    b.define_local_boundary(l)
    logger.info(f"step: {step} :: Generate Boundary")
    b.generate()
    input_generator.generate(step, REstart, INTERM)
    load_generator.generate()
    if step >= 1 + INTERM or REstart == 1:
        init = initial(step, l, g)
    os.chdir("../../")
    return l, g

def jint(step, l):
    try:
        logger.info(f"step: {step} :: Jint: {jintegral(step, l)}")
    except Exception as e:
        logger.error(f"step: {step} :: {e}")

def main():
    logger.info("SIMULATION START")
    logger.info(f"hL: {const_local_mesh.hL}")
    logger.info(f"hLy: {const_local_mesh.hLy}")
    logger.info(f"aL, HL: {const_local_mesh.aL}, {const_local_mesh.HL}")
    logger.info(f"MinhG: {const_global_mesh.hGMin}")
    logger.info(f"MaxhG: {const_global_mesh.hGMax}")
    logger.info(f"MinhGy: {const_global_mesh.hGyMin}")
    logger.info(f"nGtheta1: {const_global_mesh.nGtheta1}")
    logger.info(f"hGz: {15. / const_global_mesh.nGz1}")
    logger.info(f"hG/HL > 6: {const_global_mesh.hGMin/const_local_mesh.hL}")
    logger.info(f"aL/hG>2.83: {const_local_mesh.aL*const_local_mesh.hL/const_global_mesh.hGMax}")
    logger.info(f"lL/hL>15: {const_local_mesh.lL}")
    logger.info(f"HL/hG>1.2: {const_local_mesh.HL*const_local_mesh.hLy/const_global_mesh.hGyMin}")
    logger.info(os.getcwd())

    argparser = argparse.ArgumentParser()
    argparser.add_argument("--test_start", type=int, default=sim_params.TEST_START, help="")
    argparser.add_argument("--test_end", type=int, default=sim_params.TEST_END, help="")
    argparser.add_argument("--step_start", type=int, default=0, help="")
    argparser.add_argument("--step_end", type=int, default=1000, help="")
    argparser.add_argument("--delete", type=bool, default=False, help="if Ture, delete all files in inputfiles and results. (default: False)")
    argparser.add_argument("--debugmode", type=bool, default=False, help="if True, run in debug mode. (default: False)")
    argparser.add_argument("--particular", type=bool, default=False, help="if True, run in particular step. (default: False)")
    argparser.add_argument("--is_jonly", type=bool, default=False, help="if True, only culculate jintegral. (default: False)")
    argparser.add_argument("--is_meshonly", type=bool, default=False, help="if True, only generate mesh. (default: False)")
    args = argparser.parse_args()
    logger.info(f"args: {args}")

    test_start = args.test_start
    test_end = args.test_end
    user_name = sim_params.USER_NAME

    step_start = args.step_start
    step_last = args.step_end
    INTERM = 0
    STEP_LAST = 300
    REstart = 0 if step_start == INTERM else 1  # 1:restart from dynamic analysis
    logger.info(f"REstart: {REstart}")

    for test in range(test_start, test_end + 1):
        if args.debugmode:
            if args.particular:
                step = int(input("---- step: "))
                l, g = makemeshs(step, REstart, INTERM)
                if args.is_jonly:
                    if step >= 3:
                        Jlist = jintegral(step, l)
                        logger.info(f"Jlist: {Jlist}")
                else:
                    if step >= 1 + INTERM or REstart == 1:
                        init = initial(step, l)
                    if args.is_meshonly:
                        break
                    linux_command.run(step)
                    if step >= 3:
                        Jlist = jintegral(step, l)
                        logger.info(f"Jlist: {Jlist}")
            else:
                step_list = range(step_start, step_last)
                for step in step_list:
                    l, g = makemeshs(step, REstart, INTERM)
                    if args.is_jonly:
                        if step >= 3:
                            Jlist = jintegral(step, l)
                            logger.info(f"Jlist: {Jlist}")
                        continue
                    if step >= 1 + INTERM or REstart == 1:
                        init = initial(step, l, g)
                    if args.is_meshonly:
                        continue
                    linux_command.run(step)
                    if step >= 3:
                        Jlist = jintegral(step, l)
                        logger.info(f"Jlist: {Jlist}")
        else:
            step_start = args.step_start
            step_last = args.step_end
            step_list = range(step_start, step_last)
            for step in step_list:
                l, g = makemeshs(step, REstart, INTERM)
                if args.is_jonly:
                    if step >= 3:
                        Jlist = jintegral(step, l)
                        logger.info(f"Jlist: {Jlist}")
                    continue
                if step >= 1 + INTERM or REstart == 1:
                    init = initial(step, l, g)
                if args.is_meshonly:
                    break
                linux_command.run(step)
                if step >= 3:
                    Jlist = jintegral(step, l)
                    logger.info(f"Jlist: {Jlist}")
                # logger.info(f"Jlist: {Jlist}")

if __name__ == "__main__":
    main()
