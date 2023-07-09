import os
import subprocess
from const import simulation_params as sim_params
from utils.step2str import step2str
from utils.logger import logger

def run(step):
    logger.info(f"step: {step} :: Run Fortrun...")
    username = sim_params.USER_NAME
    reponame = sim_params.REPO_NAME
    dirname = sim_params.UBUNTU_DIR
    dirnametest = sim_params.DIR_NAME_TEST
    day = sim_params.DAY

    str_step = step2str(step)
    logger.info(f"run: now {os.getcwd()}")
    try:
        os.mkdir(f"Newton/{dirnametest}")
    except:
        logger.error("Did NOT make directory")
    try:
        os.mkdir(f"Newton/{dirnametest}/{day}")
    except:
        logger.error("Did NOT make directory")
    try:
        os.mkdir(f"Newton/{dirnametest}/{day}/step{str_step}")
    except:
        logger.error("Did NOT make directory")
    
    subprocess.run([f"cp", "-r", f"inputfiles/step{str_step}",  f"../{reponame}/example"])
    subprocess.run([f"chmod", "+x", f"../{reponame}/bin/sfem_linear"])
    subprocess.run([f"chmod", "+r", f"../{reponame}/bin/sfem_linear"])
    os.chdir(f"../{reponame}/example/step{str_step}")
    subprocess.run([f"../../bin/sfem_linear"])
    os.chdir(f"../") # /example
    subprocess.run([f"cp", f"step{str_step}", f"../../generation/Newton/{dirnametest}/{day}", "-r"])
    subprocess.run([f"rm", f"step{str_step}", "-r"])
    os.chdir(f"../../generation")
    logger.info(f"run: end {os.getcwd()}")
