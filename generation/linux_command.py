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

    try:
        os.mkdir(f"/home/{username}/{reponame}/arrest/generation/Newton/{dirnametest}")
        os.mkdir(f"/home/{username}/{reponame}/arrest/generation/Newton/{dirnametest}/{day}")
        os.mkdir(f"/home/{username}/{reponame}/arrest/generation/Newton/{dirnametest}/{day}/step{str_step}")
    except:
        pass

    com_list = [
        f"cp -r {os.getcwd()} /home/{username}/{reponame}/example", # /home/{username}/sfem/example
        f"cd /home/{username}/{reponame}/example/step{str_step}",
        f"export OPM_NUM_THREADS={sim_params.OPENMP}",
        f"cd /home/{username}/{reponame}/example",
        f"chmod +x /home/{username}/{reponame}/bin/sfem_linear",
        f"/home/{username}/{reponame}/bin/sfem_linear",
        f"cp /home/{username}/{reponame}/example/step{str_step}/ /home/{username}/{reponame}/example/{dirname}/ -r",
        f"cp /home/{username}/{reponame}/example/step{str_step}/ /home/{username}/{reponame}/arrest/generation/Newton/{dirnametest}/{day} -r",
        f"rm /home/{username}/{reponame}/example/step{str_step} -r",
        f"cp /home/{username}/{reponame}/arrest/generation/inputfiles/step{str_step}/ /home/{username}/{reponame}/arrest/generation/Newton/{dirnametest}/{day} -r",

    ]

    for com in com_list:
        subprocess.run(com, shell=True)
        