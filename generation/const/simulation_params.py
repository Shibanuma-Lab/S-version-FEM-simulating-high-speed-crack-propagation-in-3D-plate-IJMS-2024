import datetime

TEST_NUMBER_LIST = [4]
TEST_NUMBER = TEST_NUMBER_LIST[0]
INC_LIST = [1]
DIR_NAME_ADD_LIST = ["R0.8_No.1", "R0.8_No.2", "R1.0_No.1", "R1.0_No.2", "R1.2_No.1", "R1.2_No.2"]
DIR_NAME_ADD = DIR_NAME_ADD_LIST[TEST_NUMBER-1]
DIR_NAME_TEST_LIST = ["R0_8__NO1", "R0_8__NO2", "R1_0__NO1", "R1_0__NO2", "R1_2__NO1", "R1_2__NO2"]
DIR_NAME_TEST = DIR_NAME_TEST_LIST[TEST_NUMBER-1]
# DAY = datetime.datetime.now().strftime("%Y-%m-%d")
DAY = "test"

LOCAL_01_LIST = [1]
Local01 = LOCAL_01_LIST[0]
DYNAMIC_01_LIST = [1]
HREF_LIST = [3]
HREF = HREF_LIST[0]
NFR_LIST = [1]

TEST_START = 1
TEST_END = 1
GET_CTOD = 0

UBUNTU_DIR = "PMMA_old"

OPENMP = 8      # number of openmp thread
DOS_OPEN = 2    # dos screen 0: Not open　1:Open and close　2:Keeping open

ABO = 0  # 1:abort just before calculation of solver
USER_NAME = "lab"
REPO_NAME = "sfem_linear"

REstart = 1  # 1:restart from dynamic analysis
STEP_START = 0
INTERM = 0
STEP_LAST = 300
STEP_TIME = 1.2915155311516764e-6

OUTPUT_FOLDER = "outputfolder"
CALC_STEP = 1
JOB_NAME = "step_" + str(STEP_START) + "_exp_" + str(CALC_STEP)
FOLDER_NAME = JOB_NAME

is_cross = False
isgetctod = False
