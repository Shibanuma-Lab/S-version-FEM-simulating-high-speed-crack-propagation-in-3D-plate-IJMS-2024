import numpy as np
from typing import List

def theta(cfront: List, cfront1: List):
    """ input: cfront, cfront1
        output: theta(angle of cfront and cfront1)
    """
    if len(cfront1) < len(cfront):
        cfront00 = cfront[1:-1]
        cfront01 = cfront1
    elif len(cfront1) > len(cfront):
        cfront00 = cfront
        cfront01 = cfront1[1:-1]
    elif len(cfront1) == len(cfront):
        cfront00 = cfront
        cfront01 = cfront1
    l = min(len(cfront00), len(cfront01))
    theta = [0] * l
    for i in range(l):
        theta[i] = np.arctan((cfront01[i][0] - cfront00[i][0]) / (cfront01[i][1] - cfront00[i][1]))
    return theta
