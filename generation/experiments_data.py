import numpy as np
import pandas as pd
import os
from scipy.interpolate import interp1d
from const import simulation_params as sim_params

with open(f"experiments_data/test{sim_params.TEST_NUMBER}.inp", "r") as f:
    data = f.readlines()
data = [row.split()[0] for row in data]
radius = float(data[0])
nfcoef = int(data[1])
coefhis = pd.read_excel(f"experiments_data/{data[2]}", header=None, index_col=None).to_numpy()
nimg = coefhis.shape[0]
nterm = coefhis.shape[1]
coeft = [[(coefhis[i, nterm-1], coefhis[:, j][i]) for i in range(nfcoef-1, nimg)] for j in range(nterm)]
coefcoint = [interp1d([row[0] for row in coeft[i]], [row[1] for row in coeft[i]], kind="cubic", fill_value='extrapolate') for i in range(nterm)]
velhis = pd.read_excel(f"experiments_data/{data[3]}",  header=None, index_col=None).to_numpy()
velcoint = interp1d(velhis[:,0], velhis[:, 1], kind="cubic", fill_value='extrapolate')
ctod = float(data[4]) * 0.001
cross01 = int(data[5])
hole = float(data[6])
