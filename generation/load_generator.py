import numpy as np

def generate():
    loadGdat = [0]
    np.savetxt("load.dat", loadGdat, fmt="%s")