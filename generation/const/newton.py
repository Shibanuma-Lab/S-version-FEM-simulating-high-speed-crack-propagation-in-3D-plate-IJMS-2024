import datetime

# Newton loop
day = datetime.datetime.now().date().isoformat() # get current date
zmax = 30/2*0.001
lenz = 20
calstep = 1
calsteppre = 1
error = 0
eps = 0.0005
dis = 0.00001
differential = 0 # 0 for forward, 1 for central
BFGS = 0