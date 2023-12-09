import numpy as np
import math
from scipy import interpolate
from scipy.optimize import fsolve, root_scalar, root
from typing import List
import pickle
import experiments_data
from const import simulation_params as sim_params
from const import const_local_mesh
from const import const_jintegral as j
from utils.logger import logger

nterm = experiments_data.nterm

def solve_z(zpos: float, th: float, xpos: float, a: List[float], h :float) -> float:
    # h = Decimal(str(h))
    # xpos = Decimal(str(xpos))
    # zpos = Decimal(str(zpos))
    # th = Decimal(str(th))
    if nterm == 5:
        def f(z):
            return (a[0]*z**4 + a[1]*z**3 + a[2]*z**2 + a[3]*z + a[4] - xpos)**2 + (z - zpos)**2 - h**2
    elif nterm == 4:
        def f(z):
            return (a[0]*z**4 + a[1]*z**3 + a[2]*z**2 + a[3] - xpos)**2 + (z - zpos)**2 - h**2
    elif nterm == 3:
        def f(z):
            return (a[0]*z**2 + a[1]*z + a[2] - xpos)**2 + (z - zpos)**2 - h**2
    elif nterm == 2:
        def f(z):
            # a[0] = Decimal(str(a[0]))
            # a[1] = Decimal(str(a[1]))
            # z = Decimal(str(z))
            return (a[0]*z**2 + a[1] - xpos)**2 + (z - zpos)**2 - h**2
    else:
        logger.error("--nterm is 2,3,4,5")
    sol = root_scalar(f, bracket=[zpos, th], method="bisect")
    # solb = get_bisect(f, zpos, th)
    # logger.info(f"sol: {sol.root}, solb: {solb}")
    return sol

def solve_zb(a:list, b:list, p:list):
    if nterm == 5 or nterm == 4:
        def f(z):
            return ((p[0] - z) * (4 * a[0] * p[0]**3 + 3 * a[1] * p[0]**2 + 2 * a[2] * p[0]) + 
                    (b[0] * z**4 + b[1] * z**3 + b[2] * z**2 + b[3] * z) - p[1]) + 1
    elif nterm == 3:
        def f(z):
            return ((p[0] - z) * (2 * a[0] * p[0] + a[1]) + (b[0] * z**2 + b[1] * z + b[2]) - p[1]) + 1
    elif nterm == 2:
        def f(z):
            return (p[1] - (b[0] * z**2 + b[1] )) * (2* a[0] * p[0]) + (p[0] - z)
    else:
        logger.error("nterm is 2,3,4,5")
    try:
        sol = root_scalar(f, bracket=[0, p[0]], method="bisect")
    except:
        logger.error('solver didnot work.')
        logger.error(f'a={a}, b={b}, p={p}')

    return sol

def solve_zv(a, b, p):
    if nterm == 5:
        def f(z):
            return (1/(p[0] - z)) * (p[1] - (b[0]*z**4 + b[1]*z**3 + b[2]*z**2 + b[3]*z + b[4])) * (4*a[0]*p[0]**3 + 3*a[1]*p[0]**2 + 2*a[2]*p[0] + a[3]) + 1
    elif nterm == 4:
        def f(z):
            return ((p[1] - (b[0]*z**4 + b[1]*z**3 + b[2]*z**2 + b[3]*z))/(p[0] - z)) * (4*a[0]*p[0]**3 + 3*a[1]*p[0]**2 + 2*a[2]*p[0]) + 1
    elif nterm == 3:
        def f(z):
            return ((p[1] - (b[0]*z**2 + b[1]*z + b[2]))/(p[0] - z)) * (2*a[0]*p[0] + a[1]) + 1
    elif nterm == 2:
        def f(z):
            return (p[1] - (b[0]*z**2+b[1]))*(2*a[0]*p[0]) + p[0]-z
    else:
        logger.error("nterm is 2,3,4,5")
    sol = root_scalar(f, bracket=[p[0], 100], method="bisect")
    return sol

def solve_0(a, b, p):
    if nterm == 5:
        def f(z):
            return ((p[1] + 35.)/(p[0] - z)) * (4*a[0]*p[0]**3 + 3*a[1]*p[0]**2 + 2*a[2]*p[0] + a[3]) + 1
    elif nterm == 4:
        def f(z):
            return ((p[1] + 35.)/(p[0] - z)) * (4*a[0]*p[0]**3 + 3*a[1]*p[0]**2 + 2*a[2]*p[0]) + 1
    elif nterm == 3:
        def f(z):
            return ((p[1] + 35.)/(p[0] - z)) * (2*a[0]*p[0] + a[1]) + 1
    elif nterm == 2:
        def f(z):
            return (((p[1] - z) * (2 * a[0] * p[1])) / (p[1] + 35)) + 1
    else:
        logger.error("nterm is 2,3,4,5")
    sol = root_scalar(f, bracket=[0, p[0]], method="bisect")
    return sol

def solve_zi(a, p):
    if nterm == 5:
        def f(z):
            return np.sqrt((p[0] - z)**2 + (p[1] - (a[0] * z**4 + a[1] * z**3 + a[2] * z**2 + a[3] * z + a[4]))**2)
    elif nterm == 4:
        def f(z):
            return ((p[1] - z) * (4 * a[0] * z**3 + 3 * a[1] * z**2 + 2 * a[2] * z) 
                    + (p[1] * (a[0] * z**4 + a[1] * z**3 + a[2] * z**2 + a[3] * z) - p[1] * p[2]))**2
    elif nterm == 3:
        def f(z):
            return ((p[1] - z) * (2 * a[0] * z + a[1]) 
                    + (p[1] * (a[0] * z**2 + a[1] * z + a[2]) - p[1] * p[2]))**2
    elif nterm == 2:
        def f(z):
            return ((p[1]-(a[0]*z**2+a[1])))*(2*a[0]*z) + (p[0]-z)
    else:
        logger.error("nterm is 2,3,4,5")
    sol = root_scalar(f, bracket=[-100, 100], method="bisect")
    return sol

def solve_zcf2(p, co, elesizeL, disinf):
    def f(z):
        return disinf(z, p[1] + (p[0] - z)/(2*co[0]*p[0]))[0] - elesizeL
    # logger.info(root(f, x0=p[0]).x)
    z = root(f, x0=p[0]).x[0]
    return z

def solve_zcf22(pz, elesizeL, disinf, postip2x):
    def f(x):
        return disinf([pz], x)[0] - elesizeL
    x = root(f, x0=postip2x).x[0]
    return x

def fa(a: List[float], z: float) -> float:
    # a[0] = Decimal(str(a[0]))
    # a[1] = Decimal(str(a[1]))
    # z = Decimal(str(z))
    if nterm == 5:
        return a[0]*z**4 + a[1]*z**3 + a[2]*z**2 + a[3]*z + a[4]
    elif nterm == 4:
        return a[0]*z**4 + a[1]*z**3 + a[2]*z**2 + a[3]
    elif nterm == 3:
        return a[0]*z**2 + a[1]*z + a[2]
    elif nterm == 2:
        return a[0]*z**2 + a[1]
    else:
        logger.error("nterm is 2,3,4,5")

def write_data(path, data):
    with open(path, "w") as f:
        f.write(data)

def coeff(t):
    nterm = experiments_data.nterm
    coefcoint = experiments_data.coefcoint
    return list([float(coefcoint[i](t)) for i in range(nterm)])

def levelset(a, pos, tip_x):
    if math.isclose(pos[0], 0.):
        sol = abs(pos[1] - tip_x)
    else:
        z = solve_zi(a, pos).root
        sol = ((pos[0] - z)**2 + (pos[1] - fa(a, z))**2)**0.5
        # logger.info(f"\nz: {z}, \npos: {pos}, \nsol: {sol}")
    if pos[1] < fa(a, pos[0]):
        sol = -sol
    return sol

def morethan(co:List[float]) -> List[float]:
    th = 15.
    nco = len(co)
    for i in range(1, nco):
        if co[i][0] >= th:
            nc = i+1
            break
    return co[:nc+1]

class LocalMesh:
    def __init__(self, step):
        self.step = step
        self.node = []
        self.num_node = len(self.node)
        self.elem = []
        self.num_elem = len(self.elem)

        self.tip_x = -35. + const_local_mesh.elesizeL*step
        coefhis = experiments_data.coefhis.tolist()
        nfcoef = experiments_data.nfcoef
        self.nterm = experiments_data.nterm

        if coefhis[nfcoef-1][self.nterm-1] <= self.tip_x:
            if self.tip_x <= coefhis[-1][self.nterm-1]:
                self.cotip = coeff(self.tip_x)
            else:
                self.cotip = np.concatenate((coefhis[-1][:self.nterm-1], [self.tip_x]))
        else:
            self.cotip = np.concatenate((coefhis[nfcoef-1][:self.nterm-1], [self.tip_x]))

        self.postip1x = const_local_mesh.elesizeL * (self.step + 1) - 35.
        if coefhis[nfcoef-1][self.nterm-1] <= self.postip1x:
            if self.postip1x <= coefhis[-1][self.nterm-1]:
                self.cotip1 = coeff(self.postip1x)
            else:
                self.cotip1 = np.concatenate((coefhis[-1][:self.nterm-1], [self.postip1x]))
        else:
            self.cotip1 = np.concatenate((coefhis[nfcoef-1][:self.nterm-1], [self.postip1x]))

    def _calc_crack_front_1_node(self):
        z = 0.0
        x = self.tip_x
        self.crack_front_1_node = [[z, x]]

        checker = 1
        elesizezL = const_local_mesh.elesizezL

        while checker < 5:
            th = 20
            sol = solve_z(z, th, x, self.cotip, elesizezL)
            z = sol.root
            x = fa(self.cotip, z)
            self.crack_front_1_node.append([z, x])
            if z >= 15: 
                checker += 1

    def _calc_back_node(self):
        elesizeL = const_local_mesh.elesizeL
        self.back_node_list = []
        for i in range(1,const_local_mesh.nback+1):
            if i == 1:
                _back_node = np.vstack([np.array([0.0, self.tip_x - i*elesizeL]), self.crack_front_1_node[1:]])
            elif i > 1:
                _back_node = np.vstack([np.array([0.0, self.tip_x - i*elesizeL]), _back_node[1:]])
            cotipb = np.hstack([self.cotip[:self.nterm-1], self.tip_x - i*elesizeL])
            checker = 1
            idx = 0
            th = 20
            while checker < 5:
                idx += 1
                if idx < len(_back_node):
                    z = solve_zb(self.cotip, cotipb, _back_node[idx]).root
                    x = fa(cotipb, z)
                    _back_node[idx][0] = z
                    _back_node[idx][1] = x
                else:
                    hdis = np.linalg.norm(np.array([_back_node[idx-1][0], _back_node[idx-1][1]]) - np.array([_back_node[idx-2][0], _back_node[idx-2][1]]))
                    z = solve_z(_back_node[idx-1][0], th, _back_node[idx-1][1], cotipb, hdis).root
                    x = fa(cotipb, z)
                    _back_node = np.vstack([_back_node, [z, x]])
                if z >= 15:
                    checker += 1
            self.back_node_list.insert(0, _back_node)

    def _calc_crack_front_2_node(self):
        z = 0.
        x = self.postip1x
        self.crack_front_2_node = []
        
        checker = 1
        idx_crack_front_1 = 0

        _crack_front_1_node = self.crack_front_1_node
        while checker < 5:
            th = 20
            if idx_crack_front_1 <= len(_crack_front_1_node):
                _z = _crack_front_1_node[idx_crack_front_1][0]
                _x = _crack_front_1_node[idx_crack_front_1][1]
                z = solve_zv(self.cotip, self.cotip1, [_z, _x]).root
                x = fa(self.cotip1, z)
            else:
                z = solve_z(z, th, x, self.cotip1, const_local_mesh.elesizeL).root
                x = fa(self.cotip1, z)
            idx_crack_front_1 += 1
            if z >= 0:
                self.crack_front_2_node.append([z, x])
            if z >= 15:
                checker += 1

    def _calc_crack_front_3_node(self):
        verin = const_local_mesh.verin
        horin = const_local_mesh.horin
        intin = const_local_mesh.intin
        verinmin = self.postip1x - verin
        horp = int(np.ceil(horin / intin)) + 1
        verp = int(np.ceil((2 * verin) / intin)) + 1
        NP = horp * verp

        posI = np.array([[intin * i, verinmin + intin * j] for j in range(-verp, verp+1) for i in range(-horp, horp+1)])
        x_inp = posI[:, 0]
        y_inp = posI[:, 1]
        z_inp = [levelset(self.cotip1, posI[i], self.postip1x) for i in range(len(posI))]

        self.disinf = interpolate.Rbf(x_inp, y_inp, z_inp, kind='cubic')

        z = 0.
        postip2x = const_local_mesh.elesizeL*(self.step + 2) - 35.
        self.postip2x = postip2x
        x = postip2x
        self.crack_front_3_node = np.array([[z, x]] + self.crack_front_2_node[1:])
        n = 1
        
        while n < len(self.crack_front_2_node):
            solcf2 = solve_zcf2(self.crack_front_2_node[n], self.cotip1, const_local_mesh.elesizeL, self.disinf)
            self.crack_front_3_node[n][0] = solcf2
            self.crack_front_3_node[n][1] = solve_zcf22(solcf2, const_local_mesh.elesizeL, self.disinf, postip2x)
            # logger.info(f"n: {n}, {self.crack_front_3_node[n]}")
            n += 1

    def _calc_front_node(self):
        self.front_node_list = []

        elesizeL = const_local_mesh.elesizeL
        nfront = const_local_mesh.nfront
        nlc = const_local_mesh.nlc
        delx = const_local_mesh.delx
        delz = const_local_mesh.delz
        lc = 1/nlc
        lco = elesizeL*lc

        for i in range(3, nfront+1):
            if i == 3: _front_node = np.vstack([[0, self.tip_x + i*elesizeL], self.crack_front_3_node[1:]])
            else: _front_node = np.vstack([[0, self.tip_x + i*elesizeL], _front_node[1:]])
            for idx_front_node in range(1, len(_front_node)):
                for k in range(nlc):
                    vec = np.array([(1/(2*delz)) * (self.disinf(_front_node[idx_front_node][0] + delz, _front_node[idx_front_node][1]) - self.disinf(_front_node[idx_front_node][0] - delz, _front_node[idx_front_node][1])),
                                    (1/(2*delx)) * (self.disinf(_front_node[idx_front_node][0], _front_node[idx_front_node][1] + delx) - self.disinf(_front_node[idx_front_node][0], _front_node[idx_front_node][1] - delx))])
                    _front_node[idx_front_node] += lco*(vec/np.linalg.norm(vec)).flatten()
            self.front_node_list.append(_front_node)

    def make_local_mesh(self):
        username = sim_params.USER_NAME
        reponame = sim_params.REPO_NAME
        dirname = sim_params.UBUNTU_DIR

        self._calc_crack_front_1_node()
        self._calc_back_node()
        self._calc_crack_front_2_node()
        self._calc_crack_front_3_node()
        self._calc_front_node()
        

        # count up nodes by line
        node_by_line = []
        num_node_by_line = []
        for line in self.back_node_list:
            node_by_line.append(morethan(line))
        node_by_line.append(morethan(self.crack_front_1_node))
        node_by_line.append(morethan(self.crack_front_2_node))
        node_by_line.append(morethan(self.crack_front_3_node))
        for line in self.front_node_list:
            node_by_line.append(morethan(line))
            
        for line in node_by_line:
            num_node_by_line.append(len(line))

        if num_node_by_line[0] != num_node_by_line[1]:
            node_by_line[0] = node_by_line[0][:num_node_by_line[1]]
            num_node_by_line[0] = len(node_by_line[0])

        self.node_by_line = node_by_line # posLzx0
        self.num_node_by_line = num_node_by_line # nposLzx0

        logger.info(f"self.cotip: {self.cotip}, self.cotip1: {self.cotip1}")
        
        self.node_zx = np.array([node for nodes in node_by_line for node in nodes])
        self.num_node_zx = len(self.node_zx)

        enLzx = []
        nback = const_local_mesh.nback
        nfront = const_local_mesh.nfront

        for j in range(1, nback+nfront+1):
            for i in range(1, num_node_by_line[j]):
                enLzx.append([
                            sum(num_node_by_line[:j-1])+i,
                            sum(num_node_by_line[:j-1])+i+1,
                            sum(num_node_by_line[:j-1])+num_node_by_line[j-1]+i+1, 
                            sum(num_node_by_line[:j-1])+num_node_by_line[j-1]+i
                            ])
        self.enLzx = enLzx
        num_elem_zx = len(enLzx)

        elesizeyL = const_local_mesh.elesizeyL
        ndyL = const_local_mesh.ndyL

        self.node = [[self.node_zx[idx][1]*10**-3, elesizeyL*j*10**-3, self.node_zx[idx][0]*10**-3] for j in range(ndyL+1) for idx in range(self.num_node_zx)]
        self.num_node = len(self.node)
        self.elem = [np.hstack([np.array(enLzx[idx])+self.num_node_zx*j, np.array(enLzx[idx])+self.num_node_zx*(j+1)]) for j in range(ndyL) for idx in range(num_elem_zx)]
        self.num_elem= len(self.elem)
        self.node_dat = [[idx+1] + list(self.node[idx]) for idx in range(self.num_node)]

    def generate(self):
        node = self.node
        num_node = self.num_node
        elem = self.elem
        num_elem = self.num_elem
        node_by_line = self.node_by_line

        node_dat = [[idx+1] + list(node[idx]) for idx in range(num_node)]
        elem_dat = [[idx+1] + list(elem[idx]) for idx in range(num_elem)]

        with open("node.l.dat", mode="w") as f:
            f.write(f"{num_node}\n")
            for line in node_dat:
                _ = []
                for l in line[1:]:
                    if l == 0.:
                        _.append("0.")
                    else:
                        _.append(float(l))
                f.write(f"{int(line[0])}\t{_[0]}\t{_[1]}\t{_[2]}\n")
            f.close()

        with open("elem.l.dat", mode="w") as f:
            f.write(f"{num_elem}\t8")
            tab = "\t"
            for line in elem_dat:
                f.write(f"\n{tab.join(list(map(str, map(int, line))))}")
        f.close()

        if sim_params.LOCAL_01_LIST[0] == 1:
            qcal = [[j.Rj0], [j.Rj1], [j.Wj0], [j.Wj1]]
            lL = const_local_mesh.lL
            behinddat = [len([j for i in node_by_line[:lL-1] for j in i])]
            np.savetxt("qcalculation.dat", qcal, fmt="%s")
            np.savetxt("behind.dat", behinddat, fmt="%s")

        with open("local_mesh.pickel", mode="wb") as f:
            pickle.dump(self, f)
