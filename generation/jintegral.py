import numpy as np
import math
from typing import List

from const import const_local_mesh
from const import const_jintegral as jint
from const import simulation_params as sim_params
from utils.step2str import step2str
from utils.theta import theta
from utils.logger import logger

def q0(i: List[float]) -> float:
    """ input: i = [R, W]
        output: q0
    """
    qR = 0
    if i[0] < jint.Rj0:
        qR = 1
    elif jint.Rj0 <= i[0] < jint.Rj1:
        qR = (jint.Rj0 - i[0])/(jint.Rj1-jint.Rj0)
    elif jint.Rj1 <= i[0]:
        qR = 0
    
    qW = 0
    if i[1] < jint.Wj0:
        qW = 1
    elif jint.Wj0 <= i[1] < jint.Wj1:
        qW = (jint.Wj0 - i[1])/(jint.Wj1-jint.Wj0)
    elif jint.Wj1 <= i[1]:
        qW = 0
    
    return qR * qW

def GP(ngp: int)->List[float]:
    """ input: ngp
        output: GP
    """
    if ngp == 1:
        return [0.]
    elif ngp == 2:
        return [-0.577350269189626, 0.577350269189626]
    elif ngp == 3:
        return [0., -0.774596669241483, 0.774596669241483]
    elif ngp == 4:
        return [-0.339981043584856, 0.339981043584856, -0.861136311594053, 0.861136311594053]
    elif ngp == 5:
        return [0., -0.538469310105683, 0.538469310105683, -0.906179845938664, 0.906179845938664]
    elif ngp == 6:
        return [-0.238619186083197, 0.238619186083197, -0.661209386466265, 0.661209386466265, -0.932469514203152, 0.932469514203152]
    
def GW(ngp: int) -> List[float]:
    """" input: ngp
        output: GW = [GW1, GW2, GW3, GW4, GW5, GW6]
    """
    if ngp == 1:
        return [2.]
    elif ngp == 2:
        return [1., 1.]
    elif ngp == 3:
        return [0.888888888888889, 0.555555555555556, 0.555555555555556]
    elif ngp == 4:
        return [0.652145154862546, 0.652145154862546, 0.347854845137454, 0.347854845137454]
    elif ngp == 5:
        return [0.568888888888889, 0.478628670499366, 0.478628670499366, 0.236926885056189, 0.236926885056189]
    elif ngp == 6:
        return [0.467913934572691, 0.467913934572691, 0.360761573048139, 0.360761573048139, 0.171324492379170, 0.171324492379170]
    
def shp(psietazeta: List[float]) -> List[float]:
    """ input: psi, eta, zeta
        output: N = [N1, N2, N3, N4, N5, N6, N7, N8]
    """
    # TODO: check
    psi = psietazeta[0]
    eta = psietazeta[1]
    zeta = psietazeta[2]
    N = [0.] * 8
    N[0] = 0.125*(1-psi)*(1-eta)*(1-zeta)
    N[1] = 0.125*(1+psi)*(1-eta)*(1-zeta)
    N[2] = 0.125*(1+psi)*(1+eta)*(1-zeta)
    N[3] = 0.125*(1-psi)*(1+eta)*(1-zeta)
    N[4] = 0.125*(1-psi)*(1-eta)*(1+zeta)
    N[5] = 0.125*(1+psi)*(1-eta)*(1+zeta)
    N[6] = 0.125*(1+psi)*(1+eta)*(1+zeta)
    N[7] = 0.125*(1-psi)*(1+eta)*(1+zeta)
    return N

def Dshp(psietazeta: List[float]) -> List[List[float]]:
    """ input: psi, eta, zeta
        output: Dshp = [[dN/dpsi], [dN/deta], [dN/dzeta]]
    """
    psi = psietazeta[0]
    eta = psietazeta[1]
    zeta = psietazeta[2]
    N = []
    N.append([
        -0.125*(1-eta)*(1-zeta),
        0.125*(1-eta)*(1-zeta),
        0.125*(1+eta)*(1-zeta),
        -0.125*(1+eta)*(1-zeta),
        -0.125*(1-eta)*(1+zeta),
        0.125*(1-eta)*(1+zeta),
        0.125*(1+eta)*(1+zeta),
        -0.125*(1+eta)*(1+zeta),
    ])
    N.append([
        -0.125*(1-psi)*(1-zeta),
        -0.125*(1+psi)*(1-zeta),
        0.125*(1+psi)*(1-zeta),
        0.125*(1-psi)*(1-zeta),
        -0.125*(1-psi)*(1+zeta),
        -0.125*(1+psi)*(1+zeta),
        0.125*(1+psi)*(1+zeta),
        0.125*(1-psi)*(1+zeta)
    ])
    N.append([
        -0.125*(1-psi)*(1-eta),
        -0.125*(1+psi)*(1-eta),
        -0.125*(1+psi)*(1+eta),
        -0.125*(1-psi)*(1+eta),
        0.125*(1-psi)*(1-eta),
        0.125*(1+psi)*(1-eta),
        0.125*(1+psi)*(1+eta),
        0.125*(1-psi)*(1+eta)
    ])
    return N

def calc_eps(bb: List[List[float]], disp: List[List[float]]) -> List[List[float]]:
    """" Calculate strain tensor """
    _ = []
    for b in bb:
        _.append(
            [
                np.dot(b[0], np.transpose(disp)[0]),
                np.dot(b[1], np.transpose(disp)[1]),
                np.dot(b[2], np.transpose(disp)[2]),
                np.dot(b[0], np.transpose(disp)[1]) + np.dot(b[1], np.transpose(disp)[0]),
                np.dot(b[0], np.transpose(disp)[2]) + np.dot(b[2], np.transpose(disp)[0]),
                np.dot(b[1], np.transpose(disp)[2]) + np.dot(b[2], np.transpose(disp)[1]),
            ]
        )

    return _

def calc_Du(bb, disp):
    """" Calculate displacement gradient """
    _ = []
    for b in bb:
        _.append(np.dot(b[0], disp))
    return _

def calc_Dq(bb, qi):
    """" Calculate displacement gradient """
    _ = []
    for b1 in bb:
        _.append([np.dot(b2, qi) for b2 in b1])
    return _

def flatten(array):
    return [j for i in array for j in i]

def mapthread(func, array):
    n = len(array[0])
    for i in array:
        if n != len(i):
            raise ValueError("len(a) != len(b)")
    return [func(*array[:][i]) for i in range(n)]

def jintegral(step, l):
    logger.info(f"step: {step} :: CALCULATE J-INTEGRAL")
    posLzx0 = l.node_by_line # nodes of crack line

    posLzx0X = [0] * len(posLzx0) # nodes of crack line (List[List[float]]])
    for i in range(len(posLzx0)):
        revd_posLzx0 = posLzx0[i][1:][::-1]
        sym_posLzx0 = np.concatenate([revd_posLzx0, posLzx0[i]])
        _ = [0] * len(sym_posLzx0)
        for j in range(len(sym_posLzx0)):
            if j < (len(sym_posLzx0) + 1) / 2:
                _[j] = [-sym_posLzx0[j][0], sym_posLzx0[j][1]]
            else:
                _[j] = [sym_posLzx0[j][0], sym_posLzx0[j][1]]
        posLzx0X[i] = _
    
    nposLzxX = [len(posLzx0X[i]) for i in range(len(posLzx0X))] # number of nodes of crack line
    posLzxX = flatten(posLzx0X) # all nodes
    nnpLzxX = len(posLzxX) # number of all nodes
    
    nback = const_local_mesh.nback
    nfront = const_local_mesh.nfront
    enLzxX = []

    for j in range(nback+nfront): # nback + nfront = len(nposLzxX): number of crack lines
        for i in range(1, nposLzxX[j+1]):
            if j == 0:
                enLzxX.append([i, i + 1, nposLzxX[j] + i + 1, nposLzxX[j] + i])
            else:
                s = sum(nposLzxX[:j])
                if nposLzxX[j+1] < nposLzxX[j]:
                    enLzxX.append([
                        s+1+i,
                        s+1+i+1,
                        s+1+nposLzxX[j]+i,
                        s+1+nposLzxX[j]+i-1])
                else:
                    enLzxX.append([
                        s+i,
                        s+i+1,
                        s+nposLzxX[j]+i+1,
                        s+nposLzxX[j]+i])
    ndyL = const_local_mesh.ndyL
    nelLzxX = len(enLzxX)
    enLX = [
        [enLzxX[i][k]+nnpLzxX*(j-1) for k in range(4)]+[enLzxX[i][k]+nnpLzxX*j for k in range(4)] for j in range(1, ndyL+1) for i in range(nelLzxX)
        ] # all elements

    elesizeyL = const_local_mesh.elesizeyL
    newposLX = [[posLzxX[i][1]*0.001, elesizeyL*j*0.001, posLzxX[i][0]*0.001] for j in range(ndyL+1) for i in range(nnpLzxX)] # all nodes
    nnpLX = len(newposLX) # number of all nodes

    aL = const_local_mesh.aL # -1: from number to index
    hLz = const_local_mesh.hLz
    cotip = l.cotip

    numposLzx0X = [0] * len(posLzx0X)
    key = 0
    for i in range(len(posLzx0X)):
        shift = [0] * len(posLzx0X[i])
        for j in range(len(posLzx0X[i])):
            key += 1
            shift[j] = [key] + posLzx0X[i][j]
        numposLzx0X[i] = shift

    if [idx for idx, val in enumerate(numposLzx0X[aL]) if math.isclose(val[2], cotip[1])] == []:
        logger.error(f"""
        Error: cotip is not in numposLzx0X[{aL}].
        cotip: {cotip[1]}
        numposLzx0X[{aL-1}]: {numposLzx0X[aL-1]}
        numposLzx0X[{aL}]: {numposLzx0X[aL]}
        numposLzx0X[{aL+1}]: {numposLzx0X[aL+1]}
        """)
        return
    poscotip = numposLzx0X[aL][[idx for idx, val in enumerate(numposLzx0X[aL]) if math.isclose(val[2], cotip[1])][0]][0]
    front = range(poscotip-int(15/hLz), poscotip+int(15/hLz)+1)

    poscotip1 = numposLzx0X[aL-1][(int((len(numposLzx0X[aL-1])+1)//2))-1][0]
    poscotip2 = numposLzx0X[aL-2][(int((len(numposLzx0X[aL-2])+1)//2))-1][0]
    poscotip01 = numposLzx0X[aL+1][(len(numposLzx0X[aL+1])+1)//2-1][0]
    poscotip02 = numposLzx0X[aL+2][(len(numposLzx0X[aL+2])+1)//2-1][0]

    str_step = step2str(step)
    username = sim_params.USER_NAME
    reponame = sim_params.REPO_NAME
    dirname = sim_params.UBUNTU_DIR
    dirnametest = sim_params.DIR_NAME_TEST
    day = sim_params.DAY
    path = f"/home/{username}/{reponame}/arrest/generation/Newton/{dirnametest}/{day}/step{str_step}"
    # indexを削除する
    disLG = np.delete(np.loadtxt(fname=path+"/log/u_gl.l.dat", skiprows=1), obj=0, axis=1)
    acceLG = np.delete(np.loadtxt(fname=path+"/log/a.l.dat", skiprows=1), obj=0, axis=1)

    ndoflist = np.loadtxt(fname=path+"/ndof_list.txt")

    ndoflistX = [0] * nnpLX
    ndof6 = len([x for x in ndoflist if x == 6])
    nodeLdat = l.node_dat
    _ = [idx for idx in range(1, len(nodeLdat)) if math.isclose(nodeLdat[idx][3], 0.)]
    ndof6z0 = len([x for x in [ndoflist[i] for i in _] if x == 6])
    disLGX = [[0,0,0]] * (nnpLX+ndof6*2-ndof6z0)
    acceLGX = [[0,0,0]] * (nnpLX+ndof6*2-ndof6z0)
    numposLzx0 = [0] * len(posLzx0)
    
    key = 0
    for i in range(len(posLzx0)):
        shift = [0] * len(posLzx0[i])
        for j in range(len(posLzx0[i])):
            key += 1
            shift[j] = [key] + [posLzx0[i][j]]
        numposLzx0[i] = shift
    idcorzx = [0] * len(numposLzx0)
    for i in range(len(numposLzx0)):
        idcorzx[i] = [x[0] for x in numposLzx0[i][1:][::-1]] + [x[0] for x in numposLzx0[i]]
    
    def convert(xyz):
        return [xyz[0], xyz[1], -xyz[2]]

    nnpLzx = l.num_node_zx
    HL = const_local_mesh.HL
    idcor = [num + nnpLzx*i for i in range(HL+1) for _ in idcorzx for num in _]
    
    num = 0
    for i in range(nnpLX):
        if newposLX[i][2] < 0.:
            is_newnode = True
        else:
            is_newnode = False

        id = idcor[i] - 1
        idx = int(sum(ndoflist[:id])//3)
        if ndoflist[id] == 3:
            if is_newnode == False:
                disLGX[num] = disLG[idx]
                acceLGX[num] = acceLG[idx]
            else:
                disLGX[num] = convert(disLG[idx])
                acceLGX[num] = convert(acceLG[idx])

            ndoflistX[i] = 3
            num += 1
        elif ndoflist[id] == 6:
            if is_newnode == False:
                disLGX[num] = disLG[idx]
                disLGX[num+1] = disLG[idx+1]
                acceLGX[num] = acceLG[idx]
                acceLGX[num+1] = acceLG[idx+1]
            else:
                disLGX[num] = convert(disLG[idx])
                disLGX[num+1] = convert(disLG[idx+1])
                acceLGX[num] = convert(acceLG[idx])
                acceLGX[num+1] = convert(acceLG[idx+1])

            ndoflistX[i] = 6
            num += 2
        else:
            ValueError("ndof is not 3 or 6")

    def _try(nL):
        x0 = newposLX[nL-1][0]
        z0 = newposLX[nL-1][2]
        dx = -x0
        dz = -z0
        nl = nL - poscotip

        x1 = newposLX[poscotip2+nl][0]
        z1 = newposLX[poscotip2+nl][2]

        x01 = newposLX[poscotip02+nl][0]
        z01 = newposLX[poscotip02+nl][2]

        if x1 < -0.035 or 0.015 < z01:
            return
        else:
            pass

        idx = [_ for _ in range(len(numposLzx0X[aL])) if numposLzx0X[aL][_][0]==nL][0]
        if len(posLzx0X[aL+1]) < len(posLzx0X[aL]):
            angle = theta(posLzx0X[aL], posLzx0X[aL+1])[idx-1]
        else:
            angle = theta(posLzx0X[aL], posLzx0X[aL+1])[idx]

        xyz = [[
            math.cos(angle)*_[0] + math.sin(angle)*_[2] + math.cos(angle)*dx + math.sin(angle)*dz,
            -math.sin(angle)*_[0] + math.cos(angle)*_[2] - math.sin(angle)*dx + math.cos(angle)*dz,
            _[1]
        ] for _ in newposLX]

        lL = const_local_mesh.lL
        nnpL = l.num_node # num of node in local mesh
        posLzx0xyz = [0] * (nnpL*2-(aL+lL+1)*(HL+1))
        nposlist = [len(posLzx0[i]) for i in range(aL+lL+1)]
        n = 0
        for i in range(1, HL+2):
            for j in range(1, aL+lL+2):
                for k in range(1, nposlist[j-1]*2):
                    posLzx0xyz[n] = [j-aL-1, i-1, k-nposlist[j-1]-(nL-poscotip)]
                    n += 1
        posLzx0rz = [[
            math.sqrt(posLzx0xyz[i][0]**2 + posLzx0xyz[i][1]**2),
            abs(posLzx0xyz[i][2])
        ] for i in range(nnpL*2-(aL+lL+1)*(HL+1))]

        # calculate q
        
        qi0 = list(map(q0, posLzx0rz))
        qi0fornt = list(map(q0, [posLzx0rz[i] for i in front]))

        hL = const_local_mesh.hL
        meas = 0
        for i in range(len(qi0fornt)):
            meas += sum(qi0fornt[i:i+2])*hL*0.0005
        qe0 = list(map(sum, [[qi0[i-1] for i in l] for l in enLX]))
        nq = []

        for idx, val in enumerate(qe0):
            if val > 10e-8: nq.append(idx)
        
        elemq = [enLX[i] for i in nq]
        enode = np.array([[xyz[i-1] for i in l] for l in elemq])
        
        ngp = jint.ngp
        psietazeta = [[k, j, i] for i in GP(ngp) for j in GP(ngp) for k in GP(ngp)]
        weight = flatten([[k, j, i] for i in GW(ngp) for j in GW(ngp) for k in GW(ngp)])
        nn = np.array(list(map(shp, psietazeta)))
        Dnn = np.array(list(map(Dshp, psietazeta)))

        J = []
        detJ = []
        for e in enode:
            j = []
            for d in Dnn:
                j.append(np.dot(d, e))
            J.append(j)
            detJ.append(-np.linalg.det(j))

        EE = jint.ee
        nu = jint.Nu
        de = EE/((1.+nu)*(1.-2.*nu)) * np.array([
            [1.-nu, nu, nu, 0., 0., 0.],
            [nu, 1.-nu, nu, 0., 0., 0.],
            [nu, nu, 1.-nu, 0., 0., 0.],
            [0., 0., 0., 0.5-nu, 0., 0.],
            [0., 0., 0., 0., 0.5-nu, 0.],
            [0., 0., 0., 0., 0., 0.5-nu]
        ])

        bb = []
        for j1 in J:
            _ = []
            for j2, dnn in zip(j1, Dnn):
                _.append(np.dot(np.linalg.inv(j2), dnn))
            bb.append(_)

        def uxyz(disp):
            ux = disp[0]*math.cos(angle) + disp[2]*math.sin(angle)
            uy = -disp[0]*math.sin(angle) + disp[2]*math.cos(angle)
            uz = disp[1]
            return [ux, uy, uz]
        
        dispR = list(map(uxyz, disLGX))
        acceR = list(map(uxyz, acceLGX))
        dispqe = []
        acceqe = []
        dispqerich = []
        accerich = []

        for i in range(len(elemq)):
            d = []
            a = []
            den = []
            aen = []
            for j in range(8):
                idx = int(sum(ndoflistX[:elemq[i][j]-1])//3)
                if ndoflistX[elemq[i][j]-1] == 3:
                    d.append(dispR[idx])
                    a.append(acceR[idx])
                    den.append([0., 0., 0.])
                    aen.append([0., 0., 0.])
                elif ndoflistX[elemq[i][j]-1] == 6:
                    d.append(dispR[idx])
                    a.append(acceR[idx])
                    den.append(dispR[idx+1])
                    aen.append(acceR[idx+1])
                elif ndoflistX[elemq[i][j]-1] == 0:
                    raise ValueError("ndoflistX[elemq[i][j]] is not to be 0")
                else:
                    raise ValueError("ndoflistX is wrong")
            dispqe.append(d)
            acceqe.append(a)
            dispqerich.append(den)
            accerich.append(aen)
        accee = np.array([[np.dot(n, a)] for a in acceqe for n in nn]) + np.array([[np.dot(n, a)] for a in accerich for n in nn])
        epsbf = np.array([calc_eps(b, d) for b, d in zip(bb, dispqe)])
        epsenrich = np.array([calc_eps(b, d) for b, d in zip(bb, dispqerich)])
        eps = epsbf + epsenrich
        deltabf = np.array([[np.dot(de, e2) for e2 in e1] for e1 in epsbf])
        deltaenrich = np.array([[np.dot(de, e2) for e2 in e1] for e1 in epsenrich])
        delta = deltabf + deltaenrich

        Dubf = np.array([calc_Du(b, d) for b, d in zip(bb,dispqe)])
        Duenrich = np.array([calc_Du(b, d) for b, d in zip(bb,dispqerich)])
        Du = Dubf + Duenrich

        qie = [[qi0[i-1] for i in l] for l in elemq]
        Dq = np.array([calc_Dq(b, d) for b, d in zip(bb, qie)])

        def Jint0s(eps, delta, Du, Dq, detJ):
            _ = []
            for e, d, du, dq, dj, w in zip(eps, delta, Du, Dq, detJ, weight):
                _.append(
                    (
                    (d[0]*du[0] + d[3]*du[1] + d[4]*du[2] - 0.5*np.dot(e,d)) * np.array(dq[0]) +
                    (d[3]*du[0] + d[1]*du[1] + d[5]*du[2]) * np.array(dq[1]) +
                    (d[4]*du[0] + d[5]*du[1] + d[2]*du[2]) * np.array(dq[2])
                    ) * dj * w
                )
            return sum(_)
        
        Jintsv = []
        for e, d, du, dq, dj in zip(eps,delta,Du,Dq,detJ):
            Jintsv.append(Jint0s(e,d,du,dq,dj))
        # logger.info(f"Jintsv: {Jintsv}")
        Jints = sum(Jintsv)
        # logger.info(f"meas: {meas}")
        Jints /= meas
        # logger.info(f"Jints = {Jints}")

        rho = jint.Rho
        def Jint0d(acce, Du, detJ):
            _ = []
            for a, du, dj, w in zip(acce, Du, detJ, weight):
                _.append(
                    (
                    (a[0]*du[0] + a[1]*du[1] + a[2]*du[2]) * rho * np.array(dj) * w
                    )
                )
            return sum(_)
        
        # logger.info(f"accee: {np.array(accee).shape}")
        Jintdv = []
        for a, du, dj in zip(accee, Du, detJ):
            Jintdv.append(Jint0d(a,du,dj))

        Jintd = sum(Jintdv)

        Jint = Jints + Jintd

        return Jint
    
    Jlist = []
    # logger.info(f"range: {range(poscotip, poscotip+int(15./hLz)-1)}")
    for a in range(poscotip, poscotip+int(15./hLz)-1):
        Jval = _try(a)
        Jlist.append(Jval)
        # break
    
    dirnametest = sim_params.DIR_NAME_TEST
    day = sim_params.DAY
    path = f"/home/{username}/{reponame}/arrest/generation/Newton/{dirnametest}/{day}/step{str_step}/Jlist.dat"
    with open (path, mode="w+") as f:
        for line in Jlist:
            f.write(str(line)+"\n")
    f.close()

    path = f"/home/{username}/{reponame}/arrest/generation/Newton/{dirnametest}/{day}/Jlistall.dat"
    with open (path, mode="a") as f:
        for line in Jlist:
            if line is None:
                pass
            else:
                f.write(str(line)+" ")
        f.write("\n")
    f.close()

    return Jlist
