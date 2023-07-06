import math
from scipy import optimize
import sys
import pickle
from const import simulation_params as sim_params
from const import const_local_mesh
import local_mesh
from utils.judge import judge
from utils.step2str import step2str
from utils.logger import logger

def initial(step: int, l, g) -> None:
    REstart = sim_params.REstart
    if step == 0:
        return
    posL = l.node
    nnpL = len(posL)
    enL = l.elem
    posLzx = l.node_zx
    edgenode = g.edgenode
    nnpLzx = len(posLzx)
    positiveorangebf = [-1] * len(posL)
    negativeorangebf = [-1] * len(posL)

    for i in range(len(enL)):
        ii, jj, kk, mm = enL[i][0]-1, enL[i][1]-1, enL[i][5]-1, enL[i][4]-1
        nn, oo, pp, qq = enL[i][3]-1, enL[i][2]-1, enL[i][6]-1, enL[i][7]-1

        xi, yi, zi = posL[ii][0], posL[ii][1], posL[ii][2]
        xj, yj, zj = posL[jj][0], posL[jj][1], posL[jj][2]
        xk, yk, zk = posL[kk][0], posL[kk][1], posL[kk][2]
        xm, ym, zm = posL[mm][0], posL[mm][1], posL[mm][2]
        xn, yn, zn = posL[nn][0], posL[nn][1], posL[nn][2]
        xo, yo, zo = posL[oo][0], posL[oo][1], posL[oo][2]
        xp, yp, zp = posL[pp][0], posL[pp][1], posL[pp][2]
        xq, yq, zq = posL[qq][0], posL[qq][1], posL[qq][2]

        judge_i = judge(xi, yi, zi, edgenode)
        judge_j = judge(xj, yj, zj, edgenode)
        judge_k = judge(xk, yk, zk, edgenode)
        judge_m = judge(xm, ym, zm, edgenode)
        judge_n = judge(xn, yn, zn, edgenode)
        judge_o = judge(xo, yo, zo, edgenode)
        judge_p = judge(xp, yp, zp, edgenode)
        judge_q = judge(xq, yq, zq, edgenode)

        if (not (judge_i and judge_j and judge_k and judge_m and judge_n and judge_o and judge_p and judge_q)) and (judge_i or judge_j or judge_k or judge_m or judge_n or judge_o or judge_p or judge_q):
            if judge_i:
                negativeorangebf[ii] = ii
            else:
                positiveorangebf[ii] = ii
            if judge_j:
                negativeorangebf[jj] = jj
            else:
                positiveorangebf[jj] = jj
            if judge_k:
                negativeorangebf[kk] = kk
            else:
                positiveorangebf[kk] = kk
            if judge_m:
                negativeorangebf[mm] = mm
            else:
                positiveorangebf[mm] = mm
            if judge_n:
                negativeorangebf[nn] = nn
            else:
                positiveorangebf[nn] = nn
            if judge_o:
                negativeorangebf[oo] = oo
            else:
                positiveorangebf[oo] = oo
            if judge_p:
                negativeorangebf[pp] = pp
            else:
                positiveorangebf[pp] = pp
            if judge_q:
                negativeorangebf[qq] = qq
            else:
                positiveorangebf[qq] = qq

    negativeorange = [x for x in negativeorangebf if x >= 0]
    positiveorange = [x for x in positiveorangebf if x >= 0]
    negative = len(negativeorange)
    positive = len(positiveorange)

    Ora = positiveorangebf + negativeorangebf
    orange = sorted(list(set(positiveorange).union(set(negativeorange))))
    orangeDIM = len(orange)
    ndoflistnowL = [3] * nnpL
    for i in range(orangeDIM):
        ndoflistnowL[orange[i]] = 6

    user_name = sim_params.USER_NAME
    repo_name = sim_params.REPO_NAME
    dirnametest = sim_params.DIR_NAME_TEST
    interm = sim_params.INTERM
    str_step_pre = step2str(step-interm-1)
    path = f"/home/{user_name}/arrest/generation/inputfiles/step{str_step_pre}/"
    with open(path+"local_mesh.pickel", "rb") as f:
        l_bf = pickle.load(f)
    posLzxps = l_bf.node_zx
    nnpLzxps = len(posLzxps)
    enLzxps = l_bf.enLzx
    nelLzxps = len(enLzxps)

    day = sim_params.DAY
    str_step = step2str(step-1)
    path = f"/home/{user_name}/arrest/generation/Newton/{dirnametest}/{day}/step{str_step}/"
    ndoflistpsL = []
    with open(path+"ndof_list.txt", "r") as f:
        lines = f.readlines()
        for line in lines:
            ndoflistpsL.append(int(line))

    ndoflistpsLpart  = []
    cnt = 0
    while cnt+nnpLzxps <= len(ndoflistpsL):
        ndoflistpsLpart.append(ndoflistpsL[cnt:cnt+nnpLzxps])
        cnt += nnpLzxps

    positiveorangeps = []
    with open(path+"positiveorange.txt", "r") as f:
        lines = f.readlines()
        for line in lines:
            positiveorangeps.append(int(line))

    def load(filename):
        ans = []
        path = f"/home/{user_name}/arrest/generation/Newton/{dirnametest}/{day}/step{str_step}/"
        with open(path+filename, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                ans.append(list(map(float, line.split()))[1:])
        return ans

    try:
        disG = load("log/u.g.dat")
    except:
        logger.info(f"load log/u.g_{str_step}.dat")
        disG = load(f"log/u.g_{str_step}.dat")
    try:
        disL = load("log/u.l.dat")
    except:
        logger.info(f"load log/u.l_{str_step}.dat")
        disL = load(f"log/u.l_{str_step}.dat")
    try:
        velG = load("log/v.g.dat")
    except:
        logger.info(f"load log/v.g_{str_step}.dat")
        velG = load(f"log/v.g_{str_step}.dat")
    try:
        velL = load("log/v.l.dat")
    except:
        logger.info(f"load log/v.l_{str_step}.dat")
        velL = load(f"log/v.l_{str_step}.dat")
    try:
        acceG = load("log/a.g.dat")
    except:
        logger.info(f"load log/a.g_{str_step}.dat")
        acceG = load(f"log/a.g_{str_step}.dat")
    try:
        acceL = load("log/a.l.dat")
    except:
        logger.info(f"load log/a.l_{str_step}.dat")
        acceL = load(f"log/a.l_{str_step}.dat")

    disiniG = disG
    veliniG = velG
    acceiniG = acceG

    zxLps = [[posLzxps[j-1] for j in i] for i in enLzxps]
    zmaxL = [max([z[0] for z in zxLps[i]]) for i in range(nelLzxps)]
    zminL = [min([z[0] for z in zxLps[i]]) for i in range(nelLzxps)]
    xmaxL = [max([x[1] for x in zxLps[i]]) for i in range(nelLzxps)]
    xminL = [min([x[1] for x in zxLps[i]]) for i in range(nelLzxps)]

    def N1(xi, eta): return (1 - xi) * (1 - eta) * 0.25
    def N2(xi, eta): return (1 + xi) * (1 - eta) * 0.25
    def N3(xi, eta): return (1 + xi) * (1 + eta) * 0.25
    def N4(xi, eta): return (1 - xi) * (1 + eta) * 0.25

    def intpo(pLzx):
        elLpsabb = [i for i in range(nelLzxps)
                    if zminL[i] <= posLzx[pLzx][0] <= zmaxL[i] and xminL[i] <= posLzx[pLzx][1] <= xmaxL[i]]
        nelLpsabb = len(elLpsabb)
        if nelLpsabb == 0:
            return 0
        else:
            def solvexietaL(eL, XY):
                XY = [i for i in XY]
                xy = [i for i in zxLps[eL]]
                def func(xieta):
                    xi = xieta[0]
                    eta = xieta[1]
                    return [
                        N1(xi, eta) * xy[0][0] + N2(xi, eta) * xy[1][0] + N3(xi, eta) * xy[2][0] + N4(xi, eta) * xy[3][0] - XY[0],
                        N1(xi, eta) * xy[0][1] + N2(xi, eta) * xy[1][1] + N3(xi, eta) * xy[2][1] + N4(xi, eta) * xy[3][1] - XY[1]
                    ]
                return optimize.root(func, x0=[0., 0.], method="hybr").x
            xietaLcand = [[solvexietaL(i, posLzx[pLzx])] for i in elLpsabb]
            nxietaLcand = [len(i) for i in xietaLcand]
            xietaLabb = []
            for j in range(nelLpsabb):
                for i in range(nxietaLcand[j]):
                    a = round(xietaLcand[j][i][0], 5)
                    b = round(xietaLcand[j][i][1], 5)
                    if -1 <= a <= 1 and -1 <= b <= 1:
                        xietaLabb.append([elLpsabb[j], xietaLcand[j][i]])
            if xietaLabb == []:
                return 0
            else:
                xietaL0 = xietaLabb[0]
                return [xietaL0[0], xietaL0[1]]

    def internewenrich(array, eL, xieta, layer, type):
        if type == "not":
            f, fx = 1, 1
        elif type == "normal":
            f, fx = 1, 0
        elif type == "X":
            f, fx = 0, 1
        else:
            logger.error("type is not defined")
            sys.exit()

        n = [
            N1(xieta[0], xieta[1]),
            N2(xieta[0], xieta[1]),
            N3(xieta[0], xieta[1]),
            N4(xieta[0], xieta[1])
        ]

        ans = [0, 0, 0]
        for i in range(4):
            shift = len([x for x in ndoflistpsLpart[layer][:enLzxps[eL][i]-1] if x == 6])
            arrayLi = array[enLzxps[eL][i]+shift-1]
            if ndoflistpsLpart[layer][enLzxps[eL][i]-1] == 3:
                arrayLx = [0, 0, 0]
            else:
                if len([x for x in positiveorangeps if x == enLzxps[eL][i]]) == 1 or type == "normal" or type == "X":
                    arrayLx = [a for a in array[enLzxps[eL][i]+shift]]
                else:
                    arrayLx = [-a for a in array[enLzxps[eL][i]+shift]]
            ans = [ans[j]+n[i]*(fx*arrayLx[j]+f*arrayLi[j]) for j in range(3)]
        return ans
    
    xietaLzx = [intpo(i) for i in range(nnpLzx)]
    modli = []
    for i in range(nnpL):
        modli.append(i%nnpLzx)
    ndyL = const_local_mesh.ndyL
    count = [0] * (ndyL+1)
    for j in range(1, ndyL+2):
        for i in range(1, nnpLzxps+1):
            if nnpLzxps*(j-1)+i-1 >= len(ndoflistpsL):
                break
            if ndoflistpsL[nnpLzxps*(j-1)+i-1] == 6:
                count[j-1] += 1
    count = [c+nnpLzxps for c in count]
    disLpart = []
    velLpart = []
    acceLpart = []
    cnt = 0
    for i in range(len(count)):
        disLpart.append(disL[cnt:cnt+count[i]])
        velLpart.append(velL[cnt:cnt+count[i]])
        acceLpart.append(acceL[cnt:cnt+count[i]])
        cnt += count[i]

    disiniL = [0]*(orangeDIM+nnpL)
    veliniL = [0]*(orangeDIM+nnpL)
    acceiniL = [0]*(orangeDIM+nnpL)
    temp = 0
    for i in range(1, nnpL+1):
        idx = math.ceil(i/nnpLzx)-1
        if ndoflistnowL[i-1] == 3:
            if xietaLzx[modli[i-1]] == 0:
                disiniL[i+temp-1] = [0, 0, 0]
                veliniL[i+temp-1] = [0, 0, 0]
                acceiniL[i+temp-1] = [0, 0, 0]
            else:
                disiniL[i+temp-1] = internewenrich(disLpart[idx], xietaLzx[modli[i-1]][0], xietaLzx[modli[i-1]][1], idx, "not")
                veliniL[i+temp-1] = internewenrich(velLpart[idx], xietaLzx[modli[i-1]][0], xietaLzx[modli[i-1]][1], idx, "not")
                acceiniL[i+temp-1] = internewenrich(acceLpart[idx], xietaLzx[modli[i-1]][0], xietaLzx[modli[i-1]][1], idx, "not")
        else:
            if xietaLzx[modli[i-1]] == 0:
                disiniL[i+temp-1] = [0, 0, 0]
                disiniL[i+temp] = [0, 0, 0]
                veliniL[i+temp-1] = [0, 0, 0]
                veliniL[i+temp] = [0, 0, 0]
                acceiniL[i+temp-1] = [0, 0, 0]
                acceiniL[i+temp] = [0, 0, 0]
            else:
                disiniL[i+temp-1] = internewenrich(disLpart[idx], xietaLzx[modli[i-1]][0], xietaLzx[modli[i-1]][1], idx, "normal")
                disiniL[i+temp] = internewenrich(disLpart[idx], xietaLzx[modli[i-1]][0], xietaLzx[modli[i-1]][1], idx, "X")
                veliniL[i+temp-1] = internewenrich(velLpart[idx], xietaLzx[modli[i-1]][0], xietaLzx[modli[i-1]][1], idx, "normal")
                veliniL[i+temp] = internewenrich(velLpart[idx], xietaLzx[modli[i-1]][0], xietaLzx[modli[i-1]][1], idx, "X")
                acceiniL[i+temp-1] = internewenrich(acceLpart[idx], xietaLzx[modli[i-1]][0], xietaLzx[modli[i-1]][1], idx, "normal")
                acceiniL[i+temp] = internewenrich(acceLpart[idx], xietaLzx[modli[i-1]][0], xietaLzx[modli[i-1]][1], idx, "X")
            temp += 1
    
    
    str_step = step2str(step)
    path = f"/home/{user_name}/arrest/generation/inputfiles/step{str_step}/"
    def write(filename, array):
        with open(path+filename, "w") as f:
            f.write(f"{len(array)}\n")
            for i in range(len(array)):
                line = [str(b) if b != 0. or b != 0 else "0." for b in array[i]]
                f.write(f"{int(i)+1}\t{line[0]}\t{line[1]}\t{line[2]}\n")
        f.close()

    write("init.u.g.dat", disiniG)
    write("init.u.l.dat", disiniL)
    write("init.v.g.dat", veliniG)
    write("init.v.l.dat", veliniL)
    write("init.a.g.dat", acceiniG)
    write("init.a.l.dat", acceiniL)
    
    return
