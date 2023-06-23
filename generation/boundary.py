import numpy as np
import math
import global_mesh
import local_mesh
import experiments_data
from const import simulation_params as sim_params
from const import const_global_mesh
from const import const_local_mesh
from utils.judge import judge
from utils.logger import logger

class Boundary:
    def __init__(self, l) -> None:
        self.cotip = l.cotip
        self.g = global_mesh.GlobalMesh()

    def define_global_boundary(self):
        g = self.g
        node = g.node
        nGz1 = const_global_mesh.nGz1
        nnmGxy = g.num_node_xy
        nodeGxy = g.node_xy
        isgetctod = sim_params.isgetctod
        ctod = experiments_data.ctod

        xfixG = []
        for idx in range(len(node)):
            if np.abs(np.array(node[idx][0]) - 0.115) < 1e-4: xfixG.append(idx+1)
        nxfixG = len(xfixG)

        npG3Dxz = []
        z1Mesh = const_global_mesh.z1Mesh
        z1 = 15./z1Mesh
        for i in range(nGz1+1):
            _ = []
            for idx in range(len(node)):
                if np.linalg.norm(np.array(node[idx][1:]) - np.array([0, i*z1 * 0.001])) < 0.0001:
                    _.append(idx+1)
            npG3Dxz.append(_)
        nnpG3Dxz = len(npG3Dxz)

        def fixselect(nol):
            node = g.node
            num_node = len(node)
            ucfp = set()
            for num in nol:
                if num < num_node*(nnpG3Dxz-1)/nnpG3Dxz:
                    hoge = num + num_node//nnpG3Dxz
                else: hoge = num
                if local_mesh.fa(self.cotip, node[hoge-1][2]*1000)*0.001 <= node[num-1][0]:
                    ucfp.add(num)
                else: pass
            dcfp = set(nol).difference(ucfp)
            dcfp = sorted(list(dcfp))
            if len(dcfp) == 0:
                return sorted(list(ucfp))
            else:
                maxG = max([node[i-1][0] for i in dcfp])
                idx = [dcfp.index(i) for i in dcfp if node[i-1][0] == maxG]
                maxpos = dcfp[idx[0]]
                return np.hstack(([maxpos], sorted(list(ucfp))))

        yfixG = np.concatenate([fixselect(npG3Dxz[i]) for i in range(nnpG3Dxz)])
        nyfixG = len(yfixG)

        pinposxy = [-57., 8.]
        cgposxy = [-42, 8.]

        pinxy = np.array([i for i in range(nnmGxy) if np.linalg.norm(nodeGxy[i] - pinposxy) <= 5.0])
        cgxy = np.array([i for i in range(nnmGxy) if np.linalg.norm(nodeGxy[i] - cgposxy) <= 1.25])

        pinG = np.concatenate([pinxy + nnmGxy * i + 1 for i in range(nGz1 + 1)])
        npinG = len(pinG)
        cghc = 1
        while cgxy.size == 0:
            cghc += 1
            cgxy = np.array([i for i in range(nnmGxy) if np.linalg.norm(nodeGxy[i] - cgposxy) <= 0. + 1.25 * cghc])

        zfixG = []
        for idx in range(len(node)):
            if np.abs(node[idx][2]) < 0.00001:
                zfixG.append(idx+1)
        nzfixG = len(zfixG)

        self.nfixG = nxfixG + nyfixG + npinG + nzfixG
        self.bcG = np.concatenate([
                            np.column_stack([xfixG, np.ones(nxfixG), np.zeros(nxfixG)]), 
                            np.column_stack([yfixG, 2*np.ones(nyfixG), np.zeros(nyfixG)]), 
                            np.column_stack([pinG, 2*np.ones(npinG), np.where(isgetctod == 1, 1., ctod)*np.ones(npinG)]),
                            np.column_stack([zfixG, 3*np.ones(nzfixG), np.zeros(nzfixG)])])

    def define_local_boundary(self, l):
        num_node_by_line = l.num_node_by_line
        num_node_zx = l.num_node_zx

        posL = l.node
        nnpL = len(posL)

        ndyL = const_local_mesh.ndyL
        nback = const_local_mesh.nback
        xminfixL = [j + num_node_zx * i for i in range(ndyL + 1) for j in range(1, num_node_by_line[0] + 1)]
        xmaxfixL = [j + num_node_zx * i for i in range(ndyL + 1) for j in range(num_node_zx - num_node_by_line[-1] + 1, num_node_zx + 1)]
        yminfixL = list(range(sum(num_node_by_line[:nback]) + 1, num_node_zx + 1))
        ymaxfixL = list(range(num_node_zx * ndyL + 1, nnpL + 1))
        zminfixL = [i+1 for i, node in enumerate(posL) if math.isclose(node[2], 0.)]
        TianyuxyzfixL = set()
        edgenode = self.g.edgenode
        for i in range(len(posL)):
            if judge(posL[i][0], posL[i][1], posL[i][2], edgenode):
                TianyuxyzfixL.add(i+1)

        node_zx = l.node_zx
        node_by_line = l.node_by_line

        Tianyufixinnerbf = set()
        for i in range(len(node_by_line)):
            for j in range(len(node_by_line[i])):
                if node_by_line[i][j][0] > 15.:
                    for k in range(len(node_zx)):
                        if np.linalg.norm(node_zx[k] - node_by_line[i][j]) < 10e-8:
                            Tianyufixinnerbf.add(k+1)
                    break
        Tianyufixinner = []
        HL = const_local_mesh.HL
        for i in range(1, HL + 2):
            Tianyufixinner += [x + (i-1)*num_node_zx for x in Tianyufixinnerbf]

        diff = set(TianyuxyzfixL).difference(set(Tianyufixinner))
        xfixL = set(xminfixL).union(set(xmaxfixL)).union(set(ymaxfixL))
        xfixL = xfixL.union(diff)
        nxfixL = len(xfixL)
        yfixL = set(xminfixL).union(set(xmaxfixL)).union(set(yminfixL)).union(set(ymaxfixL))
        yfixL = yfixL.union(diff)
        nyfixL = len(yfixL)
        zfixL = set(xminfixL).union(set(xmaxfixL)).union(set(ymaxfixL)).union(set(zminfixL))
        zfixL = zfixL.union(diff)
        nzfixL = len(zfixL)
        self.xfixL = sorted(list(xfixL))
        self.yfixL = sorted(list(yfixL))
        self.zfixL = sorted(list(zfixL))
        self.nfixL = nxfixL + nyfixL + nzfixL

    def generate(self):
        with open("bc.g.dat", "w") as f:
            f.write(f"{self.nfixG}\n")
            for line in self.bcG:
                f.write(f"{int(line[0])}\t{int(line[1])}\t{line[2]}\n")
        f.close()

        boundaryLdat =([[x, 1, 0.] for x in self.xfixL] 
                    + [[y, 2, 0.] for y in self.yfixL] 
                    + [[z, 3, 0.] for z in self.zfixL])
        with open("bc.l.dat", "w") as f:
            f.write(f"{self.nfixL}\n")
            for line in boundaryLdat:
                f.write(f"{int(line[0])}\t{int(line[1])}\t{line[2]}\n")
        f.close()