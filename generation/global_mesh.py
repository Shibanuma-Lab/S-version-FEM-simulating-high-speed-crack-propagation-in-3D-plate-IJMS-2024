# reffer to docs\PMMA_variables.mdpf
import numpy as np
import math

from const import const_global_mesh
from const import simulation_params as sim_params
import experiments_data
from utils.bias import bias
from utils.get_gep import get_gep
from utils.step2str import step2str

class GlobalMesh:
    def extendsq(self, p):
        nor = const_global_mesh.nor
        r1len = const_global_mesh.r1len
        r1_mesh = const_global_mesh.r1Mesh
        r1 = r1len / r1_mesh
        s, theta= p
        if theta < np.pi / 4:
            x = ((s - nor) / ((nor + r1) - nor)) * ((nor + r1) * np.sqrt(1 + np.tan(theta)**2)) + (((nor + r1) - s) / ((nor + r1) - nor)) * nor
        else:
            x = ((s - nor) / ((nor + r1) - nor)) * ((nor + r1) * np.sqrt(1 + np.tan(0.5 * np.pi - theta)**2)) + (((nor + r1) - s) / ((nor + r1) - nor)) * nor
        return [x, theta]

    def get_car(self, p):
        nor = const_global_mesh.nor
        r, theta = p
        x = r * np.cos(theta) - 35 - nor
        y = r * np.sin(theta)
        return [x, y]
    
    def __init__(self):
        is_cross = sim_params.is_cross
        x1r = const_global_mesh.x1r
        r1len = const_global_mesh.r1len
        r1_mesh = const_global_mesh.r1Mesh
        x1Mesh = const_global_mesh.x1Mesh
        x2Mesh = const_global_mesh.x2Mesh
        x2r = const_global_mesh.x2r
        x4Mesh = const_global_mesh.x4Mesh
        y1Mesh = const_global_mesh.y1Mesh
        y1r = const_global_mesh.y1r
        z1Mesh = const_global_mesh.z1Mesh
        x9Mesh = const_global_mesh.x9Mesh
        x10Mesh = const_global_mesh.x10Mesh
        x11Mesh = const_global_mesh.x11Mesh
        hole = experiments_data.hole
        radius = experiments_data.radius

        if is_cross == False and (hole == 0 or hole == 0.):
            area_div = 8
            x1min = bias(x1Mesh, x1r, 30. - 7.0)[0][0]
            x1max = x1r * x1min
            x2min = bias(x2Mesh, x2r, 7.0 - radius)[0][0]
            x2max = x2r * x2min
            r1 = r1len / r1_mesh
            x4 = (150. - r1len) / x4Mesh
            y1min = bias(y1Mesh, y1r, 55. - radius - r1len)[1]
            y1max = y1r * y1min
            z1 = 15. / z1Mesh

        elif is_cross == True and (hole == 0 or hole == 0.):
            area_div = 12
            areadiv3D = 2
            x1min = bias(x1Mesh, x1r, 30. - 7.0)
            x1max = x1r * x1min
            x2min = bias(x2Mesh, x2r, 7.0 - radius)
            x2max = x2r * x2min
            r1 = r1len / r1_mesh
            x9 = (29.5 - r1len) / x9Mesh
            x10 = 5.0/x10Mesh
            x11 = 62.5/x11Mesh
            y1min = bias(y1Mesh, y1r, 55. - radius - r1len)
            y1max = y1r * y1min
            z1 = 15./z1Mesh
            z2min = bias(z2_mesh, z2r, 18.5)
            z2max = z2r * z2min

        nGx1 = const_global_mesh.nGx1
        bx1 = const_global_mesh.bx1
        nor = const_global_mesh.nor
        nGx2 = const_global_mesh.nGx2
        bx2 = const_global_mesh.bx2
        nGtheta1 = const_global_mesh.nGtheta1
        nGr1 = const_global_mesh.nGr1
        nGx3 = const_global_mesh.nGx3
        nGy1 = const_global_mesh.nGy1
        by1 = const_global_mesh.by1
        nGz1 = const_global_mesh.nGz1

        bias_x1 = bias(nGx1, bx1, 8.)  # area 1
        hx1_min = bias_x1[0][0]
        bias_per_element_x1 = bias_x1[1]
        hx1_max = bx1 * hx1_min

        bias_x2 = bias(nGx2, bx2, 22. - nor)  # area 2
        hx2_min = bias_x2[0][0]
        bias_per_element_x2 = bias_x2[1]
        hx2max = bx2 * hx2_min

        h_theta1a = ((2*nor**2) * (1 - math.cos((math.pi/2)/nGtheta1))) ** 0.5
        h_theta1b = (nor + r1) / (nGtheta1 / 2)

        hr1 = r1 / nGr1
        hx3 = (150. - r1) / nGx3  # r1=5, the annotated part was originally 97, subtracting notch length
        biasy1 = bias(nGy1, by1, 55. - r1 - nor)  # half of the width, originally 40
        hy1min = biasy1[0][0]
        bias_per_element_y1 = biasy1[1]
        hy1max = by1 * hy1min
        hz1 = 15. / nGz1  # thickness, originally 6.5

        # area1
        node_x1 = get_gep(-65., hx1_max, 1./bias_per_element_x1, nGx1)      # hx1max, bpex1 and nGx1 are already defined earlier
        node_y1 = get_gep(nor, hr1, 1., nGr1)                               # nor, hr1 and nGr1 are already defined earlier
        node_xy1 = [[x, y] for y in node_y1 for x in node_x1]
        num_node_xy1 = len(node_xy1)

        _ien_xy1 = np.array([(j-1) * (nGx1 + 1) + i for i in range(1, nGx1 + 1) for j in range(1, nGr1 + 1)])
        elem_xy1 = [[i, i+1, i + nGx1 + 2, i + nGx1+1] for i in _ien_xy1]
        num_elem_xy1 = len(elem_xy1)

        if sim_params.is_cross:
            node_all_xy = [node_xy1]
            elem_all_xy = [elem_xy1]
            num_node_all_xy = [num_node_xy1]
            num_elem_all_xy = [num_elem_xy1]

            node_xy = np.array(node_all_xy[0])
            elem_xy = np.array(elem_all_xy[0])
            num_node_xy = np.array(num_node_all_xy[0])

        else:
            # area2
            node_x2 = get_gep(-57., hx2max, 1./bias_per_element_x2, nGx2) # hx2max, bpex2 and nGx2 are already defined earlier
            node_y2 = node_y1
            node_xy2 = [[x, y] for y in node_y2 for x in node_x2]
            num_node_xy2 = len(node_xy2)

            _ien_xy2 = np.array([(j-1)*(nGx2+1)+i for i in range(1, nGx2+1) for j in range(1, nGr1+1)])
            elem_xy2 = [[i, i+1, i+nGx2+2, i+nGx2+1] for i in _ien_xy2]
            num_elem_xy2 = len(elem_xy2)

            # area3
            node_r3 = np.arange(nor, nor + r1 + r1/nGr1, r1/nGr1)  # r coordinates
            node_theta3 = np.arange(0.0, 0.5 * np.pi + 0.5 * np.pi / nGtheta1, 0.5 * np.pi / nGtheta1)  # theta coordinates
            node_rtheta3 = np.array([[r, theta] for theta in node_theta3 for r in node_r3])  # r-theta coordinates
            
            node_rthetasq3 = np.array(list(map(self.extendsq, node_rtheta3)))  # r-theta coordinates for square ok

            car_coords = np.array(list(map(self.get_car, node_rtheta3)))  # rectangular coordinates
            node_xy3 = [self.get_car(node) for node in node_rthetasq3]
            num_node_xy3 = len(node_xy3)
            e1g = np.setdiff1d(np.arange(1, (nGr1 + 1) * nGtheta1 + 1), (nGr1 + 1) * np.arange(1, nGtheta1 + 2))
            e2g = np.setdiff1d(np.arange(1, (nGr1 + 1) * nGtheta1 + 1), (nGr1 + 1) * np.arange(0, nGtheta1 + 2) + 1)
            e3g = np.setdiff1d(np.arange(nGr1 + 2, (nGr1 + 1) * (nGtheta1 + 1) + 1), (nGr1 + 1) * np.arange(1, nGtheta1 + 2))
            e4g = np.setdiff1d(np.arange(nGr1 + 2, (nGr1 + 1) * (nGtheta1 + 1) + 1), (nGr1 + 1) * np.arange(0, nGtheta1 + 2) + 1)

            num_elem_xy3 = nGr1 * nGtheta1 #num. of elements
            elem_xy3 = np.vstack([e1g, e2g, e4g, e3g]).T #element-node relationship

            temp = np.zeros(len(node_xy3))
            for i in range(len(node_xy3)):
                if np.linalg.norm(node_xy3[i] - np.array([-35. - nor, 0.])) - nor < 10e-8:
                    temp[i] = i + 1
            _node_xy3 = np.array(node_xy3)
            _edgenode = _node_xy3[temp.nonzero()[0]]
            edgenode = np.vstack(([-65., -nor], [-35., -nor], _edgenode, [-65., nor]))

            # area4
            node_x4 = get_gep(-35. + r1, hx3, 1., nGx3)
            node_y4 = sorted(list(set(_node_xy3[:,1][np.where(np.isclose(_node_xy3[:,0],-35. + r1))])))
            node_xy4 = [[x, y] for y in node_y4 for x in node_x4]
            num_node_xy4 = len(node_xy4)

            _ien_xy4 = np.array([j + (nGx3 + 1) * (i - 1) for j in range(1, nGx3+1) for i in range(1, nGtheta1 // 2 + 1)])
            elem_xy4 = [[i, i+1, i+nGx3+2, i+nGx3+1] for i in _ien_xy4]
            num_elem_xy4 = len(elem_xy4)

            # area5
            node_x5 = node_x1
            node_y5 = get_gep(nor + r1, hy1min, bias_per_element_y1, nGy1)
            node_xy5 = [[x, y] for y in node_y5 for x in node_x5]
            num_node_xy5 = len(node_xy5)

            _ien_xy5 = np.array([(j + (nGx1 + 1) * (i - 1)) for i in range(1, nGy1 + 1) for j in range(1, nGx1 + 1)])
            elem_xy5 = [[i, i+1, i+nGx1+2, i+nGx1+1] for i in _ien_xy5]
            num_elem_xy5 = len(elem_xy5)

            # area6
            node_x6 = node_x2
            node_y6 = node_y5
            node_xy6 = [[x, y] for y in node_y6 for x in node_x6]
            num_node_xy6 = len(node_xy6)
            _ien_xy6 = np.array([(j + (nGx2 + 1) * (i - 1)) for i in range(1, nGy1 + 1) for j in range(1, nGx2 + 1)])
            elem_xy6 = [[i, i+1, i+nGx2+2, i+nGx2+1] for i in _ien_xy6]
            num_elem_xy6 = len(elem_xy6)

            # area7
            node_x7 = sorted(list(set([i for i, j in node_xy3 if abs(j - (nor + r1)) < 0.0001])))
            node_y7 = node_y5
            node_xy7 = [[x, y] for y in node_y7 for x in node_x7]
            num_node_xy7 = len(node_xy7)
            _ien_xy7 = np.array([(j + (nGtheta1 // 2 + 1) * (i - 1)) for i in range(1, nGy1 + 1) for j in range(1, nGtheta1 // 2 + 1)])
            elem_xy7 = [[i, i+1, i+nGtheta1//2+2, i+nGtheta1//2+1] for i in _ien_xy7]
            num_elem_xy7 = len(elem_xy7)

            # area8
            node_x8 = node_x4
            node_y8 = node_y5
            node_xy8 = [[x, y] for y in node_y8 for x in node_x8]
            num_node_xy8 = len(node_xy8)
            _ien_xy8 = np.array([[j + (nGx3 + 1)*(i - 1) for j in range(1, nGx3+1)] for i in range(1, nGy1+1)]).flatten()
            elem_xy8 = [[i, i+1, i+nGx3+2, i+nGx3+1] for i in _ien_xy8]
            num_elem_xy8 = len(elem_xy8)

            node_all_xy = [node_xy1, node_xy2, node_xy3, node_xy4, node_xy5, node_xy6, node_xy7, node_xy8]
            _elem_all_xy = [elem_xy1, elem_xy2, elem_xy3, elem_xy4, elem_xy5, elem_xy6, elem_xy7, elem_xy8]
            elem_all_xy = [np.sort(e, axis=0) for e in _elem_all_xy]
            num_node_all_xy = [num_node_xy1, num_node_xy2, num_node_xy3, num_node_xy4, num_node_xy5, num_node_xy6, num_node_xy7, num_node_xy8]
            num_elem_all_xy = [num_elem_xy1, num_elem_xy2, num_elem_xy3, num_elem_xy4, num_elem_xy5, num_elem_xy6, num_elem_xy7, num_elem_xy8]

            node_xy = node_all_xy[0]
            elem_xy = np.sort(elem_all_xy[0], axis=0)
            num_node_xy = num_node_all_xy[0]

            for area in range(1, area_div):
                com_coord_area = []
                idx_com_coord_present = []
                com_coord_next_area = []
                idx_com_coord_additional = []
                for idx_new_xy in range(len(node_all_xy[area])):
                    for idx_xy in range(len(node_xy)):
                        if np.linalg.norm(np.array(node_xy[idx_xy]) - np.array(node_all_xy[area][idx_new_xy])) <= 10**(-8):
                            xy = list(node_xy[idx_xy])
                            new_xy = list(node_all_xy[area][idx_new_xy])
                            if xy not in com_coord_area:
                                com_coord_area.append(xy)
                                idx_com_coord_present.append(idx_xy)
                            if new_xy not in com_coord_next_area:
                                com_coord_next_area.append(new_xy)
                                idx_com_coord_additional.append(idx_new_xy)
                idx_additional_node = [] # additional node number
                num_node_xy_present = num_node_xy # number of nodes in present nodes
                for node_number in range(1, num_node_all_xy[area]+1):    
                    idx_com_node = np.where(np.array(idx_com_coord_additional) == node_number-1)[0] # common node number in additional nodes
                    if len(idx_com_node) == 0:
                        idx_additional_node.append(node_number + num_node_xy_present)
                    else:
                        nn2 = idx_com_coord_present[idx_com_node[0]] + 1 # idx -> node num
                        idx_additional_node.append(nn2)
                        num_node_xy_present -= 1
                node_xy = np.concatenate((node_xy, np.delete(node_all_xy[area], idx_com_coord_additional, axis=0)))
                
                elem_xy = np.concatenate((elem_xy, [[idx_additional_node[e-1] for e in elem] for elem in elem_all_xy[area]]), axis=0)
                num_node_xy = num_node_xy + num_node_all_xy[area] - len(com_coord_area)

        node_z = get_gep(0., hz1, 1., nGz1)
        num_node_z = len(node_z)

        self.node = [[node_xy[idx_xy][0]*0.001, node_xy[idx_xy][1]*0.001, node_z[idx_z]*0.001] for idx_z in range(num_node_z) for idx_xy in range(num_node_xy)]
        num_node = len(self.node)
        self.num_node = num_node
        self.num_node_xy = num_node_xy
        self.node_xy = node_xy

        num_elem_xy = len(elem_xy)
        self.elem = [np.concatenate([np.array(elem_xy[j]) + num_node_xy*(i-1), np.array(elem_xy[j])+num_node_xy*i]) for i in range(1, nGz1+1) for j in range(num_elem_xy)]
        num_elem = len(self.elem)
        self.num_elem = num_elem

        self.node_dat = [[idx+1] + list(self.node[idx]) for idx in range(num_node)]
        self.elem_dat = [[idx+1] + list(self.elem[idx]) for idx in range(num_elem)]

        self.edgenode_dat = []
        self.edgenode = edgenode
        for en in edgenode:
            self.edgenode_dat.append(en*0.001)

    def generate(self, step):
        node_dat = self.node_dat
        num_node = self.num_node
        elem_dat = self.elem_dat
        num_elem = self.num_elem
        edgenode_dat = self.edgenode_dat

        str_step = step2str(step)
        user_name = sim_params.USER_NAME
        repo_name = sim_params.REPO_NAME
        with open("node.g.dat", mode="w") as f:
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

        with open("elem.g.dat", mode="w") as f:
            f.write(f"{num_elem}\t8")
            tab = "\t"
            for line in elem_dat:
                f.write(f"\n{tab.join(list(map(str, map(int, line))))}")
        f.close()

        with open("edge.dat", mode="w") as f:
            f.write(f"8")
            for line in edgenode_dat:
                f.write(f"\n{tab.join(list(map(str, map(float, line))))}")
        f.close()
