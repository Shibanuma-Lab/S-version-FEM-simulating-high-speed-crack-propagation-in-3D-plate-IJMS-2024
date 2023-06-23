import global_mesh

def judge(x, y, z, edgenode) -> bool:
        judge = False
        ncross = 0
        x = x * 1000
        y = y * 1000
        num = len(edgenode)
        for i in range(num):
            if i + 1 >= num:
                n, m = num-1, 0
            else:
                n, m = i-1, i
            en_n = edgenode[n][1] 
            en_m = edgenode[m][1]
            if en_n == en_m:
                continue
            else:
                if y < min(en_n, en_m):
                    continue
                else:
                    if max(en_n, en_m) <= y:
                        continue
                    else:
                        xx = edgenode[n][0] + (y - en_n) * (edgenode[m][0] - edgenode[n][0]) / (en_m - en_n)
                    if abs(y-edgenode[n][1]) < 1e-8:
                        xx = edgenode[n][0]
                    if edgenode[m][0] == edgenode[n][0]:
                        xx = edgenode[m][0]
            if xx > x:
                ncross += 1

        if ncross % 2 == 1:
            judge = True
        if z > 0.015:
            judge = True

        return judge