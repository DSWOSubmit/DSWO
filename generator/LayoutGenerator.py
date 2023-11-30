def generateGridLayout(n = 10, m = 10):
    phy_layout = []
    # n * m
    for i in range(0, m):
        for j in range(0, n):
            now = i * n + j
            if j + 1 < n:
                phy_layout.append((now, now + 1))
            if i + 1 < m:
                phy_layout.append((now, now + n))
    return n * m, phy_layout

def generateIBMQLayout(number = 15):
    # https://journals.aps.org/prx/pdf/10.1103/PhysRevX.10.011022
    count5 = (number - 1) // 4
    mod5 = number - count5 * 4 - 1
    numberOfLine = count5 * 2 + 1
    qubitid = 0
    edge = []
    for i in range(numberOfLine):
        countOfLine = count5 * 4 + 3
        if i == 0 or i + 1 == numberOfLine:
            countOfLine = number
        thisLine = []
        for j in range(countOfLine - 1):
            thisLine.append(qubitid)
            edge.append((qubitid, qubitid + 1))
            qubitid += 1
        thisLine.append(qubitid)
        if i + 1 == numberOfLine:
            break
        if i % 2 == 0:
            for j in range(count5 + 1):
                qubitid += 1
                edge.append((qubitid - countOfLine + 3 * j, qubitid))
                edge.append((qubitid, qubitid + count5 + 1 + 3 * j))
            qubitid += 1
        else:
            for j in range(count5 + 1):
                qubitid += 1
                edge.append((qubitid - countOfLine + 3 * j + 2, qubitid))
                if i + 2 == numberOfLine:
                    edge.append((qubitid, qubitid + 3 * j + count5 + 1 + mod5))
                else:
                    edge.append((qubitid, qubitid + count5 + 1 + 3 * j + 2))
            qubitid += 1
    return qubitid + 1, edge
