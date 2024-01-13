import numpy as np
import matplotlib.pyplot as plt

from matplotlib.path import Path
from matplotlib.patches import PathPatch

mpl_colors = ['blue', 'green', 'red', 'orange', 'purple',
              'cyan', 'olive', 'brown', 'pink', 'gray']

from pyarcknot.knot_matrix import xabo_xco
from pyarcknot.knot_calculate import turaev_genus, jones_polynomial


# draw arc diagram
def draw_arc(xco):
    s = np.shape(xco)[0]
    vertices = np.zeros([1, 2])

    sep = .25  # separation horizontal/vertical gap in cross

    for i in range(s):
        row, col = xco[i, :], xco[:, i]

        hverts, n_hlines = np.where(row > 0)[0], np.count_nonzero(row > 0) - 1
        for j in range(n_hlines):
            c0, c1 = hverts[j:j + 2]
            if row[c0] not in [0, 1, 2]: c0 += sep
            if row[c1] not in [0, 1, 2]: c1 -= sep
            vertices = np.concatenate((vertices,
                                       [[c0, s - i]], [[c1, s - i]]), axis=0)

        vverts, n_vlines = np.where(col > 0)[0], np.count_nonzero(col > 0) - 1
        for j in range(n_vlines):
            r0, r1 = vverts[j:j + 2]
            vertices = np.concatenate((vertices,
                                       [[i, s - r0]], [[i, s - r1]]), axis=0)

    vertices = vertices[1:, :]

    lines = round(np.shape(vertices)[0] / 2)
    codes = [Path.MOVETO] + [Path.LINETO]
    codes *= lines

    path = Path(vertices, codes)
    pathpatch = PathPatch(path, facecolor='none', edgecolor='black', lw=3)
    fig, ax = plt.subplots()
    ax.add_patch(pathpatch)

    ax.axis(False)
    ax.autoscale_view(True)
    ax.set_aspect('equal', 'box')
    plt.show()


# draw modified arc diagram
def draw_mod_arc(xabo):
    s = np.shape(xabo)[0]
    vertices = np.zeros([1, 2])

    sep = .25  # separation horizontal/vertical gap in cross

    for i in range(s):
        row, col = xabo[i, :], xabo[:, i]

        hverts, n_hlines = np.where(row > 0)[0], np.count_nonzero(row > 0) - 1
        for j in range(n_hlines):
            c0, c1 = hverts[j:j + 2]
            if row[c0] == 4: c0 += sep
            if row[c1] == 4: c1 -= sep
            vertices = np.concatenate((vertices,
                                       [[c0, s - i]], [[c1, s - i]]), axis=0)

        vverts, n_vlines = np.where(col > 0)[0], np.count_nonzero(col > 0) - 1
        for j in range(n_vlines):
            r0, r1 = vverts[j:j + 2]
            if col[r0] == 5: r0 += sep
            if col[r1] == 5: r1 -= sep
            vertices = np.concatenate((vertices,
                                       [[i, s - r0]], [[i, s - r1]]), axis=0)

    vertices = vertices[1:, :]

    lines = round(np.shape(vertices)[0] / 2)
    codes = [Path.MOVETO] + [Path.LINETO]
    codes *= lines

    path = Path(vertices, codes)
    pathpatch = PathPatch(path, facecolor='none', edgecolor='black', lw=3)
    fig, ax = plt.subplots()
    ax.add_patch(pathpatch)

    ax.axis(False)
    ax.autoscale_view(True)
    ax.set_aspect('equal', 'box')
    plt.show()


# get loop coords and loop codes
def loop_crd_cds(xabo):
    size = np.shape(xabo)[0]
    crossings = np.array(np.where(xabo > 3))

    n_crossings = np.shape(crossings)[1]
    cr_check = np.zeros((n_crossings * 2), dtype=np.int8)  # cross_check

    # 4 = A -> (up-right)(down-left)
    # 5 = B -> (up-left)(down-right)
    # cr_check = [c1u, c1d, c2u, c2d, ..., cnu, cnd]

    # CODES 0=vertex ; 1=crossing
    # p = position = [[coords in matrix],
    #                 is_horizontal_bool, is_increasing_bool]

    loop_coords, loop_codes = [], []
    n_loop = 0

    while np.any(cr_check == 0):
        loop_coords.append([])
        loop_codes.append([])
        fz = np.where(cr_check == 0)[0][0]  # first_zero
        p = [[crossings[:, int(fz / 2)][0],
              crossings[:, int(fz / 2)][1]], False, fz % 2 != 0]
        p0 = p[:]
        while True:
            loop_coords[n_loop].append(p[0])
            for i in range(1, size):
                r, c = 0, 0
                if not (p[2]): i = -i

                if p[1]:
                    c = i
                else:
                    r = i

                px = xabo[p[0][0] + r, p[0][1] + c]
                if px != 0:
                    # update p0 & p1
                    p[0], p[1] = [p[0][0] + r, p[0][1] + c], not (p[1])
                    break

                    # update p2 & codes
            if px == 1 or px == 2:
                line = xabo[p[0][0], :] if p[1] else xabo[:, p[0][1]]
                coord = p[0][1] if p[1] else p[0][0]
                p[2] = False if any(line[:coord] != 0) else True
                loop_codes[n_loop].append(0)
            elif px == 5:
                p[2] = not (p[2])

            # update cr_check & codes
            if px == 4 or px == 5:
                loop_codes[n_loop].append(1)
                for i in range(n_crossings):
                    if all(p[0] == crossings[:, i]):
                        if not (p[1]) or px == 5:
                            if p[2]:
                                cr_check[2 * i + 1] = 1
                            else:
                                cr_check[2 * i] = 1
                        else:
                            if p[2]:
                                cr_check[2 * i] = 1
                            else:
                                cr_check[2 * i + 1] = 1
                        break

            if np.all(p == p0): break

        n_loop += 1

    for i in range(n_loop):
        loop_codes[i] = [loop_codes[i][-1]] + loop_codes[i][:-1]

    return loop_coords, loop_codes


# draw loops
def draw_loops(xabo):
    radius = .5
    size = np.shape(xabo)[0]
    loop_coords, loop_codes = loop_crd_cds(xabo)

    lp_crd_path = loop_coords.copy()
    lp_cds_path = loop_codes.copy()

    n_loops = len(loop_codes)
    for i in range(n_loops):
        n_inserts = 0
        for j in range(len(loop_codes[i])):
            if loop_codes[i][j] != 0:
                i_p = j + n_inserts
                ins1, ins2 = [0, 0], [0, 0]
                len_lcpath = len(lp_crd_path[i])
                for k in range(2):
                    c0 = lp_crd_path[i][i_p - 1][k]
                    c1 = lp_crd_path[i][i_p][k]
                    c2 = lp_crd_path[i][(i_p + 1) % len_lcpath][k]
                    if c0 == c1:
                        ins1[k] = c1
                    elif c0 < c1:
                        ins1[k] = c1 - radius
                    elif c0 > c1:
                        ins1[k] = c1 + radius
                    if c2 == c1:
                        ins2[k] = c1
                    elif c2 < c1:
                        ins2[k] = c1 - radius
                    elif c2 > c1:
                        ins2[k] = c1 + radius

                lp_crd_path[i].insert(i_p, ins1)
                lp_crd_path[i].insert(i_p + 2, ins2)
                n_inserts += 2

    for i in range(n_loops):
        lp_cds_path[i] = [Path.MOVETO]
        for code in loop_codes[i]:
            if code == 0:
                lp_cds_path[i] += [Path.LINETO]
            else:
                lp_cds_path[i] += [Path.CURVE3] * 2 + [Path.LINETO]
        lp_cds_path[i].pop(-1)
        lp_cds_path[i] += [Path.CLOSEPOLY]

    for loop in lp_crd_path:
        for coord in loop:
            coord[0], coord[1] = coord[1], size - (coord[0] + 1)
        loop += [[0, 0]]

    fig, ax = plt.subplots()

    for i in range(n_loops):
        path = Path(lp_crd_path[i], lp_cds_path[i])
        pathpatch = PathPatch(path, facecolor='none',
                              edgecolor=mpl_colors[i % len(mpl_colors)], lw=3)
        ax.add_patch(pathpatch)

    ax.axis(False)
    ax.autoscale_view(True)
    ax.set_aspect('equal', 'box')
    plt.show()

    return n_loops


# draw arc diagram, sA diagram & sB diagram
def draw_diagrams(xco, t_genus=False, jones_poly=False):
    c = np.shape(np.array(np.where(xco == 3)))[1]
    print('Number of crossings =', c)
    if t_genus:
        print('Turaev genus =', int(turaev_genus(xco)))
    if jones_poly:
        print('Jones Polynomial =', jones_polynomial(xco), '\n')
    print()

    xabo_a, xabo_b = xabo_xco(xco, np.ones(c) * 4), xabo_xco(xco, np.ones(c) * 5)
    print('|s_A D| =', draw_loops(xabo_a), '\n')
    draw_arc(xco)
    print('\n')
    print('|s_B D| =', draw_loops(xabo_b), '\n')
