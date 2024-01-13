import numpy as np
from numba import njit, int8, int32
import sympy as smp
import matplotlib.pyplot as plt
from pyarcknot.knot_matrix import xabo_xco, xABo_xco


# Transform decimal numbers to binary
@njit(int8[:](int32, int32))
def dec_to_bin(n, p):
    n_bin = np.zeros(p, dtype=np.int8)
    for i in range(p):
        x = 2 ** (p - i - 1)
        if n >= x:
            n_bin[i] = 1
            n -= x
    return n_bin


# Better loop counter
@njit(int32(int8[:, :]))
def count_loops(xabo):
    n_crossings = np.count_nonzero(xabo > 3)
    crossings = np.zeros((2, n_crossings), dtype=np.int8)
    crossings[0], crossings[1] = np.where(xabo > 3)[0], np.where(xabo > 3)[1]
    cr_check = np.zeros(n_crossings * 2, dtype=np.int8)  # cross_check

    # 4 = A -> (up-right)(down-left)
    # 5 = B -> (up-left)(down-right)
    # cr_check = [c1u, c1d, c2u, c2d, ..., cnu, cnd]

    # CODES False=vertex ; True=crossing
    # pos = [coords in matrix]
    # horiz = is_horizontal_bool
    # incr = is_increasing_bool (down/right = True)

    pos = np.zeros(2, dtype=np.int8)
    n_loop = 0

    while not np.all(cr_check):
        fz = np.where(cr_check == 0)[0][0]  # first_zero
        pos[:] = [crossings[:, int(fz / 2)][0], crossings[:, int(fz / 2)][1]]
        horiz, incr = False, fz % 2 != 0
        pos0, horiz0, incr0 = np.copy(pos), horiz, incr

        while True:
            # find next point(number) in matrix
            if horiz:
                if incr:
                    col = np.argwhere(xabo[pos[0], pos[1]:] > 0)[1, 0] + pos[1]
                else:
                    col = np.argwhere(xabo[pos[0], :pos[1]] > 0)[-1, 0]
                pos[1] = col
            else:
                if incr:
                    row = np.argwhere(xabo[pos[0]:, pos[1]] > 0)[1, 0] + pos[0]
                else:
                    row = np.argwhere(xabo[:pos[0], pos[1]] > 0)[-1, 0]
                pos[0] = row

            horiz = not horiz
            px = xabo[pos[0], pos[1]]

            # update incr
            if px == 1 or px == 2:
                line = xabo[pos[0], :] if horiz else xabo[:, pos[1]]
                coord = pos[1] if horiz else pos[0]
                incr = True if np.all(line[:coord] == 0) else False
            elif px == 5:
                incr = not incr

            # update cr_check
            if px == 4 or px == 5:
                r, c = crossings[0, :] == pos[0], crossings[1, :] == pos[1]
                i = np.argwhere(np.logical_and(r, c))[0, 0]
                if not horiz or px == 5:
                    if incr:
                        cr_check[2 * i + 1] = 1
                    else:
                        cr_check[2 * i] = 1
                else:
                    if incr:
                        cr_check[2 * i] = 1
                    else:
                        cr_check[2 * i + 1] = 1

            if np.all(pos == pos0) and horiz == horiz0 and incr == incr0:
                break
        n_loop += 1

    return n_loop


# Genus of the Turaev Surface of the diagram
def turaev_genus(xco):
    n_crossings = np.count_nonzero(xco==3)
    sa, sb = count_loops(xabo_xco(xco, 4)), count_loops(xabo_xco(xco, 5))
    genus = (2 + n_crossings - sa - sb)*.5
    return genus


# Data to calculate Kauffman bracket
@njit(int32[:, :](int8[:, :]))
def kauffman_data(xco):
    n_crossings = np.count_nonzero(xco == 3)
    n_states = 2 ** n_crossings
    kauffman_data = np.zeros((n_states, 2), dtype=np.int32)
    ab_base = np.ones(n_crossings, dtype=np.int8) * np.int8(4)
    for i in range(n_states):
        ab_crossings = ab_base + dec_to_bin(i, n_crossings)
        kauffman_data[i, 0] = np.count_nonzero(ab_crossings == 4)
        xabo_state = xABo_xco(xco, ab_crossings)
        kauffman_data[i, 1] = count_loops(xabo_state)

    return kauffman_data


# Display Kauffman data
def kauffman_heatmap(xco, figsize=None, lps=5, labs=3, tics=2):
    kauf = kauffman_data(xco)
    k, c = np.unique(kauf, axis=0, return_counts=True)
    n = k[-1, 0] + 1
    k_mat = np.zeros((n, n), dtype=int)
    for i in range(len(c)):
        k_mat[k[i, 0], k[i, 1] - 1] = c[i]

    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(k_mat, cmap='GnBu')

    min_figsize = min(fig.get_size_inches())
    lp, lab_fs, tic_fs = lps * min_figsize, labs * min_figsize, tics * min_figsize

    ax.set_xticks(np.arange(n) + .5, minor=True)
    ax.set_yticks(np.arange(n) + .5, minor=True)
    ax.grid(which="minor", color="k", linestyle='-', linewidth=1)

    ax.tick_params(top=False, left=False, bottom=False,
                   labeltop=True, labelbottom=False)
    ax.set_xticks(np.arange(n), labels=np.arange(1, n + 1, 1), fontsize=tic_fs)
    ax.set_yticks(np.arange(n), labels=np.arange(n), fontsize=tic_fs)

    ax.xaxis.set_label_position('top')
    ax.set_xlabel('|sD|', labelpad=lp * .5, fontsize=lab_fs)
    ax.set_ylabel('a(s)', labelpad=lp, fontsize=lab_fs, rotation='horizontal')

    ax.set_title('States with a(s) & |sD|',
                 pad=lp * .5, fontsize=lab_fs * 1.5)

    sq_fs = 5 * min_figsize / len(str(max(c)))
    for i in range(n):
        for j in range(n):
            if k_mat[i, j] != 0:
                col = 'w' if k_mat[i, j] > .75 * max(c) else 'k'
                text = ax.text(j, i, k_mat[i, j], ha="center", va="center",
                               color=col, fontsize=sq_fs)
    plt.show()


# Kauffman bracket polynomial
def kauffman_bracket(xco):
    A, D = smp.symbols('A'), 0
    n_crossings = np.count_nonzero(xco == 3)
    k_d = kauffman_data(xco)
    na_nloops, count = np.unique(k_d, axis=0, return_counts=True)
    for i in range(len(count)):
        n_as, n_bs = na_nloops[i, 0], n_crossings - na_nloops[i, 0]
        n_loops_x = na_nloops[i, 1]
        D += count[i] * (A ** (n_as - n_bs) * (-A ** -2 - A ** 2) ** (n_loops_x - 1))
    return D.expand()


# Writhe
def writhe(xco):
    # 1 -> 2 in vertical, down=True
    # 2 -> 1 in horizontal, right=True
    c_coords = np.array(np.where(xco == 3))
    n_crossings = np.shape(c_coords)[1]
    w = 0
    for i in range(n_crossings):
        row, col = c_coords[0, i], c_coords[1, i]
        h = True if np.count_nonzero(xco[row, col:] == 1) == 1 else False
        v = True if np.count_nonzero(xco[row:, col] == 2) == 1 else False
        w += 1 if h == v else -1
    return w


# Jones polynomial
def jones_polynomial(xco):
    A, t = smp.symbols('A'), smp.symbols('t')
    D = kauffman_bracket(xco)
    w = writhe(xco)
    V_A = (-A) ** (-3 * w) * D
    V_t = V_A.xreplace({A: t ** (-1 / 4)})
    pol = str(V_t.expand()).replace('.0', '')
    if pol == '1':
        return 1
    else:
        return smp.Poly(pol).as_expr()
