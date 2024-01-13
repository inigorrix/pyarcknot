import numpy as np
from numba import njit, int8
import re


# cleans str and returns an arc array
def clean_k_arc(k_str):
    matches = re.compile(r'\d+').findall(k_str)
    nums = np.array([int(match) for match in matches])
    s = int(nums.shape[0] / 2)
    arc = nums.reshape(s, 2)
    while len(np.unique(arc[:, 0])) < s:
        for i in range(s):
            if np.count_nonzero(arc[:, 0] == arc[i, 0]) > 1:
                arc[i] = np.flip(arc[i])
    return arc


# turns an arc array into a xo array
def xo_arc(arc):
    s = np.shape(arc)[0]
    xo = np.zeros([s, s], dtype=np.int8)
    for i in range(s):
        for j in range(2):
            ik = i + 1
            ij = arc[i, j] - 1
            xo[-ik, ij] = j + 1
    return xo


# turns a xo array into a xco array
def xco_xo(xo):
    xco = xo.astype(np.int8)
    s = np.shape(xo)[0]
    r = range(1, s - 1)
    for i in r:
        for j in r:
            if xo[i, j] == 0:
                hlz = np.count_nonzero(xo[i, :j])
                hrz = np.count_nonzero(xo[i, j + 1:])
                vuz = np.count_nonzero(xo[:i, j])
                vdz = np.count_nonzero(xo[i + 1:, j])
                if hlz == 1 and hrz == 1 and vuz == 1 and vdz == 1:
                    xco[i, j] = 3
    return xco


# turns an arc array into a xco array
def xco_arc(arc):
    return (xco_xo(xo_arc(arc)))


# turns a xo or xco array into an arc array
def arc_xoco(xoco):
    s = np.shape(xoco)[0]
    arc = np.zeros([s, 2], dtype=np.int32)
    for i in range(s):
        ik = i + 1
        j1 = int(np.where(xoco[-ik] == 1)[0]) + 1
        j2 = int(np.where(xoco[-ik] == 2)[0]) + 1
        arc[i, 0], arc[i, 1] = j1, j2
    return arc


# turns a xco array into a a xabo array
def xabo_xco(xco, ab_crossings):
    xabo = xco.copy()
    xabo[xabo == 3] = ab_crossings
    return xabo


# turns a xco array into a a xabo array
@njit(int8[:, :](int8[:, :], int8[:]))
def xABo_xco(xco, ab_crossings):
    xabo = np.copy(xco)
    nc = np.count_nonzero(xabo == 3)
    crossings = np.zeros((2, nc), dtype=np.int8)
    crossings[0], crossings[1] = np.where(xabo == 3)[0], np.where(xabo == 3)[1]
    for i in range(nc):
        xabo[crossings[0, i], crossings[1, i]] = ab_crossings[i]
    return xabo
