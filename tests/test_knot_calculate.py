import numpy as np
from pyarcknot.knot_calculate import (dec_to_bin, count_loops, turaev_genus,
                                   kauffman_data, writhe)


def test_dec_to_bin():
    assert (dec_to_bin(1023, 12) ==
            [0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]).all()


def test_count_loops():
    assert (count_loops(np.array(
        [[0, 0, 0, 0, 0, 1, 0, 2],
         [0, 0, 1, 0, 0, 4, 2, 0],
         [1, 0, 4, 0, 2, 0, 0, 0],
         [0, 1, 4, 0, 4, 2, 0, 0],
         [0, 0, 0, 2, 4, 0, 4, 1],
         [0, 0, 2, 4, 1, 0, 0, 0],
         [2, 4, 0, 1, 0, 0, 0, 0],
         [0, 2, 0, 0, 0, 0, 1, 0]],
        dtype=np.int8)) == 5)


def test_turaev_genus():
    assert (turaev_genus(np.array(
        [[0, 0, 0, 0, 0, 1, 0, 2],
         [0, 0, 1, 0, 0, 3, 2, 0],
         [1, 0, 3, 0, 2, 0, 0, 0],
         [0, 1, 3, 0, 3, 2, 0, 0],
         [0, 0, 0, 2, 3, 0, 3, 1],
         [0, 0, 2, 3, 1, 0, 0, 0],
         [2, 3, 0, 1, 0, 0, 0, 0],
         [0, 2, 0, 0, 0, 0, 1, 0]],
        dtype=np.int8)) == 2)


def test_kauffman_data():
    assert (kauffman_data(np.array(
        [[2, 0, 0, 1, 0],
         [0, 0, 1, 3, 2],
         [0, 1, 3, 2, 0],
         [1, 3, 2, 0, 0],
         [0, 2, 0, 0, 1]],
        dtype=np.int8)) ==
            [[3, 3],
             [2, 2],
             [2, 2],
             [1, 1],
             [2, 2],
             [1, 1],
             [1, 1],
             [0, 2]]
            ).all()


def test_writhe():
    assert (writhe(np.array(
        [[0, 0, 0, 0, 0, 1, 0, 2],
         [0, 0, 1, 0, 0, 3, 2, 0],
         [1, 0, 3, 0, 2, 0, 0, 0],
         [0, 1, 3, 0, 3, 2, 0, 0],
         [0, 0, 0, 2, 3, 0, 3, 1],
         [0, 0, 2, 3, 1, 0, 0, 0],
         [2, 3, 0, 1, 0, 0, 0, 0],
         [0, 2, 0, 0, 0, 0, 1, 0]],
        dtype=np.int8)) == -4)
