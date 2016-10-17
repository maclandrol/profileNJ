# This file is part of profileNJ testing module

__author__ = "Emmanuel Noutahi"

import unittest
from ..TreeLib import ClusterUtils as C
from ..TreeLib import TreeClass
from ..tests import dirname
import numpy as np
import os

distmatfilename1 = os.path.join(dirname, "distmat/distmat1.dist")


class TestCluster(unittest.TestCase):

    def setUp(self):
        self.distmat1 = np.array([[0, 5, 9, 9, 8], [5, 0, 10, 10, 9],
                                  [9, 10, 0, 8, 7], [9, 10, 8, 0, 3],
                                  [8, 9, 7, 3, 0]], dtype=float)
        self.upgmamat = np.array([[0, 20, 60, 100, 90], [20, 0, 50, 90, 80],
                                  [60, 50, 0, 40, 50], [100, 90, 40, 0, 30],
                                  [90, 80, 50, 30, 0]], dtype=float)
        self.condmat1 = np.array(
            [[0, 7, 7, 6], [7, 0, 8, 7], [7, 8, 0, 3], [6, 7, 3, 0]])

        self.qmat1 = np.array([[0, -50, -38, -34, -34], [-50, 0, -38, -34, -34],
                               [-38, -38, 0, -40, -40], [-34, -34, -40, 0, -48],
                               [-34, -34, -40, -48, 0]])
        self.node_names = ['a', 'b', 'c', 'd', 'e']
        self.first_smallest_ind = [0, 1]

    def test_distMatProcessor(self):
        loadmat1 = C.distMatProcessor(distmatfilename1)
        assert np.array_equal(self.distmat1, loadmat1[0])
        assert self.node_names == loadmat1[1]

    def test_q_matrix(self):
        tmpqmat1 = C.calculate_Q_matrix(self.distmat1)
        np.fill_diagonal(tmpqmat1, 0)
        assert np.array_equal(self.qmat1, tmpqmat1)

    def test_smallest_index(self):
        s = C.find_smallest_index(self.qmat1)
        assert sorted(s) == self.first_smallest_ind

    def test_paired_node_distance(self):
        i, j = sorted(self.first_smallest_ind)
        assert C.paired_node_distance(self.distmat1, (i, j)) == (2, 3)

    def test_condense_matrix(self):
        assert np.array_equal(self.condmat1, C.condense_matrix(
            self.distmat1, self.first_smallest_ind, 'nj'))

    def test_nj(self):
        node_order = [TreeClass("a:1;"), TreeClass("b:1;"), TreeClass(
            "c:1;"), TreeClass("d:1;"), TreeClass("e:1;")]
        t2_nj, final_array, smallest_index = C.treeCluster(
            self.distmat1, node_order, depth=None, method='nj')
        t1_nj = TreeClass("(((a,b),c),(d,e));")
        rf = t1_nj.robinson_foulds(t2_nj, unrooted_trees=True)
        assert rf[0] == 0

    def test_upgma(self):
        node_order = [TreeClass("a:1;"), TreeClass("b:1;"), TreeClass(
            "c:1;"), TreeClass("d:1;"), TreeClass("e:1;")]
        t2_upgma, final_array, smal_pos = C.treeCluster(
            self.upgmamat, node_order, depth=None, method='upgma')
        t1_upgma = TreeClass('((a,b),((d,e),c));')
        rf = t1_upgma.robinson_foulds(t2_upgma)
        assert rf[0] == 0
