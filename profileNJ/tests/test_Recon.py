# This file is part of profileNJ testing module

__author__ = "Emmanuel Noutahi"

import unittest
from ..TreeLib import ClusterUtils as C
from ..TreeLib import TreeClass
from ..TreeLib.TreeUtils import *
from ..TreeLib import params
from ..tests import dirname
import numpy as np
import os


genefile1 = os.path.join(dirname, "genetree/gtree1.nw")
genefile2 = os.path.join(dirname, "genetree/gtree2.nw")
specfile = os.path.join(dirname, "specietree/spec1.nw")


class TestReconciliation(unittest.TestCase):

    def setUp(self):
        self.gtree1 = TreeClass(genefile1)
        self.gtree2 = TreeClass(genefile2)
        self.stree = TreeClass(specfile)
        self.tree1_spec_map = {}
        for node in self.gtree1:
            self.tree1_spec_map[node.name] = node.name.split('_')[0]

    def test_set_species(self):
        self.gtree1.set_species(pos="prefix")
        for node in self.gtree1:
            self.assertEqual(self.tree1_spec_map[node.name], node.species)

    def test_retrict_to_spec(self):
        self.gtree1.set_species(pos="prefix")
        treecopy = self.gtree1.copy()
        treecopy.restrict_to_species(['dsim', 'dsec', 'dmel', 'dyak', 'dere'])
        restrict_tree = TreeClass(
            "(((dsim_10, dsec_9), dmel_7), (dyak_13, dere_12));")
        rf = restrict_tree.robinson_foulds(treecopy)[0]
        assert(rf == 0)

    def test_compute_dl(self):
        # tree1
        self.gtree1.set_species(pos="prefix")
        stree1 = self.stree.copy()
        stree1.prune(self.gtree1.get_leaf_species())
        lcamap1 = lcaMapping(self.gtree1, stree1)

        # tree 2
        self.gtree2.set_species(pos="prefix")
        stree2 = self.stree.copy()
        stree2.prune(self.gtree2.get_leaf_species())
        lcamap2 = lcaMapping(self.gtree2, stree2)

        nd1, nl1 = computeDL(self.gtree1)
        nd2, nl2 = computeDL(self.gtree2)
        assert nd1 == 0 and nl1 == 0 and nd2 == 6 and nl2 == 3

        # do not prune the specie tree to the genetree before
        lcamap2 = lcaMapping(self.gtree2, self.stree)
        lcamap1 = lcaMapping(self.gtree1, self.stree)
        d1, l1 = computeDL(self.gtree1)
        d2, l2 = computeDL(self.gtree2)
        assert d1 == 0 and l1 == 2 and d2 == 6 and l2 == 4

    def test_computeCost(self):
        import hashlib
        dupcost = {}
        losscost = {}
        set_to_half = ['dper', 'dwil']
        for n in self.stree.get_leaf_names():
            hashed_n = params.get_hash([n])
            dupcost[hashed_n] = 1
            losscost[hashed_n] = 1
            if n in set_to_half:
                losscost[hashed_n] = 1.5

        params.set(dupcost, losscost, internal_mode='mean')
        self.gtree1.set_species(pos="prefix")
        lcamap1 = lcaMapping(self.gtree1, self.stree)
        d, l = computeDLScore(self.gtree1)
        assert d == 0
        self.assertAlmostEqual(l, 2.75)
