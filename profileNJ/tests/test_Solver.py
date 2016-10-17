# This file is part of profileNJ testing module

__author__ = "Emmanuel Noutahi"

import unittest
from ..TreeLib import ClusterUtils as C
from ..TreeLib import TreeClass
from ..TreeLib.TreeUtils import *
from ..TreeLib import params
from ..PolytomySolver.Multipolysolver import polytomyPreprocess, polySolver
from ..PolytomySolver.Multipolysolver import solvePolytomy, computePolytomyReconCost
from ..PolytomySolver import *
from ..tests import dirname

import numpy as np
import os

genefile = os.path.join(dirname, "genetree/tree1.nw")
g2solve = os.path.join(dirname, "genetree/500.nw")
s2solve = os.path.join(dirname, "specietree/500.nw")
specfile = os.path.join(dirname, "specietree/spec2.nw")
ntgenefile = os.path.join(dirname, "genetree/notung_tree.nw")
distmat = os.path.join(dirname, "distmat/distmat2.dist")


class TestSolver(unittest.TestCase):

    def setUp(self):
        self.star = TreeClass("(a_1,a_2,b_1,b_2,b_3,c_1);")
        self.nostar = TreeClass(
            "((a_1, b_2), (d_1, d_2), a_2, b_1, b_3, c_1);")
        self.sptreestar = TreeClass("((c,d)e,(a,b)f)g;", format=1)
        self.gtree = TreeClass(genefile)
        self.stree = TreeClass(specfile)
        self.nttree = TreeClass(ntgenefile, format=1)
        self.matrix = {
            'a': np.asarray([0.5, 0, 1, 2], dtype=float),
            'b': np.asarray([1, 0.5, 0, 1], dtype=float),
            'c': np.asarray([0, 1, 2, 3], dtype=float),
            'd': np.asarray([1.5, 3, 4.5, 6], dtype=float),
            'e': np.asarray([1.5, 2.75, 4, 5.25], dtype=float),
            'f': np.asarray([1, 0.5, 1, 2], dtype=float),
            'g': np.asarray([2.5, 3.25, 4.375, 5.5], dtype=float)
        }

        self.nostarmat = {
            'a': np.asarray([0, 1, 2], dtype=float),
            'b': np.asarray([1, 0, 1], dtype=float),
            'c': np.asarray([0, 1, 2], dtype=float),
            'd': np.asarray([0, 1, 2], dtype=float),
            'f': np.asarray([2, 1, 1], dtype=float),
            'e': np.asarray([0, 1, 2], dtype=float),
            'g': np.asarray([2, 2, 3], dtype=float)
        }

        self.distance_matrix, self.node_order = C.distMatProcessor(distmat)
        self.star.set_species(pos="prefix")
        self.nostar.set_species(pos="prefix")

    def test_polysolver_star(self):
        dup_cost = {}
        loss_cost = {}
        for node in self.sptreestar:
            dup_cost[params.get_hash(node.name)] = 1
            loss_cost[params.get_hash(node.name)] = 1
            if node.name == 'd':
                loss_cost[params.get_hash(node.name)] = 1.5
            if node.name in ['a', 'b']:
                dup_cost[params.get_hash(node.name)] = 0.5
        params.set(dup_cost, loss_cost, internal_mode='mean')

        for node in self.sptreestar.traverse():
            print(node.name, ' L : ', params.getloss(
                node), ' D : ', params.getdup(node))
        matrix, row_node = polySolver(treeHash(
            self.star, addinfos='demo'), self.star, self.sptreestar, None, [], 1, verbose=False, mode="none")

        for key in row_node:
            assert np.array_equal(matrix[key, :], self.matrix[
                                  row_node[key].name])
        solutions = polySolver(treeHash(self.star), self.star, self.sptreestar,
                               self.distance_matrix, self.node_order, -1, verbose=False, cluster_method='nj')
        # only one solution
        assert len(solutions) == 1
        resolved_tree = solutions[0]
        rf = resolved_tree.robinson_foulds(
            TreeClass("(((a_2, b_3), ((b_1,b_2), a_1)), c_1);"))[0]
        assert rf == 0

    def test_polysolver_no_star(self):
        dup_cost = {}
        loss_cost = {}
        for node in self.sptreestar:
            dup_cost[params.get_hash(node.name)] = 1
            loss_cost[params.get_hash(node.name)] = 1
        params.set(dup_cost, loss_cost)

        node_order = self.nostar.get_leaf_names()
        dist_mat = C.makeFakeDstMatrice(len(node_order), 1, 10)
        dist_mat, node_order = polytomyPreprocess(
            self.nostar, self.sptreestar, dist_mat, node_order, method='nj')
        matrix, row_node = polySolver(treeHash(self.nostar, addinfos='demo'),
                                      self.nostar, self.sptreestar, None, [], 1, verbose=False, mode="none")
        for key in row_node:
            assert np.array_equal(matrix[key, :], self.nostarmat[
                                  row_node[key].name])

        solutions = polySolver(treeHash(self.nostar), self.nostar, self.sptreestar,
                               dist_mat, node_order, -1, verbose=False, cluster_method='nj')
        # only one solution
        assert len(solutions) == 1
        tree = solutions[0]
        true_recon_cost = self.nostarmat['g'][0]
        lcamap = lcaMapping(tree, self.sptreestar)
        dup, loss = computeDLScore(tree)
        # because we have one more duplication already added in the original
        # tree
        self.assertEqual(dup + loss, true_recon_cost + 1)
        self.assertEqual(
            dup + loss, computePolytomyReconCost(self.nostar, self.sptreestar))

    def test_solvepolytomy(self):
        self.gtree.set_species(pos="prefix")
        self.stree.label_internal_node()
        node_order = self.gtree.get_leaf_names()
        dist_mat = C.makeFakeDstMatrice(len(node_order), 1, 10)
        recon_score = computePolytomyReconCost(self.gtree, self.stree)
        notung_score = (6, 9)
        assert recon_score == notung_score[0] + notung_score[1]
        solutions = solvePolytomy(
            self.gtree, self.stree, dist_mat, node_order, False, -1, 'nj', -1)
        tree = solutions[0]
        lcamap = lcaMapping(tree, self.stree)
        dup, loss = computeDL(tree)
        assert (dup, loss) == notung_score
        # rf = tree.robinson_foulds(self.nttree)
        # print tree
        # print self.nttree
        # print rf[0]

    def test_polySolver(self):
        genetree = TreeClass(g2solve)
        genetree.set_species(use_fn=lambda x: x.name)
        specietree = TreeClass(s2solve)
        lcamap = lcaMapping(genetree, specietree, False)
        solver = PolySolver.GeneTreeSolver(genetree, specietree, lcamap, 1, 1)
        solver.labelInternalNodes(genetree)
        solver.labelInternalNodes(specietree)
        solver.use_dp = False
        self.solvePolytomies(solver, 'psolver')

    def test_notungSolver(self):
        genetree = TreeClass(g2solve)
        genetree.set_species(use_fn=lambda x: x.name)
        specietree = TreeClass(s2solve)
        specietree.label_internal_node()
        lcamap = lcaMapping(genetree, specietree, False)
        solver = SingleSolver.NotungSolver(genetree, specietree, lcamap)
        self.solvePolytomies(solver)

    def test_zheng_Solver(self):
        genetree = TreeClass(g2solve)
        genetree.set_species(use_fn=lambda x: x.name)
        specietree = TreeClass(s2solve)
        specietree.label_internal_node()
        lcamap = lcaMapping(genetree, specietree, False)
        solver = SingleSolver.LinPolySolver(genetree, specietree, lcamap)
        self.solvePolytomies(solver)

    def test_dynSolver(self):
        genetree = TreeClass(g2solve)
        genetree.set_species(use_fn=lambda x: x.name)
        specietree = TreeClass(s2solve)
        specietree.label_internal_node()
        lcamap = lcaMapping(genetree, specietree, False)
        solver = SingleSolver.Dynamiq2(genetree, specietree, lcamap)
        self.solvePolytomies(solver)

    def solvePolytomies(self, solver, solver_type='defsolver'):

        sols = []
        if solver_type == 'psolver':
            sols = [x + ";" for x in solver.solvePolytomies(1)]
        else:
            sols = [solver.reconstruct()]
        t = TreeClass(sols[0])
        # should have 750 leaves
        assert len(t) == 750
        t.set_species(use_fn=lambda x: x.name)
        specietree = TreeClass(s2solve)
        lcamap = lcaMapping(t, specietree, False)
        dup, loss = computeDL(t, lcamap)
        self.assertAlmostEqual((dup + loss) / 7605, 1)
