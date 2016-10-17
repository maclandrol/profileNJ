# This file is part of profileNJ testing module

__author__ = "Emmanuel Noutahi"

import unittest
from ..TreeLib import ClusterUtils as C
from ..TreeLib import TreeClass
from ..tests import dirname

import numpy as np
import os

treefile = os.path.join(dirname, "genetree/tree1.nw")


class TestTreeClass(unittest.TestCase):

    def setUp(self):
        self.tree1 = TreeClass(treefile)
        self.tree2 = TreeClass("(((a,b)e,c)f,d)h;", format=1)

    def test_copy(self):
        treecopy = self.tree1.copy("simplecopy")
        rf = treecopy.robinson_foulds(self.tree1, unrooted_trees=True)
        assert rf[0] == 0
        for node in treecopy:
            ori_node = self.tree1 & node.name
            for feat in ori_node.features:
                assert getattr(node, feat) == getattr(ori_node, feat)

    def test_get_child(self):
        node = self.tree2.get_children()[0]
        with self.assertRaises(IndexError):
            child = node.get_child_at(2)

    def test_has_ancestor(self):
        node = self.tree2 & 'a'
        self.assertTrue(node.has_ancestor(['e', 'f', 'h']))
        self.assertFalse(node.has_ancestor(['e', 'd']))

    def test_has_descendant(self):
        node = self.tree2 & 'f'
        self.assertTrue(self.tree2.has_descendant(['e', 'b', 'd']))
        self.assertFalse(node.has_descendant(['d']))

    def test_euler_visit(self):
        euler_tour = self.tree2._euler_visit([])
        example_tour = "hfeaebefcfhdh"
        euler_tour = "".join([node.name for node in euler_tour] * 2)
        self.assertTrue(example_tour in euler_tour)

    def test_get_nodes_between(self):
        node1 = self.tree2 & 'a'
        node2 = self.tree2 & 'c'
        t_nodes_between = 'ef'
        computed_nodes_between = "".join(
            [n.name for n in node1.get_nodes_between(node2)])
        self.assertEqual(t_nodes_between, computed_nodes_between)

    def test_to_polytomy(self):
        polytomy1 = TreeClass("((a,b)e, c, d)h;", format=1)
        polytomy2 = TreeClass("(a,b,c,d)h;", format=1)
        t1 = self.tree2.copy()
        t2 = self.tree2.copy()
        (t1 & 'f').to_polytomy()
        (t2 & 'f').to_star()
        rf1 = polytomy1.robinson_foulds(t1, unrooted_trees=True)[0]
        rf2 = polytomy2.robinson_foulds(t2, unrooted_trees=True)[0]
        assert rf1 == rf2 == 0

    def test_is_binary(self):
        polytomy1 = TreeClass("((a,b)e, c, d)h;", format=1)
        for node in self.tree2.iter_internal_node():
            self.assertTrue(node.is_binary())
        self.assertFalse(polytomy1.is_binary())

    def test_is_internal(self):
        self.assertTrue((self.tree2 & 'f').is_internal())
        self.assertFalse(self.tree2.is_internal())
        self.assertFalse((self.tree2 & 'a').is_internal())

    def test_reroot(self):
        rerooted = list(self.tree2.reroot())
        edge_reroot = list(self.tree2.edge_reroot())
        for i in xrange(len(rerooted)):
            for j in xrange(i, len(rerooted)):
                rf = rerooted[i].robinson_foulds(
                    rerooted[j], unrooted_trees=True)[0]
                assert rf == 0
        for i in xrange(len(edge_reroot)):
            for j in xrange(i, len(edge_reroot)):
                rf = edge_reroot[i].robinson_foulds(
                    edge_reroot[j], unrooted_trees=True)[0]
                assert rf == 0
