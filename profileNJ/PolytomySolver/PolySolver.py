"""
Author: Manuel Lafond
Date: 07/2015
Functions for the the dynamic programming scheme of polytomy resolution
"""
import os
import sys

import copy

from ..TreeLib import *
from ..TreeLib import TreeUtils, TreeClass


class PolytomySolver:

    """
    Usage example :
    lcamap = TreeUtils.lcaMapping(genetree, specietree, multspeciename=False)

    ps = PolytomySolver(genetree, specietree, lcamap)
    ps.computeCostsTable()

    print "COST=", ps.getTableValue(specietree, 1)

    r = ps.getResolution()
    print r
    """

    def __init__(self, polytomy, speciestree, lcaMapping):
        self.polytomy = polytomy
        self.speciestree = speciestree
        self.lcaMapping = lcaMapping
        self.multiplicities = {}
        self.cup_values = {}
        self.debug = False
        self.dupcost = 1
        self.losscost = 1
        self.use_dp = False
        self.dp_values = {}
        self.special_species_dupcost = {}
        self.special_species_losscost = {}

    def setDupLossCosts(self, dupcost, losscost):
        self.dupcost = dupcost
        self.losscost = losscost

    def setSpecialDupLossCosts(self, species_node, dupcost, losscost):
        self.special_species_dupcost[species_node] = dupcost
        self.special_species_losscost[species_node] = losscost

    def computeMultiplicities(self):
        """ Compute multiplicities for each species s, i.e. the number of children x of self.polytomy such that lcaMapping[x] = s
        """

        self.multiplicities = {}
        for g in self.polytomy.get_children():
            s = self.lcaMapping[g]

            if s not in self.multiplicities:
                self.multiplicities[s] = 0

            self.multiplicities[s] += 1

    def computeCostsTable(self):
        """ Compute costs table, in time O(S)
        """
        if self.dupcost != self.losscost:
            self.use_dp = True

        # dp stands for dynamic-programming, which is used when DL costs are unequal
        # print "new call"
        # print  len(self.polytomy.get_children())
        # print len(self.polytomy.children)
        # print self.polytomy
        # EN changed this
        dpsize = len(self.polytomy.get_children()) + 2
        dp_values = {}

        self.computeMultiplicities()

        # for a species s, cup_values[s] = ( bottom plateau value,  breakpt
        # left, breakpt right )
        self.cup_values = {}

        for s in self.speciestree.traverse("postorder"):

            s_dupcost = self.dupcost
            if s in self.special_species_dupcost:
                s_dupcost = self.special_species_dupcost[s]

            s_losscost = self.losscost
            if s in self.special_species_losscost:
                s_losscost = self.special_species_losscost[s]

            mult = 0
            if s in self.multiplicities:
                mult = self.multiplicities[s]

            if(s.is_leaf()):

                if mult == 0:
                    self.cup_values[s] = (s_losscost, 1, 1)
                else:
                    self.cup_values[s] = (0, mult, mult)

                # ----------------------------------
                # DP MODE
                if self.use_dp:
                    dp_values[s] = dpsize * [0]
                    for k in range(1, mult):
                        dp_values[s][k] = (mult - k) * s_dupcost
                    for k in range(mult + 1, dpsize):
                        dp_values[s][k] = (k - mult) * s_losscost
                # ----------------------------------

                if self.debug:
                    print("CUP=", self.cup_values[s])
                    if self.use_dp:
                        print("DP=", dp_values[s])
            else:
                sl = s.get_children()[0]
                sr = s.get_children()[1]

                cl = self.cup_values[sl]
                cr = self.cup_values[sr]

                if sl in self.special_species_losscost:
                    sploss_l = self.special_species_losscost[sl]
                else:
                    sploss_l = self.losscost

                if sr in self.special_species_losscost:
                    sploss_r = self.special_species_losscost[sr]
                else:
                    sploss_r = self.losscost

                # ---------------------------------------------------
                # SPECIAL CASE : we are NOT in use_dp mode, but some child has
                # a special cost
                if (not self.use_dp) and sploss_l > self.losscost:

                    self.cup_values[s] = (
                        self.getTableValue(sr, 1) + sploss_l, 1, 1)

                elif (not self.use_dp) and sploss_r > self.losscost:

                    self.cup_values[s] = (
                        self.getTableValue(sl, 1) + sploss_r, 1, 1)

                # ---------------------------------------------------
                else:

                    lmin = cl[0]
                    l1 = cl[1]
                    l2 = cl[2]

                    rmin = cr[0]
                    r1 = cr[1]
                    r2 = cr[2]

                    # here we go, magic
                    if l1 < r1 and l2 < r1:
                        smin = lmin + rmin + r1 - l2
                        b1 = l2
                        b2 = r1
                    elif l1 < r1 and r1 <= l2 and l2 <= r2:
                        smin = lmin + rmin
                        b1 = r1
                        b2 = l2
                    elif l1 < r1 and l2 > r2:
                        smin = lmin + rmin
                        b1 = r1
                        b2 = r2
                    elif r1 <= l1 and l1 <= r2 and r1 <= l2 and l2 <= r2:
                        smin = lmin + rmin
                        b1 = l1
                        b2 = l2
                    elif r1 <= l1 and l1 <= r2 and l2 > r2:
                        smin = lmin + rmin
                        b1 = l1
                        b2 = r2
                    elif l1 > r2 and l2 > r2:
                        smin = lmin + rmin + l1 - r2
                        b1 = r2
                        b2 = l1
                    else:
                        raise ValueError("CASE NOT COVERED")

                    b1 += mult  # todo : verify this
                    b2 += mult

                    self.cup_values[s] = (smin, b1, b2)

                    # ---------------------------------------------------
                    # DP MODE

                    MAXINT = sys.maxint
                    if self.use_dp:
                        dp_values[s] = dpsize * [MAXINT]
                        bb1 = MAXINT
                        bb2 = MAXINT
                        bmin = MAXINT
                        s1 = s.get_children()[0]
                        s2 = s.get_children()[1]

                        # first pass : compute values from above (when they
                        # exist)
                        for k in range(1, dpsize):

                            if k >= mult + 1:
                                dp_values[s][k] = dp_values[s1][
                                    k - mult] + dp_values[s2][k - mult]
                                if dp_values[s][k] < bmin:
                                    bmin = dp_values[s][k]
                            else:
                                dp_values[s][k] = MAXINT

                        # second pass : locate the minimums on the row -> bb1
                        # and bb2 are their positions
                        for k in range(1, dpsize):
                            if bb1 == MAXINT:
                                if dp_values[s][k] == bmin:
                                    bb1 = k
                            elif bb1 != MAXINT and bb2 == MAXINT:
                                if dp_values[s][k] != bmin:
                                    bb2 = k - 1

                        # third pass : flatten
                        for k in range(1, bb1):
                            dp_values[s][
                                bb1 - k] = min(dp_values[s][bb1 - k], dp_values[s][bb1 - k + 1] + s_dupcost)
                        for k in range(bb2 + 1, dpsize):
                            dp_values[s][k] = min(dp_values[s][k], dp_values[
                                                  s][k - 1] + s_losscost)
                    # ---------------------------------------------------
                # end if, else of special case
                if self.debug:
                    print("CUP=", self.cup_values[s])
                    if self.use_dp:
                        print("DP=", dp_values[s])

        self.dp_values = dp_values

    def getTableValue(self, s, k):
        """ Returns the equivalent of table[s, k], assuming that computeCostsTable has been called previously
        """
        if self.use_dp:
            return self.dp_values[s][k]
        else:
            c = self.cup_values[s]  # format is a list with c = (min, b1, b2)

            if k < c[1]:
                return c[0] + c[1] - k
            elif k > c[2]:
                return c[0] + k - c[2]
            else:
                return c[0]

    def getAllResolutions(self, limit=100):
        """ Returns the first resolution found in Newick format, assuming that computeCostsTable has been called previously
        """
        r = self.getResolutions(self.speciestree, 1, limit)

        # let's hope this doesn't slow down
        return [row[0] for row in r]

    def getResolutions(self, s, k, limit=1):
        """ Returns all possible ways of having k subtrees rooted at s.
                Return value is an array of arrays.  Each entry of the main array (i.e. the first dimension)
                is a k-resolution, which are represented as arrays.
                These subarrays contain k elements in newick form, one for each subtree.
        """
        v = self.getTableValue(s, k)
        vright = self.getTableValue(s, k + 1)
        vleft = self.getTableValue(s, k - 1)

        s_dupcost = self.dupcost
        if s in self.special_species_dupcost:
            s_dupcost = self.special_species_dupcost[s]

        s_losscost = self.losscost
        if s in self.special_species_losscost:
            s_dupcost = self.special_species_losscost[s]

        # ------------------------------------------
        # the leaf case
        # ------------------------------------------
        if s.is_leaf():

            trees_to_return = []

            if k == 0:
                return [trees_to_return]

            if self.debug:
                print("At leaf s =", s.name, " k =",
                      k, "v =", v, " vright =", vright)

            if v == 0:
                trees_to_return = [s.name] * k
            # get k + 1 guys, merge 2
            # TODO : we just take one way of doing this
            elif v == vright + s_dupcost:

                # NOTE the [0] which assumes that there's only one possible k+1
                # resolution
                resz = self.getResolutions(s, k + 1, 1)[0]
                r1 = '(' + resz[0] + ',' + resz[1] + ')'
                del resz[0]
                resz[0] = r1
                trees_to_return = resz
            # get k - 1 guys, add a loss
            else:  # v == vleft + s_losscost
                # NOTE the [0] which assumes that there's only one k-1
                # resolution
                resz = self.getResolutions(s, k - 1, 1)[0]

                # resz.append(s.name + '_LOSS')

                trees_to_return = resz

            return [trees_to_return]
        # ------------------------------------------
        # the non-leaf case
        # ------------------------------------------
        else:

            if self.debug:
                print("At internal species s =", s.name, "k =", k)

            mult = 0
            if s in self.multiplicities:
                mult = self.multiplicities[s]

            s1 = s.get_children()[0]
            s2 = s.get_children()[1]
            vup1 = self.getTableValue(s1, k - mult)
            vup2 = self.getTableValue(s2, k - mult)

            all_solutions = []  # this is the return value we'll be filling, until we reach the limit

            # speciation path, we must make k joins.
            # Note that the speciation path is prioritized
            if k - mult > 0 and v == vup1 + vup2:

                all_resz_s1 = self.getResolutions(s1, k - mult, limit)
                all_resz_s2 = self.getResolutions(s2, k - mult, limit)

                for s1_iter in range(0, len(all_resz_s1)):
                    resz_s1 = all_resz_s1[s1_iter]

                    for s2_iter in range(0, len(all_resz_s2)):
                        resz_s2 = all_resz_s2[s2_iter]

                        if limit > 0:

                            # TODO : here we should also make each s1 - s2 possible combination
                            # but for now, we just make the obvious joins
                            resz = []

                            for i in range(0, k):
                                if i <= k - mult - 1:
                                    # this can happen when losses were inserted
                                    if len(resz_s1) > i and len(resz_s2) > i:
                                        res = '(' + resz_s1[i] + \
                                            ',' + resz_s2[i] + ')'
                                        resz.append(res)
                                    elif len(resz_s1) <= i:
                                        resz.append(resz_s2[i])
                                    elif len(resz_s2) <= i:
                                        resz.append(resz_s1[i])
                                else:
                                    # we go here when a child of the polytomy
                                    # was internal and had species s (mult > 0)
                                    resz.append(s.name)

                            all_solutions.append(resz)
                            limit -= 1

            # The duplication case.  We take the first two subtrees from the
            # right and join them
            if v == vright + s_dupcost and limit > 0:

                all_resz_dup = self.getResolutions(s, k + 1, limit)

                for resz_iter in range(0, len(all_resz_dup)):

                    if limit > 0:
                        resz = all_resz_dup[resz_iter]

                        r1 = '(' + resz[0] + ',' + resz[1] + ')'
                        del resz[0]
                        resz[0] = r1

                        limit -= 1
                        all_solutions.append(resz)

            # The loss case.  Add a loss in s
            if v == vleft + s_losscost and limit > 0:

                # print "k=", k, " v=",v, " vup1=", vup1, " vup2=", vup2, "
                # vright=", vright, " vleft=",vleft

                all_resz_loss = self.getResolutions(s, k - 1, limit)

                for resz_iter in range(0, len(all_resz_loss)):

                    if limit > 0:
                        resz = all_resz_loss[resz_iter]

                        # resz.append('LOSS')

                        limit -= 1
                        all_solutions.append(resz)

            return all_solutions


class GeneTreeSolver:
    """
    Usage example :
    gts = GeneTreeSolver(genetree, speciestree, lcamap)

    r = gts.solvePolytomies(10000)

    print "NBSOLS=", len(r)
    print r
    """
    # EN changed this
    # def __init__(self, genetree, speciestree, lcaMapping):

    def __init__(self, genetree, speciestree, lcaMapping, dupcost, losscost):
        self.genetree = genetree
        self.speciestree = speciestree
        self.lcaMapping = lcaMapping
        self.debug = False
        # EN changed this : set dupcost and losscost using args values
        self.dupcost = dupcost
        self.losscost = losscost
        self.use_dp = False
        self.solutions_per_gene = {}

    def labelInternalNodes(self, tree):
        """ Gives a name to every unlabeled internal node of tree.  Names are of the form [i],
                where i starts at 1 and is incremented at each name given.
        """
        cpt = 1
        for t in tree.traverse("postorder"):

            if t.name == "" or t.name == "NoName":
                t.name = "[" + str(cpt) + "]"
                cpt += 1

    def solvePolytomies(self, limit=100):
        """ Solves each polytomy, combines each solution and returns an array of newick strings.
        """
        self.labelInternalNodes(self.speciestree)

        self.solutions_per_gene = {}

        # TreeUtils.lcaMapping(self.genetree, self.speciestree)
        s_images = TreeUtils.getImageTreeNode(
            self.genetree, self.speciestree, self.lcaMapping)

        remainingLimit = limit

        for g in self.genetree.traverse("postorder"):

            if not g.is_leaf():

                stree = s_images[g]

                cpttmp = 1
                snames_map = {}
                special_cost_species = {}  # key = species, val = cost
                # print g
                # print "BEFORE=", stree.get_ascii(show_internal=True)

                for s in stree.traverse("postorder"):
                    snames_map[s.name] = s
                    if not s.is_leaf():

                        if len(s.get_children()) == 1:
                            snew = s.add_child()
                            snew.add_features(depth=s.depth + 1)

                        s1 = s.get_children()[0]
                        s2 = s.get_children()[1]

                        if s1.depth - s.depth > 1:
                            spcost = s1.depth - s.depth
                            sptmp = s.add_child()
                            sptmp.name = "[TMP" + str(cpttmp) + "]"
                            s.remove_child(s1)
                            sptmp.add_child(s1)
                            sptmp_leaf = sptmp.add_child()
                            sptmp_leaf.name = "[TMP" + str(cpttmp + 1) + "]"
                            special_cost_species[sptmp_leaf] = spcost - 1
                            cpttmp += 2

                        if s2.depth - s.depth > 1:  # copy pasta :(
                            spcost = s2.depth - s.depth
                            sptmp = s.add_child()
                            sptmp.name = "[TMP" + str(cpttmp) + "]"
                            s.remove_child(s2)
                            sptmp.add_child(s2)
                            sptmp_leaf = sptmp.add_child()
                            sptmp_leaf.name = "[TMP" + str(cpttmp + 1) + "]"
                            special_cost_species[sptmp_leaf] = spcost - 1
                            cpttmp += 2

                # print "AFTER=", stree.get_ascii(show_internal=True)
                # print
                # "------------------------------------------------------------------"

                # for s in stree.traverse("postorder"):
                #  if not s.is_leaf():
                #    for schild in s.get_children():
                #      if schild.depth - s.depth:
                tmpMap = {}
                for gchild in g.get_children():
                    tmpMap[gchild] = snames_map[self.lcaMapping[gchild].name]

                ps = PolytomySolver(g, stree, tmpMap)

                for special_sp in special_cost_species:
                    ps.setSpecialDupLossCosts(
                        special_sp, self.dupcost, self.losscost * special_cost_species[special_sp])

                ps.debug = self.debug
                ps.setDupLossCosts(self.dupcost, self.losscost)

                # or len(special_cost_species) > 0)  TODO : ML says THIS IS
                # WRONG !
                ps.use_dp = self.use_dp
                # EN add this
                # print ps.use_dp
                ps.computeCostsTable()

                if self.debug:
                    print "COST=", ps.getTableValue(self.speciestree, 1)

                # print remainingLimit

                r = ps.getAllResolutions(remainingLimit)

                # if len(g.get_children()) > 2:
                #  print "G CHILDREN=", len(g.get_children()), " RLEN=", len(r)
                #  xstr = ""
                #  for c in g.get_children():
                #    xstr += self.lcaMapping[c].name + "  "
                #  print xstr

                # TODO : limit is not applied so nicely here
                # TODO : handle the case ( (a,a,a), b, b ) -> here a,a,a makes
                # a leaf label
                for gchild in g.get_children():

                    # for each sol of r, we replace the i-th child by each possible solution
                    # More precisely, if rchild has species [x], then the first occurrence
                    # of substring [x] in each solution of r gets replaced
                    # In this manner, the internal children of r get replaced
                    # incrementally.
                    if not gchild.is_leaf():

                        r2 = []
                        for sol in r:

                            if limit > 0:
                                s = self.lcaMapping[gchild]
                                pos = sol.find(s.name)

                                rchild = self.solutions_per_gene[gchild]

                                for sol_child in rchild:

                                    # sol2 = sol[:pos] + "(" + sol_child + ")" + sol[pos + len(s.name):]
                                    sol2 = sol[:pos] + sol_child + \
                                        sol[pos + len(s.name):]
                                    r2.append(sol2)

                                    # limit -= 1

                        # limit = max(limit, 1)

                        r = r2

                if len(r) > limit:
                    r = r[0:limit]
                self.solutions_per_gene[g] = r
                remainingLimit = limit - len(r) + 2
                # remainingLimit = max(remainingLimit, 1)

        return self.solutions_per_gene[self.genetree]
