from ..TreeLib import TreeClass, TreeUtils, params
from collections import defaultdict as ddict
import numpy as np
from functools import partial
import operator


def simpleConstruct(nodelist, nodemap, K, T, W, P, root):

    # nodelist ==> nl
    # K ==> ingene
    # T ==> outgene
    # W ==> multiplicity
    # Fid : Father id
    # image_nodes ==> P
    # print(nodelist)
    # parent = [(node.name, node.parent) for node in nodemap.values() if node.has_feature("parent")]
    # print(K)
    # print(T)
    # print(W)
    # print(parent)
    # print(P)

    def insertNode(node, nodes):
        if len(nodes) == 0:
            return
        tmp_node = TreeClass()
        tmp_node.replace_by(node)
        node.children = []
        for n in nodes[0:-1]:
            new_node = TreeClass()
            new_node.add_child(tmp_node)
            new_node.add_child(n)
            tmp_node = new_node

        node.add_child(tmp_node)
        node.add_child(nodes[-1])

    # store multiple copy of a node in a list
    nodes = ddict(partial(np.empty))

    for name in K.keys():
        nodes[name] = np.empty(int(K[name]), dtype=object)

    for name in reversed(nodelist):

        for j in xrange(K[name]):
            nodes[name][j] = TreeClass()

        node = nodemap[name]

        if W[name] > 0:
            if node.is_leaf() and K[name] != W[name] - 1:
                raise Exception("This shouldn't happen")
            elif not node.is_leaf() and K[name] < W[name]:
                raise Exception("This shouldn't happen either")

            if node.is_leaf():
                for j, n in enumerate(nodes[name]):
                    n.replace_by(P[name][j + 1])

                node.replace_by(P[name][0])

            else:
                for j, pnode in enumerate(P[name]):
                    nodes[name][K[name] - 1 - j].replace_by(pnode)

        if node.up is not None:
            j = 0
            while j < T[name] and j < K[name]:
                nodes[node.parent][j].add_child(nodes[name][j])
                j += 1

        if T[name] < K[name]:
            insertNode(node, nodes[name][T[name]:K[name]])

    for name, node in nodemap.items():
        if len(node.children) == 1:
            node.replace_by(node.children[0])

    nodemap[root].delete_single_child_internal()

    return nodemap[root]


class Solver(object):

    def __init__(self, genetree, specietree, lcamap):
        self.genetree = genetree
        self.specietree = specietree
        self.lcamap = lcamap

    def _compute_mult(self, node):
        W = ddict(int)
        rW = ddict(list)
        for g in node.get_children():
            s = self.lcamap[g]
            rW[s.name].append(g)
            W[s.name] += 1
        return W, rW

    def reconstruct(self):
        images_trees = TreeUtils.getImageTreeNode(
            self.genetree, self.specietree, self.lcamap)

        for node in self.genetree.traverse("postorder"):
            if (node.is_root() or node.is_internal()) and not node.is_binary():
                for n in images_trees[node].traverse():
                    if(n.up is not None):
                        n.add_features(parent=n.up.name)
                W, rW = self._compute_mult(node)
                solution = self.compute_table(images_trees[node], W, rW)
                node.replace_by(solution)
        return self.genetree.write(format=9)

    def get_solution(self, image_tree, multiplicities, reverse_node_map, outgene, ingene):
        """Reconstruct a genetree solution"""
        nodelist = [node.name for node in image_tree.traverse("postorder")]
        nodemap = dict((node.name, node)
                       for node in image_tree.traverse("postorder"))
        resolution = simpleConstruct(
            nodelist, nodemap, ingene, outgene, multiplicities, reverse_node_map, image_tree.name)
        # resolution.delete_single_child_internal()
        return resolution


class LinPolySolver(Solver):

    def __init__(self, genetree, specietree, lcamap):
        super(self.__class__, self).__init__(genetree, specietree, lcamap)

    def compute_table(self, image_tree, W, rW):
        """Compute cost table"""
        A = ddict(int)
        B = ddict(int)
        ingene = ddict(int)
        outgene = ddict(int)

        def update(a, b, degree_diff):
            if degree_diff > 2:
                return 0, 0
            elif degree_diff == 2:
                return 0, a
            elif degree_diff == 1:
                return a, b
            else:
                raise ValueError(
                    "Level difference incorrect, something is wrong")

        def getAB(node_mult, *args):
            args = sorted(args)
            return node_mult + args[1], node_mult + args[2]

        def get_med(i, a, b):
            if i <= a:
                return a
            elif i >= b:
                return b
            else:
                return i

        # first step (1) : traverse Image node in post-order
        for node in image_tree.traverse("postorder"):
            mult_node = W[node.name]
            a1, b1, a2, b2 = 0, 0, 0, 0

            # case where node is a leaf:
            if node.is_leaf():
                A[node.name] = mult_node - 1
                B[node.name] = A[node.name]

            else:
                # node has 2 children
                if(len(node.get_children()) == 2):
                    node2 = node.get_children()[1]
                    d2 = node2.depth - node.depth
                    a2, b2 = A[node2.name], B[node2.name]
                    a2, b2 = update(a2, b2, d2)
                    # case where node has only one child

                node1 = node.get_children()[0]
                a1, b1 = A[node1.name], B[node1.name]
                d1 = node1.depth - node.depth
                a1, b1 = update(a1, b1, d1)

                # A[node.name], B[node.name] = get_A_B(mult_node, a1, a2, b1, b2)
                A[node.name] = mult_node + min(max(a1, a2), min(b1, b2))
                B[node.name] = mult_node + max(max(a1, a2), min(b1, b2))

        # step 2 : traversing in preorder
        # in(u) and out(u) denotes the nos. of genes
        # flowing into and out of the branch (p(u), u)

        nil = list(image_tree.traverse("postorder"))
        outgene[nil[-1].name] = 0
        ingene[nil[-1].name] = get_med(0, A[nil[-1].name], B[nil[-1].name])
        for i in range(len(nil) - 2, -1, -1):

            node = nil[i]
            father = node.up
            t = ingene[father.name] - W[father.name]
            if node.depth - father.depth > 2:
                t = 0

            outgene[node.name] = t
            ingene[node.name] = get_med(t, A[node.name], B[node.name])

        return self.get_solution(image_tree, W, rW, outgene, ingene)


class DynPolySolver(Solver):

    def __init__(self, genetree, specietree, lcamap, dupcost=1, losscost=1):
        super(self.__class__, self).__init__(genetree, specietree, lcamap)
        self.dupcost = dupcost + 2 * losscost
        self.losscost = losscost

    def compute_table(self, image_tree, W, rW):
        """ Compute table for dynamique programmation"""
        U = ddict(partial(np.zeros))
        M = ddict(int)
        I = ddict(partial(np.zeros))
        # print( image_tree

        outgene = ddict(int)
        ingene = ddict(int)

        for node in image_tree.traverse("postorder"):
            mult_node = W[node.name]
            if node.is_leaf():
                M[node.name] = mult_node - 1
            elif len(node.get_children()) == 1:
                M[node.name] = M[node.get_child_at(0).name] + mult_node

            elif node.is_binary():
                v1 = M[node.get_child_at(0).name]
                v2 = M[node.get_child_at(1).name]
                M[node.name] = mult_node + max(v1, v2)

            else:
                raise ValueError(
                    "Shouldn't have a polytomy in the image tree, something wrong")

        if image_tree.is_root() and image_tree.is_leaf():
            outgene = ddict(int)
            ingene = ddict(int)
            outgene[image_tree.name] = 0
            ingene[image_tree.name] = W[image_tree.name] - 1
            return self.get_solution(image_tree, W, rW, outgene, ingene)

        mul = M[image_tree.name]
        C = np.zeros(mul + 1)

        def getCost(inp, out, d):
            m = min(float(d) * self.losscost, self.dupcost)
            if inp >= out:
                return m * float(out)
            else:
                return m * float(inp) + self.dupcost * float(out - inp)

        def getMin(node, inp, d):
            w = W[node.name]
            l = C
            p, m = w, getCost(inp, w, d) + l[w]
            for out in xrange(w + 1, M[node.name] + 1):
                t = getCost(inp, out, d) + l[out]
                if t < m:
                    p, m = out, t
            U[node.name][inp] = m
            I[node.name][inp] = p

        def getC(node, w):
            if len(node.get_children()) == 1:
                left_u = U[node.get_child_at(0).name]
                for j in xrange(w, M[node.name] + 1):
                    C[j] = left_u[j - w]

            elif len(node.get_children()) == 2:
                left_u = U[node.get_child_at(0).name]
                right_u = U[node.get_child_at(1).name]
                for j in xrange(w, M[node.name] + 1):
                    C[j] = left_u[j - w] + right_u[j - w]

        for node in image_tree.traverse("postorder"):
            mult_node = W[node.name]

            if node.is_leaf():
                U[node.name] = np.zeros(M[node.parent] - W[node.parent] + 1)
                I[node.name] = np.zeros(
                    M[node.parent] - W[node.parent] + 1, dtype=int)
                depth_diff = node.depth - node.up.depth
                for j in xrange(0, M[node.parent] - W[node.parent] + 1):
                    U[node.name][j] = getCost(j, mult_node - 1, depth_diff)
                    I[node.name][j] = mult_node - 1

            elif node.is_root():
                U[node.name] = np.zeros(1)
                I[node.name] = np.zeros(1, dtype=int)
                getC(node, mult_node)
                getMin(node, 0, 0)

            else:
                U[node.name] = np.zeros(M[node.parent] - W[node.parent] + 1)
                I[node.name] = np.zeros(
                    M[node.parent] - W[node.parent] + 1, dtype=int)
                getC(node, mult_node)
                depth_diff = node.depth - node.up.depth
                for j in xrange(0, M[node.parent] - W[node.parent] + 1):
                    getMin(node, j, depth_diff)

        ingene[image_tree.name] = I[image_tree.name][0]
        outgene[image_tree.name] = 0
        for node in reversed(image_tree.get_descendants("postorder")):
            outgene[node.name] = ingene[node.parent] - W[node.parent]
            depth_diff = node.depth - node.up.depth

            if (depth_diff * self.losscost) > self.dupcost:
                outgene[node.name] = 0

            ingene[node.name] = I[node.name][outgene[node.name]]

        # print( 'outgene ===> ', outgene)
        # print( 'ingene ===> ', ingene_
        # print( 'I == >', I_
        # print( 'W == >', W_
        # print( 'M == >', M)

        return self.get_solution(image_tree, W, rW, outgene, ingene)


class NotungSolver(Solver):

    def __init__(self, genetree, specietree, lcamap, dupcost=1, losscost=1):
        super(self.__class__, self).__init__(genetree, specietree, lcamap)
        # current version is not working if you want more that one solution
        # will let it that way for now, because we don't really need that in the simulation
        # comparision
        self.dupcost = dupcost
        self.losscost = losscost

    def reconstruct(self):
        for node in self.genetree.traverse("postorder"):

            if not node.is_leaf() and not node.is_binary():
                node_hash = TreeUtils.treeHash(node)
                W, rW = self._compute_mult(node)
                self.mincost = ddict(partial(np.zeros))
                self.cost = ddict(partial(np.zeros))
                self.out = ddict(int)
                self.losses = {}
                self.dups = {}
                self.changed = {}
                self.currDup = ddict(int)
                self.currLoss = ddict(int)
                self.currSpec = ddict(int)
                for resolved_node in self.Reconstruct(self.lcamap[node], W, rW, len(node.get_children())):
                    node.replace_by(resolved_node)

        return self.genetree.write(format=9)

    def Reconstruct(self, specT, W, rW, n, k=1):
        """Reconstruct (see notung paper)"""
        maxW = n  # max(W.values())
        self.Ascend(specT, maxW, W)
        reset = True
        i = 1
        while True:
            self.changed[specT.name] = self.Descend(specT, reset, 1, W, maxW)
            if not (self.changed[specT.name]):
                break
            if(i > k):
                break
            geneT, final_event = self.Construct(specT, rW)
            i += 1
            yield geneT
            # print( geneT.get_ascii(show_internal=True, attributes=['name']))
            reset = False

    def Ascend(self, node, maxW, W):
        mv = W[node.name]
        self.mincost[node.name] = np.zeros(maxW, dtype=int)
        if node.is_leaf():
            for i in xrange(1, maxW + 1):
                self.mincost[node.name][i - 1] = self.dupcost * \
                    max(mv - i, 0) + self.losscost * max(i - mv, 0)
        else:
            self.cost[node.name] = np.zeros((maxW, maxW), dtype=int)

            lchild = node.get_child_at(0)
            rchild = node.get_child_at(1)
            self.Ascend(lchild, maxW, W)
            self.Ascend(rchild, maxW, W)
            for i in xrange(1, maxW + 1):
                for j in xrange(1, maxW + 1):
                    self.cost[node.name][i - 1, j - 1] = self.dupcost * max(j - i + mv, 0) + self.losscost * max(
                        i - mv - j, 0) + self.mincost[lchild.name][j - 1] + self.mincost[rchild.name][j - 1]
            for i in xrange(1, maxW + 1):
                self.mincost[node.name][
                    i - 1] = np.min(self.cost[node.name][i - 1, :])

    def Descend(self, node, reset, i, W, maxW):
        mv = W[node.name]
        if node.is_leaf():
            self.out[node.name] = 0
            self.losses[node.name] = max(i - mv, 0)
            self.dups[node.name] = max(mv - i, 0)
            return reset

        lchild = node.get_child_at(0)
        rchild = node.get_child_at(1)
        if not reset:
            # Check each child - if either has another solution
            # return True
            self.changed[lchild.name] = self.Descend(
                lchild, False, self.out[node.name], W, maxW)
            if self.changed[lchild.name]:
                return True
            self.changed[rchild.name] = self.Descend(
                rchild, False, self.out[node.name], W, maxW)
            if self.changed[rchild.name]:
                return True

            # In the last case, neither child has another solution, so we choose the next
            # optimal value of out[node.name]
            self.out[node.name] += 1
            while not (self.out[node.name] > maxW or self.cost[node.name][i - 1, self.out[node.name] - 1] == self.mincost[node.name][i - 1]):
                self.out[node.name] += 1

            if self.out[node.name] > maxW:
                return False
        else:
            self.out[node.name] = 1
            while self.cost[node.name][i - 1, self.out[node.name] - 1] != self.mincost[node.name][i - 1]:
                self.out[node.name] += 1

        # We need the values for a new out[node.name]
        self.changed[lchild.name] = self.Descend(
            lchild, True, self.out[node.name], W, maxW)
        self.changed[rchild.name] = self.Descend(
            rchild, True, self.out[node.name], W, maxW)
        self.losses[node.name] = max(i - self.out[node.name] - mv, 0)
        self.dups[node.name] = max(self.out[node.name] - i + mv, 0)
        return True

    def Construct(self, node, rW):
        g = TreeClass()
        g.name = node.name
        event = None
        if self.currDup[node.name] < self.dups[node.name]:
            self.currDup[node.name] += 1
            event = "dup"
            g_lchild, el = self.Construct(node, rW)
            g_rchild, er = self.Construct(node, rW)
            g.add_features(type=TreeClass.AD)
            # this is a little strange, it's like we have a duplication
            # then one of the child go extinct ==> this just look like a
            # speciation in the tree imo
            if el == 'lost' and er == 'lost':
                raise ValueError("This case should not happen")
                self.currDup[node.name] -= 1

            if el == 'lost':
                g = g_rchild
            elif er == 'lost':
                g = g_lchild
            else:
                g.add_child(g_lchild)
                g.add_child(g_rchild)

        elif self.currLoss[node.name] < self.losses[node.name]:
            self.currLoss[node.name] += 1
            event = "lost"
            g.add_features(type=TreeClass.LOST)

        elif self.currSpec[node.name] < self.out[node.name]:
            self.currSpec[node.name] += 1
            event = "spec"
            g_lchild, el = self.Construct(node.get_child_at(0), rW)
            g_rchild, er = self.Construct(node.get_child_at(1), rW)

            if er == 'lost' and el == 'lost':
                event = 'lost'
                g.add_child(g_lchild)
                g.add_child(g_rchild)
                g.add_features(type=TreeClass.LOST)

            elif el == 'lost':
                g = g_rchild
            elif er == 'lost':
                g = g_lchild

            else:
                g.add_child(g_lchild)
                g.add_child(g_rchild)
                g.add_features(type=TreeClass.SPEC)

        else:
            event = 'None'
            g = rW[node.name].pop()
            g.name = node.name

        return g, event


class Dynamiq2(Solver):

    def __init__(self, genetree, specietree, lcamap, dupcost=1, losscost=1):
        super(self.__class__, self).__init__(genetree, specietree, lcamap)
        self.dupcost = dupcost
        self.losscost = losscost

    def reconstruct(self):
        images_trees = TreeUtils.getImageTreeNode(
            self.genetree, self.specietree, self.lcamap)

        for node in self.genetree.traverse("postorder"):
            if (node.is_root() or node.is_internal()) and not node.is_binary():
                for n in images_trees[node].traverse():
                    if(n.up is not None):
                        n.add_features(parent=n.up.name)
                W, rW = self._compute_mult(node)
                solution = self.compute_table(images_trees[node], node, W, rW)
                node.replace_by(solution)
        return self.genetree.write(format=9)

    def compute_table(self, image_tree, node, W, rW):
        """Compute cost table"""

        child_len = len(node.get_children())
        costTable = ddict(partial(np.zeros))
        ingene = ddict(int)
        outgene = ddict(int)
        # print( image_tree)

        def getCost(w_u, k):
            return self.dupcost if w_u > k else self.losscost

        def getWbar(c_u):
            return min(self.dupcost + self.losscost, c_u * self.losscost)

        for node in image_tree.traverse("postorder"):
            costTable[node.name] = np.zeros(child_len)
            ingene[node.name] = +np.inf
            c_u = 0
            if(node.up is not None):
                c_u = node.depth - node.up.depth - 1
            w_bar = getWbar(c_u)
            w_u = W[node.name]
            for k in xrange(1, child_len + 1):
                w = getCost(w_u, k)

                if len(node.get_children()) == 0:
                    costTable[node.name][k - 1] = w_bar * \
                        min(k, w_u) + w * abs(k - w_u)
                    ingene[node.name] = w_u - 1

                elif len(node.get_children()) == 1:
                    child = node.get_children()[0]
                    val_on_kp = []
                    for k_p in xrange(1, child_len + 1):
                        dt = w_bar * min(k, w_u) + w * abs(k - w_u - k_p)
                        val_on_kp.append(
                            dt + k_p * self.losscost + costTable[child.name][k_p - 1])
                    min_kp, costTable[node.name][
                        k - 1] = min(enumerate(val_on_kp), key=operator.itemgetter(1))
                    ingene[node.name] = min(min_kp + w_u, ingene[node.name])

                elif len(node.get_children()) == 2:
                    child1 = node.get_children()[0]
                    child2 = node.get_children()[1]
                    val_on_kp = []
                    for k_p in xrange(1, child_len + 1):
                        dt = w_bar * min(k, w_u) + w * abs(k - w_u - k_p)
                        val_on_kp.append(
                            dt + costTable[child1.name][k_p - 1] + costTable[child2.name][k_p - 1])
                    min_kp, costTable[node.name][
                        k - 1] = min(enumerate(val_on_kp), key=operator.itemgetter(1))
                    ingene[node.name] = min(min_kp + w_u, ingene[node.name])

        outgene[image_tree.name] = 0

        for node in reversed(image_tree.get_descendants("postorder")):
            outgene[node.name] = ingene[node.parent] - W[node.parent]
            depth_diff = node.depth - node.up.depth

            if (depth_diff * self.losscost) > self.dupcost:
                outgene[node.name] = 0

        return self.get_solution(image_tree, W, rW, outgene, ingene)
