"""
Author: Emmanuel Noutahi
Date: 05/2014
MultiPolysolver is a python module for polytomy solving
"""

import numpy
import random
from pprint import pprint
import copy

from ..TreeLib import *
"""
Gene matrix are represented by a numpy array
"""

DUP = 'd'
LOST = 'l'
SPEC = 's'
PARTIAL_RESOLUTION_ITERATOR = 1
numpy.set_printoptions(threshold='nan', precision=10)


@memorize
def polySolver(genetree, specietree, gene_matrix, node_order, limit=-1, cluster_method='upgma', verbose=False, mode="solve"):
    """This assume we that we are using the correct specietree for this genetree
    the specie tree root is the latest common ancestor of all the specie in genetree"""
    count = TreeUtils.getSpecieCount(
        genetree)  # number of specie in the genetree
    max_y = max(count.values()) + 1
    # assigning a correspondance between each row and a node
    # the assignment is done in level order
    polytomy_specie_set, row_node_corr = findMaxX(genetree, specietree)
    max_x = len(polytomy_specie_set)
    # cost cost_table to fill
    cost_table = numpy.zeros((max_x, max_y), dtype=float)
    # table to save the possible path
    path_table = numpy.ndarray((max_x, max_y), dtype='object')

    # fill the cost_table and the path_table

    for n in xrange(0, max_x):
        node = row_node_corr[n]
        zeropos = count[node.name] - 1
        # We have zeropos when the number of node from a specie is the same as the column number
        # Fill the table, using the next/previous case cost
        # The node is a leaf, just fill with dupcost and losscost
        if(node.is_leaf()):
            # find the column with a cost of zero (zeropos) and fill the table
            # according to this position
            # by default, all the position in the table are 0
            i = zeropos - 1
            while(i >= 0):
                cost_table[n, i] = cost_table[n, i + 1] + params.getdup(node)
                path_table[n, i] = DUP
                i -= 1
            i = zeropos + 1
            while(i < max_y):
                cost_table[n, i] = cost_table[n, i - 1] + params.getloss(node)
                path_table[n, i] = LOST
                i += 1
            # We should take into account the special case here
        # Here we have an internal node (not a leaf in the genetree)
        else:
            l_child_id = [key for key, value in row_node_corr.iteritems() if(
                value == node.get_child_at(0))][0]
            r_child_id = [key for key, value in row_node_corr.iteritems() if(
                value == node.get_child_at(1))][0]
            # Fill the table using only the speciation cost(sum of the
            # children's cost of this node)
            for k in xrange(0, max_y):
                if(k > zeropos):
                    cost_table[n, k] = cost_table[
                        l_child_id, k - zeropos - 1] + cost_table[r_child_id, k - zeropos - 1]
                    path_table[n, k] = SPEC
                else:
                    cost_table[n, k] = numpy.inf

            # Find all the min score position and try to minimize the score of its
            # neighborhood by lost/dup cost
            minpos = numpy.where(cost_table[n, :] == cost_table[n, :].min())
            for pos in minpos[0]:
                i = pos - 1
                while(i >= 0):
                    if(cost_table[n, i] == cost_table[n, i + 1] + params.getdup(node)):
                        if DUP not in path_table[n, i]:
                            path_table[n, i] += DUP
                    elif (cost_table[n, i] > cost_table[n, i + 1] + params.getdup(node)):
                        cost_table[n, i] = min(
                            cost_table[n, i + 1] + params.getdup(node), cost_table[n, i])
                        path_table[n, i] = DUP
                    i -= 1

                i = pos + 1
                while(i < max_y):
                    if(cost_table[n, i] == cost_table[n, i - 1] + params.getloss(node)):
                        if LOST not in path_table[n, i]:
                            path_table[n, i] += LOST
                    elif (cost_table[n, i] > cost_table[n, i - 1] + params.getloss(node)):
                        cost_table[n, i] = min(
                            cost_table[n, i - 1] + params.getloss(node), cost_table[n, i])
                        path_table[n, i] = LOST
                    i += 1

    # find the shape of the cost_table
    xsize, ysize = cost_table.shape
    if(mode is not "solve"):
        return cost_table, row_node_corr
    else:
        paths = findPathFromTable(
            path_table, row_node_corr, count, xsize - 1, 0)
        solution = []

        if(verbose):
            print("Matrix M: \n")
            print(cost_table)
            print()
            print("Path table for tree construction: \n")
            print(path_table)
            print()
            print("Correspondance: \n")
            pprint(row_node_corr)
            print()
            print("Gene Tree:\n")
            print(genetree.get_ascii(attributes=['species', 'name']))
            print()
            print("Specie Tree:\n")
            print(specietree.get_ascii())
            print("\nNumber of Tree found : ", len(paths), "\n")
            print("List of possible path: ")
            for path in paths:
                print(path)
            print()

        i = 1
        for path in paths:
            if(limit > 0 and i > limit):
                break
            solution.append(constructFromPath(path, genetree, specietree, numpy.copy(gene_matrix), node_order[
                            :], verbose=verbose, method=cluster_method, cost=cost_table[xsize - 1, 0]))
            i += 1

        return solution


def findSpeciationPathFromTable(path_table, row_node_corr, count, xpos, ypos):
    """DEBUG, choose the path that privilegie speciation only"""

    chemin = []
    if(row_node_corr[xpos].is_leaf() and (ypos < 0 or path_table[xpos, ypos] is None)):
        case = row_node_corr[xpos].name + ':%i' % (ypos + 1)
        chemin.append(case)
    else:
        # each case can have multiple path
        if 's' in path_table[xpos, ypos]:
            # this is bad for perfomance
            spec_pos_1 = [x for x in row_node_corr.keys() if row_node_corr[
                xpos].get_child_at(0) == row_node_corr[x]][0]
            spec_pos_2 = [x for x in row_node_corr.keys() if row_node_corr[
                xpos].get_child_at(1) == row_node_corr[x]][0]
            snode = row_node_corr[xpos].name
            nb_node = count[snode]
            spec_1 = findSpeciationPathFromTable(
                path_table, row_node_corr, count, spec_pos_1, ypos - nb_node)
            spec_2 = findSpeciationPathFromTable(
                path_table, row_node_corr, count, spec_pos_2, ypos - nb_node)
            # add all possible path from the children
            chemin.extend([",".join([row_node_corr[xpos].name + ':%i' %
                                     (ypos + 1), path1, path2]) for path1 in spec_1 for path2 in spec_2])

        else:
            c = path_table[xpos, ypos][0]
            if c == 'd':
                dup = findSpeciationPathFromTable(
                    path_table, row_node_corr, count, xpos, ypos + 1)
                # add possible path of the case that lead to this duplication
                chemin.extend(
                    [",".join([row_node_corr[xpos].name + ':%i' % (ypos + 1), path1]) for path1 in dup])

            # instead we found a lost
            elif c == 'l':
                lost = findSpeciationPathFromTable(
                    path_table, row_node_corr, count, xpos, ypos - 1)
                # add possible path of the case that lead to this lost
                chemin.extend(
                    [",".join([row_node_corr[xpos].name + ':%i' % (ypos + 1), path1]) for path1 in lost])

    # return a list of path, the path are specially formated string
    return chemin


def findPathFromTable(path_table, row_node_corr, count, xpos, ypos):
    """ Find all the possible path from the lower left case to the leaves"""

    chemin = []
    # Case 1: current position correspond to a leaf
    if(row_node_corr[xpos].is_leaf() and (ypos < 0 or path_table[xpos, ypos] is None)):
        case = row_node_corr[xpos].name + ':%i' % (ypos + 1)
        chemin.append(case)

    # Case 2 : this a internal node
    else:
        # each case can have multiple path
        for c in path_table[xpos, ypos]:

            # we found a speciation
            if c == 's':
                # this is bad for perfomance
                spec_pos_1 = [x for x in row_node_corr.keys() if row_node_corr[
                    xpos].get_child_at(0) == row_node_corr[x]][0]
                spec_pos_2 = [x for x in row_node_corr.keys() if row_node_corr[
                    xpos].get_child_at(1) == row_node_corr[x]][0]
                snode = row_node_corr[xpos].name
                nb_node = count[snode]
                spec_1 = findPathFromTable(
                    path_table, row_node_corr, count, spec_pos_1, ypos - nb_node)
                spec_2 = findPathFromTable(
                    path_table, row_node_corr, count, spec_pos_2, ypos - nb_node)
                # add all possible path from the children
                for path1 in spec_1:
                    for path2 in spec_2:
                        chemin.append(
                            ",".join([row_node_corr[xpos].name + ':%i' % (ypos + 1), path1, path2]))

            # we found a duplication
            elif c == 'd':
                dup = findPathFromTable(
                    path_table, row_node_corr, count, xpos, ypos + 1)
                # add possible path of the case that lead to this duplication
                for path1 in dup:
                    chemin.append(
                        ",".join([row_node_corr[xpos].name + ':%i' % (ypos + 1), path1]))

            # instead we found a lost
            elif c == 'l':
                lost = findSpeciationPathFromTable(
                    path_table, row_node_corr, count, xpos, ypos - 1)
                # add possible path of the case that lead to this lost
                for path1 in lost:
                    chemin.append(
                        ",".join([row_node_corr[xpos].name + ':%i' % (ypos + 1), path1]))

    # return a list of path, the path are specially formated string
    return chemin


def constructFromPath(chemin, genetree, specietree, gene_matrix, node_order, verbose=False, method='upgma', cost=0):
    """Construct tree from a path using the clustering method"""
    # get the node order in the path
    node_list = list(reversed(chemin.split(',')))
    # find the node to show in the tree construction
    node_in_tree = genetree.get_children()
    leaf_list = [x.name for x in specietree.get_leaves()]
    gene_tree_desc_species = genetree.get_descendant_species()
    species_list = [x.name for x in specietree.traverse()]
    gene_root_species = specietree.get_common_ancestor(
        [x for x in gene_tree_desc_species if(x in species_list)])
    # leaf_list = genetree.get_children_species()
    # keep this in case a lost node is in the path
    lost_nodes = []
    # total number of node
    tot_node = len(node_list)
    deja_vu = []
    # traverse the list of node in the path and construct the tree
    for indice in xrange(tot_node):  # 0 - len(node_list)-1

        # find the current node and its position
        [node, pos] = node_list[indice].split(":")
        # find the next node and its position
        [n_node, n_pos] = [None, -1]
        if(indice + 1 != tot_node):
            [n_node, n_pos] = node_list[indice + 1].split(":")
        pos = int(pos)
        n_pos = int(n_pos)

        # list of node which specie is the same as the current node
        node_structs = [n for n in node_in_tree if n.species == node]
        # index to keep in node_order for the join
        ind_to_keep = []

        # simple case, the node is a leaf
        if(node in leaf_list):
            # we have a duplication here,
            if(n_node == node and n_pos < pos):
                # ind to keep is empty on purpose, we have to find the good index
                # and remove the row we don't want the clustering algorithm to
                # work with
                cluster = ClusterUtils.treeCluster(getMatrix(
                    node_structs, gene_matrix, node_order, ind_to_keep), node_structs, 1, method=method)
                # find the resulting duplication tree
                dup_tree = cluster[0]
                # set the name of the duplication, i'll try noName next time
                dup_tree.name = "-".join(
                    [dup_tree.get_child_at(0).name, dup_tree.get_child_at(1).name])
                dup_tree.add_features(type="DUP")
                # find the merged index and remove them from the distance
                # matrice
                merged_index = map(lambda x: ind_to_keep[x], cluster[2])
                merged_index.sort(reverse=True)
                dup_tree.add_features(species=node)
                # Update the gene matrix and the node order
                gene_matrix = ClusterUtils.condense_matrix(
                    gene_matrix, merged_index, method=method)
                node_order[merged_index[1]] = dup_tree.name
                node_order.pop(merged_index[0])

                for n in dup_tree.get_children():
                    node_in_tree.remove(n)
                node_in_tree.append(dup_tree)

            # the leaf is a lost leaf, add it to the list of lost_nodes
            # which will be checked when we want to construct a internal node
            # with a lost child
            elif(n_node == node and n_pos > pos):
                lost = TreeClass()
                lost.name = node + '_%i' % (n_pos)
                lost.species = node
                lost.add_features(type='LOST')
                lost_nodes.append(lost)

        else:  # internal node case
            if(n_node == node and pos < n_pos):
                "Lost internal node"
                lost = TreeClass()
                lost.name = node + '_%i' % (n_pos)
                lost.species = node
                lost.add_features(type='LOST')
                lost_nodes.append(lost)

            elif(pos >= 1 and node not in deja_vu):
                i = pos
                deja_vu.append(node)
                # Case where we have this node n times, n>=1. We should then
                # construct the node all the n times
                i = len([x for x in node_in_tree if x.name == node])
                while(i < pos):
                    # Speciation case:
                    s_node = specietree & (node)
                    right_child = s_node.get_child_at(0)
                    left_child = s_node.get_child_at(1)

                    # Case where the child we're looking for is not even in the
                    # geneTree
                    r_child_not_there, l_child_not_there = (False, False)

                    if(right_child.name not in gene_root_species.get_descendant_name()):
                        r_child_not_there = True

                    elif(left_child.name not in gene_root_species.get_descendant_name()):
                        l_child_not_there = True

                    s_node.add_features(species=node)
                    child_name = s_node.get_children_name()
                    # Now we should find the index in the matrix of the best
                    # node to join
                    ind_to_keep = findSpeciationBestJoin(
                        gene_matrix, node_order, s_node, node_in_tree, method=method, verbose=verbose)
                    # the two nodes to join are found
                    if(ind_to_keep and len(ind_to_keep) == 2):
                        # node_structs= [ node_s for node_s in node_in_tree if (node_s.name in map(lambda x:node_order[x],ind_to_keep ) and ((node_s.is_leaf() and node_s.species == genetree.search_nodes(name=node_s.name)[0].species) or not(node_s.is_leaf())))]
                        # carefully choose the true node, not the added(in case
                        # of lost)
                        node_structs = [node_s for node_s in node_in_tree if (node_s.name in map(lambda x:node_order[
                                                                              x], ind_to_keep) and (not node_s.has_feature('lostnode') or not node_s.lostnode == 1))]
                        cluster = ClusterUtils.treeCluster(getMatrix(
                            node_structs, gene_matrix, node_order, ind_to_keep, got_ind=True), node_structs, method=method)
                        spec_tree = cluster[0]
                        spec_tree.name = "-".join(
                            [spec_tree.get_child_at(0).name, spec_tree.get_child_at(1).name])
                        spec_tree.add_features(type="SPEC")
                        merged_index = map(
                            lambda x: ind_to_keep[x], cluster[2])
                        merged_index.sort(reverse=True)

                        spec_tree.add_features(species=node)
                        gene_matrix = ClusterUtils.condense_matrix(
                            gene_matrix, merged_index, method=method)
                        node_order[merged_index[1]] = spec_tree.name
                        node_order.pop(merged_index[0])
                        # Remove child from path_table and add parent
                        possible_path_name = [x.name for x in node_in_tree]
                        for n in spec_tree.get_children_name():
                            while(n in possible_path_name):
                                node_in_tree.pop(possible_path_name.index(n))
                                possible_path_name.pop(
                                    possible_path_name.index(n))
                        node_in_tree.append(spec_tree)

                        # too much repetition, should make a function to do
                        # this part and the duplication part

                    # we certainly have a case of a lost node here
                    else:
                        # retrieve list of child and list of lost child
                        r_child_list = [
                            x for x in node_in_tree if x.species == right_child.name]
                        r_child_lst_list = [
                            x for x in lost_nodes if x.species == right_child.name]

                        l_child_list = [
                            x for x in node_in_tree if x.species == left_child.name]
                        l_child_lst_list = [
                            x for x in lost_nodes if x.species == left_child.name]

                        # Control 2, list of children not present in genetree
                        # check which node is lost and make the correct choice
                        if(not r_child_list and r_child_lst_list):
                            "Cas d'un noeud avec une perte droite"
                            for n in l_child_list:
                                n_copy = n.copy()
                                n_copy.add_features(species=s_node.species)
                                n_copy.add_features(lostnode=1)
                                already_in = [nd for nd in node_in_tree if(
                                    nd.name == n_copy.name and nd.species == n_copy.species and nd.has_feature('lostnode'))]
                                if(not already_in):
                                    node_in_tree.append(n_copy)

                        elif(not l_child_list and l_child_lst_list):
                            "Cas d'un noeud avec une perte gauche"
                            for n in r_child_list:
                                n_copy = n.copy()
                                n_copy.add_features(species=s_node.species)
                                n_copy.add_features(lostnode=1)
                                already_in = [nd for nd in node_in_tree if(
                                    nd.name == n_copy.name and nd.species == n_copy.species and nd.has_feature('lostnode'))]
                                if(not already_in):
                                    node_in_tree.append(n_copy)

                        if(r_child_not_there or l_child_not_there):
                            "Cas de l'absence de l'espece dans le geneTree"
                            try:
                                n = right_child if r_child_not_there else left_child
                                match_node = [nd for nd in node_in_tree if (
                                    nd.species == s_node.name and nd.has_feature('lostnode') and nd.lostnode == 1)]
                                for nd in match_node:
                                    node_in_tree.remove(nd)
                            except:
                                pass
                    i += 1

            # we have a duplication here, construct the duplicated node
            if(n_node == node and pos > n_pos):
                node_structs = [x for x in node_in_tree if x.species == node]
                node_struct_name = [x.name for x in node_structs]
                node_structs.extend(
                    [x for x in node_in_tree if(x.name in node_struct_name and x not in node_structs)])
                node_structs = [node_s for node_s in node_structs if (
                    not node_s.has_feature('lostnode') or not node_s.lostnode == 1)]
                ind_to_keep = []
                cluster = ClusterUtils.treeCluster(getMatrix(
                    node_structs, gene_matrix, node_order, ind_to_keep), node_structs, 1, method=method)
                dup_tree = cluster[0]
                dup_tree.name = "-".join(
                    [dup_tree.get_child_at(0).name, dup_tree.get_child_at(1).name])
                dup_tree.add_features(type="DUP")

                merged_index = map(lambda x: ind_to_keep[x], cluster[2])
                merged_index.sort(reverse=True)
                dup_tree.add_features(species=node)
                gene_matrix = ClusterUtils.condense_matrix(
                    gene_matrix, merged_index, method=method)
                node_order[merged_index[1]] = dup_tree.name
                del node_order[merged_index[0]]

                possible_path_name = [x.name for x in node_in_tree]
                for n in dup_tree.get_children_name():
                    while(n in possible_path_name):
                        node_in_tree.pop(possible_path_name.index(n))
                        possible_path_name.pop(possible_path_name.index(n))
                node_in_tree.append(dup_tree)

    # the tree is constructed, we should only have one tree in the
    # node_in_tree list
    if(node_in_tree and len(node_in_tree) == 1):
        if(verbose):
            print("Total number of node = ", len(node_in_tree[
                  0]), " root specie = ", node_in_tree[0].species, "\n\n")
        node_in_tree[0].add_features(cost=cost)
        return node_in_tree[0]
    else:
        # Display error or raise exception ??
        if verbose:
            print('Specie')
            print(specietree.get_ascii(show_internal=True))
            print("\n Node in Genetree\n\n")
            print(chemin)
            for node in node_in_tree:
                print(node.get_ascii(
                    show_internal=True, attributes=['species']))
        raise ValueError(
            "Cannot construct your tree, %i node still not used !\n" % len(node_in_tree))


def findSpeciationBestJoin(matrice, node_order, parent_node, node_in_tree, method='upgma', verbose=False):
    """ findSpeciationBestJoin find the best node to joint in case of speciation"""

    child_0 = parent_node.get_child_at(0).name  # find the left child specie
    child_1 = parent_node.get_child_at(1).name

    # list of node with the same specie as child_0/child_1
    child_1_list = []
    child_0_list = []

    for x in node_in_tree:
        if x.species == child_1:
            if x.name in node_order:
                child_1_list.append(node_order.index(x.name))
            else:
                child_1_list.append(-1)

        elif x.species == child_0:

            if(x.name in node_order):
                child_0_list.append(node_order.index(x.name))
            else:
                child_0_list.append(-1)

    # find the best node to join (minimal cost)
    min_val = numpy.inf
    join_index = []

    if(verbose):
        print("Using %s as clustering method" % (method))

    if(method == 'upgma'):
        for x_0 in child_0_list:
            for x_1 in child_1_list:
                if(x_0 >= 0 and x_1 >= 0 and matrice[x_0, x_1] < min_val):
                    min_val = matrice[x_0, x_1]
                    join_index = [x_0, x_1]

    elif(method == "nj"):
        mat_size = matrice.shape[0]
        for x_0 in child_0_list:
            for x_1 in child_1_list:
                if(x_0 >= 0 and x_1 >= 0):
                    val = ClusterUtils.calculate_Q_ij(
                        matrice, (x_0, x_1), mat_size)
                    if val < min_val:
                        min_val = val
                        join_index = [x_0, x_1]

    # this is the case we have rand as method
    else:
        list_0 = [x for x in child_0_list if x > -1]
        list_1 = [x for x in child_1_list if x > -1]
        join_index = [list_0[0], list_1[0]] if (
            len(list_0) > 0 and len(list_1) > 0) else []

    return join_index


def getMatrix(node_struct, gene_matrix, node_order, ind_to_keep, got_ind=False):
    """Extract the correct position of the gene_matrix for the next cluster join"""

    if(not got_ind):
        for node in node_struct:
            ind_to_keep.append(getIndex(node_order, node))
    matrix = gene_matrix  # same reference here
    return numpy.take(numpy.take(matrix, ind_to_keep, axis=0), ind_to_keep, axis=1)


def getIndex(node_order, node):
    """Get the index of a node in a node list"""
    return node_order.index(node.name)


def polytomyPreprocess(polytomy, specietree, gene_matrix, node_order, method='upgma'):
    """Preprocessing of a polytomy """
    lcamap = TreeUtils.lcaMapping(polytomy, specietree, multspeciename=False)
    for node in polytomy.traverse("postorder"):
        if(node.is_binary()):
            ind_order = []
            if(	node.name == TreeClass.DEFAULT_NAME or node.name == "NoName"):
                global PARTIAL_RESOLUTION_ITERATOR
                node.name = "%s_I_%i" % (
                    node.species, PARTIAL_RESOLUTION_ITERATOR)
                while node.name in node_order:
                    PARTIAL_RESOLUTION_ITERATOR += 1
                    node.name = "%s_I_%i" % (
                        node.species, PARTIAL_RESOLUTION_ITERATOR)
                PARTIAL_RESOLUTION_ITERATOR += 1

            for child_node in node.get_children():
                ind_order.append(getIndex(node_order, child_node))

            ind_order.sort(reverse=True)
            gene_matrix = ClusterUtils.condense_matrix(
                gene_matrix, ind_order, method=method)
            node_order[ind_order[1]] = node.name
            del node_order[ind_order[0]]
    return gene_matrix, node_order


def _compute_mult(polytomy, lcamap):
    W = ddict(int)
    for g in polytomy.get_children():
        s = s.lcamap[g]
        W[s.name] += 1
    return W


def findMaxX(polytomy, specietree):
    """Find Number of Specie and the specie list in order to create and fill the dup/cost matrix"""
    if not polytomy.has_feature('species'):
        lcamap = TreeUtils.lcaMapping(
            polytomy, specietree, multspeciename=False)

    polytomy_specie_ancestor = (specietree & polytomy.species)
    polytomy_name_set = set(polytomy.get_children_species())

    for leaf in polytomy_specie_ancestor.traverse("postorder"):
        parent = leaf.up
        if(not leaf.is_leaf()):
            if(len(set(leaf.get_descendant_name()).intersection(polytomy_name_set)) == 0):
                parent.remove_child(leaf)

            else:
                polytomy_name_set.add(leaf.name)

        else:
            if(parent is not None) and (len(set(parent.get_descendant_name()).intersection(polytomy_name_set)) == 0):
                parent.remove_child(leaf)
            else:
                polytomy_name_set.add(leaf.name)

    row_node_corr = {}
    n_row = len(polytomy_name_set) - 1

    for node in specietree.traverse("levelorder"):
        if(node.name in polytomy_name_set):
            row_node_corr[n_row] = node
            n_row -= 1

    return polytomy_name_set, row_node_corr


def solvePolytomy(genetree, specietree, gene_matrix, node_order, verbose=False, path_limit=-1, method='upgma', sol_limit=-1):

    # Start with only one polytomy

    nb_polytomy = 0
    polysolution = [genetree]
    while True:
        next_tree_solution = []  # next list of partially resolved polytomies
        for tree in polysolution:
            for polytomy in tree.iter_polytomies(strategy="postorder"):
                nb_polytomy += 1
                # copying the input for each step, necessary in order to not
                # modify by reference
                matrice = numpy.copy(gene_matrix)
                sptree = specietree.copy("newick")
                ptree = polytomy.copy()
                order = node_order[:]
                poly_parent = polytomy.up
                node_to_replace = polytomy
                matrice, order = polytomyPreprocess(
                    ptree, sptree, matrice, order, method=method)
                solution = polySolver(TreeUtils.treeHash(ptree, addinfos=str(
                    path_limit) + method), ptree, sptree, matrice, order, path_limit, cluster_method=method, verbose=verbose)
                # solution=polySolver(ptree,sptree, matrice, order,path_limit, cluster_method=method, verbose=verbose)
                if(poly_parent is None):
                    # Here we have the root. Complete solution are here
                    next_tree_solution.extend(solution)

                else:
                    # This is one of the internal polytomy.
                    for sol in solution:
                        poly_parent.replace_child(node_to_replace, sol)
                        node_to_replace = sol
                        next_tree_solution.append(tree.copy())

                break  # only solve one polytomy per iteration

        if not next_tree_solution:
            break
        polysolution = next_tree_solution
        if(sol_limit > 0 and sol_limit < len(polysolution)):
            path_limit = 1

    if(nb_polytomy < 1):
        raise ValueError("Polytomy not found in your gene tree")

    # deepcopy is not working, neither is newick-copy
    f_sol = polysolution[0:sol_limit] if sol_limit > 0 else polysolution
    return [t.copy("simplecopy") for t in f_sol]


def computePolytomyReconCost(genetree, specietree, verbose=False):
    """This is a copy pasta from the solvePolytomy function that return only the cost of a node
    """
    recon_cost = 0
    # genetree = origene.copy(method="simplecopy")
    lcamap = TreeUtils.lcaMapping(genetree, specietree, multspeciename=False)

    for node in genetree.iter_internal_node(strategy="postorder", enable_root=True):

        if (node.is_binary()):

            # if node is binary, we just compute the recon cost
            try:
                cost = TreeUtils.binaryRecScore(node, lcamap)
                recon_cost += cost[0]

            except Exception as e:
                print(e)

        elif(node.is_polytomy()):
            # here, the node is a polytomy, so we compute the table and
            # find the solution cost at [-1, 0]
            sptree = specietree.copy("simplecopy")
            mat_table, row_node = polySolver(TreeUtils.treeHash(
                node), node, sptree, None, [], 1, verbose=verbose, mode="none")
            if(verbose):
                print(node)
                pprint(mat_table)
                print(
                    "%s : -------------------------------------------------------\n" % (recon_cost))
            recon_cost += mat_table[-1, 0]

        else:
            raise Exception("Internal node with only one child in your tree")
    return recon_cost
