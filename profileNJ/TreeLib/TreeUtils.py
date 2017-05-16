# This file is part of profileNJ
#
# Date: 02/2014
# TreeUtils is a python class that offer function related to phylogeny
# tree, using TreeClass

__author__ = "Emmanuel Noutahi"

import ClusterUtils as clu
import hashlib
import os
import params
import re
import string
import sys
import urllib2
import numpy as np
import random
from TreeClass import TreeClass
from collections import defaultdict as ddict
from ete3 import Phyloxml, Tree
from ete3 import orthoxml
from ete3.parser.newick import NewickError
import itertools


# TreeUtils:

class MatrixRep():
    def __init__(self, genetree, speciestree, defval=0):
        self.gtree = genetree
        self.stree = speciestree
        # keeping tree as key (in case the name was not set for internal nodes)
        self.gmap = dict((gn, i) for i, gn in enumerate(genetree.traverse("postorder")))
        self.smap = dict((sn, i) for i, sn in enumerate(speciestree.traverse("postorder")))
        self.matrix = np.empty((len(self.gmap), len(self.smap)))
        self.matrix.fill(defval)
        self.shape = self.matrix.shape
        self.inlist = self.smap.keys() + self.gmap.keys()
    
    def __len__(self):
        """Return the len of the longuest axis"""
        return max(self.shape)
    
    def __contains__(self, item):
        return item in self.inlist
    
    def __iter__(self):
        for (g,i_g) in self.gmap.items():
            for (s,i_s) in self.smap.items():
                yield (g,s, self.matrix[i_g, i_s])
    
    def _reformat_slice(self, slice, map):
        return id if isinstance(id, int) else map.get(id, None)

    def _get_new_index(self, index, map):
        start =  index.start
        stop =  index.stop
        step = index.step
        start = self._reformat_slice(start, map)
        stop = self._reformat_slice(stop, map)
        step = step if isinstance(step, int) else None
        return slice(start, stop, step)
        
    def __getitem__(self, index):
        """Indexing with int, string or slice"""
        # the following will return a whole row
        if isinstance(index, Tree):
            return self.matrix[self.gmap[index]]
        elif isinstance(index, int):
            return self.matrix[index]
        # in th folowing, we are returning a slice
        elif isinstance(index, slice):
            index = self._get_new_index(index, self.gmap)
            return self.matrix[index] 
        # we are accepting two slices here, no more
        elif len(index) != 2:
            raise TypeError("Invalid index type.")
        # Handle double indexing
        row_index, col_index = index
        if isinstance(row_index, Tree):
            row_index = self.gmap.get(row_index, None)
        if isinstance(col_index, Tree):
            col_index = self.smap.get(col_index, None)        
        elif isinstance(row_index, slice):
            row_index =  self._get_new_index(row_index, self.gmap)
        elif isinstance(col_index, slice):
            col_index =  self._get_new_index(col_index, self.smap)
        # let numpy manage the exceptions
        return self.matrix[row_index, col_index]
    
    def __setitem__(self, index, val):
        """Indexing with int, string or slice"""
        # the following will return a whole row
        if isinstance(index, Tree):
            self.matrix[self.gmap[index]]  =  val
        elif isinstance(index, int):
            self.matrix[index] = val
        # in th folowing, we are returning a slice
        elif isinstance(index, slice):
            index = self._get_new_index(index, self.gmap)
            self.matrix[index] = val
        # we are accepting two slices here, no more
        elif len(index) != 2:
            raise TypeError("Invalid index type.")
        # Handle double indexing
        row_index, col_index = index
        if isinstance(row_index, Tree):
            row_index = self.gmap.get(row_index, None)
        if isinstance(col_index, Tree):
            col_index = self.smap.get(col_index, None)        
        elif isinstance(row_index, slice):
            row_index =  self._get_new_index(row_index, self.gmap)
        elif isinstance(col_index, slice):
            col_index =  self._get_new_index(col_index, self.smap)
        # let numpy manage the exceptions
        self.matrix[row_index, col_index] = val

    
def fetch_ensembl_genetree_by_id(treeID=None, aligned=0, sequence="none", output="nh", nh_format="full"):
    """Fetch genetree from ensembl tree ID
    :argument treeID: the ensembl tree ID, this is mandatory
    :argument aligned: boolean (0/1), used with sequence to retrieve aligned sequence
    :argument sequence: none / protein /cdna /gene, should we retrieve sequence also?, work only with phyloxml nh_format
    :argument output: nh / phyloxml, type of output we are looking for!
    :argument nh_format: full / display_label_composite / simple / species / species_short_name / ncbi_taxon / ncbi_name / njtree / phylip, The format of the nh output, only useful when the output is set to nh
    """
    if not treeID:
        raise valueError('Please provide a genetree id')
    else:
        # http = httplib2.Http(".cache")
        server = "http://rest.ensembl.org"
        ext = "/genetree/id/%s?sequence=%s;aligned=%i" % (
            treeID, sequence, aligned)
        if(output == "nh"):
            ext = ext + ";nh_format=%s" % nh_format
        output = "text/x-" + output
        request = urllib2.Request(
            server + ext, headers={"Content-Type": output})
        resp = urllib2.urlopen(request)
        content = resp.read()
        # resp, content = http.request(server+ext, method="GET", headers={"Content-Type":output})
        if not resp.getcode() == 200:
            print("Invalid response: ", resp.getcode())
            raise ValueError('Failled to process request!')

        if(output.lower() != "text/x-phyloxml"):
            return TreeClass(content)
        else:
            return getTreeFromPhyloxml(content)


def _calc_KTB_rate(starting_rate, duration, roeotroe):
    """
    Returns a simulated rate for the head node of a tree when:
        * the tail node has rate ``starting_rate``
        * the time duration of the edge is ``duration``
        * the rate of evolution of the rate of evolution is ``roeotroe`` (this is
            the parameter nu in Kishino, Thorne, and Bruno 2001)
    ``rng`` is a random number generator.
    The model used to generate the rate is the one described by Kishino, Thorne,
    and Bruno 2001.  The descendant rates or lognormally distributed.
    The mean rate returned will have an expectation of ``starting_rate``
    The variance of the normal distribution for the logarithm of the ending rate
        is the product of ``duration`` and ``roeotroe``
    THIS FUNCTION AND ITS DESCRIPTION ARE FROM THE DENDROPY PROJECT
    """
    if starting_rate <= 0.0:
        raise ValueError("starting_rate must be positive in the KTB model")
    rate_var = duration*roeotroe
    if rate_var > 0.0:
        # Kishino, Thorne and Bruno corrected the tendency for the rate to
        #   increase seen in teh TKP, 1998 model
        mu = np.log(starting_rate) - (rate_var/2.0)
        return random.lognormvariate(mu, np.sqrt(rate_var))
    return starting_rate

def _calc_KTB_rates_crop(node, **kwargs):
    """Returns a descendant rate and mean rate according to the Kishino, Thorne,
    Bruno model.  Assumes that the min_rate <= starting_rate <= max_rate if a max
    and min are provided.
    rate is kept within in the [min_rate, max_rate] range by cropping at these
    values and acting is if the cropping occurred at an appropriate point
    in the duration of the branch (based on a linear change in rate from the
    beginning of the random_variate drawn for the end of the branch).
    THIS FUNCTION AND ITS DESCRIPTION ARE FROM THE DENDROPY PROJECT
    """
    duration = node.dist
    if node.up != None and node.up.has_feature('rate'):
        starting_rate = node.up.rate
    else:
        starting_rate = kwargs.get('starting_rate', None)
    roeotroe = kwargs.get('roeotroe', None)
    min_rate = kwargs.get('min_rate', None)
    max_rate = kwargs.get('max_rate', None)
    
    if starting_rate is None or roeotroe is None:
        raise ValueError("starting_rate and roeotroe should not be None")
    if roeotroe*duration <= 0.0:
        if (min_rate and starting_rate < min_rate) or (max_rate and starting_rate > max_rate):
            raise ValueError("Parent rate is out of bounds, but no rate change is possible")
    
    r = _calc_KTB_rate(starting_rate, duration, roeotroe)
    mr = (starting_rate + r)/2.0
    if max_rate and r > max_rate:
        assert(starting_rate <= max_rate)
        p_changing =  (max_rate - starting_rate)/(r - starting_rate)
        mean_changing = (starting_rate + max_rate)/2.0
        mr = p_changing*mean_changing + (1.0 - p_changing)*max_rate
        r = max_rate
    elif min_rate and r < min_rate:
        assert(starting_rate >= min_rate)
        p_changing = (starting_rate - min_rate)/(starting_rate - r)
        mean_changing = (starting_rate + min_rate)/2.0
        mr = p_changing*mean_changing + (1.0 - p_changing)*min_rate
    node.add_features(rate=r)
    node.add_features(length=node.dist)
    return mr


def relax_molecular_clock(tree, rate_fn=_calc_KTB_rates_crop, **kwargs):
    """Relax molecular clock, using rate_fn as the function (be it autocorrelated)
    or uncorrelated. Default is Kishino Thorne Bruno"""
    for t in tree.traverse():
        t.dist = rate_fn(t, **kwargs)
    return tree

def make_clock_like(tree):
    """Transform a tree into a one that follow 
    the molecular clock (ultrametric)
    """

    height = 0.0
    cs = tree.get_children()
    if len(cs) == 0:
        return 0
    child_height = []
    for i, child in enumerate(cs):
        child_val = make_node_clock_like(child) + child.dist
        child_height.append(child_val)
        height += child_val

    height /= len(cs)
    for i, child in enumerate(cs):
        scale_subtree_branches(child, height/child_height[i])
    
    return height

def scale_subtree_branches(node, factor):
    """Scale the branches lenght to its parent of a node 
    by factor
    """
    old_dist = node.dist
    node.dist = old_dist*factor
    for child in node.get_children():
        scale_subtree_branches(child, factor)
    

def fetch_ensembl_genetree_by_member(memberID=None, species=None, id_type=None, output="nh", nh_format="full"):
    """Fetch genetree from a member ID
    :argument memberID: the ensembl gene ID member of the tree to fetch, this is mandatory! EX: ENSG00000157764
    :argument species: Registry name/aliases used to restrict searches by. Only required if a stable ID is not unique to a species (not the case with Ensembl databases) EX: human, homo_sapiens
    :argument id_type: Object type to restrict searches to. Used when a stable ID is not unique to a single class. EX: gene, transcript
    :argument output: nh / phyloxml, type of output we are looking for!
    :argument nh_format: full / display_label_composite / simple / species / species_short_name / ncbi_taxon / ncbi_name / njtree / phylip, The format of the nh output, only useful when the output is set to nh
    """
    if not memberID:
        raise valueError('Please provide a genetree id')
    else:
        http = httplib2.Http(".cache")
        server = "http://rest.ensembl.org"
        ext = "/genetree/member/id/%s?" % (memberID)
        if species:
            ext = ext + "species=" + species + ";"
        if id_type:
            ext = ext + "object_type=" + id_type + ";"
        if(output == "nh"):
            ext = ext + "nh_format=%s;" % nh_format
        output = "text/x-" + output
        resp, content = http.request(
            server + ext, method="GET", headers={"Content-Type": output})
        if not resp.status == 200:
            print("Invalid response: ", resp.status)
            raise ValueError('Failled to process request!')
        if(output.lower() != "text/x-phyloxml"):
            return TreeClass(content)
        else:
            return getTreeFromPhyloxml(content)


def lcaPreprocess(tree):
    """Make an euler tour of this tree"""
    # root = tree.get_tree_root()
    tree.add_features(depth=0)
    tree.del_feature('euler_visit')
    for node in tree.iter_descendants("levelorder"):
        node.del_feature('euler_visit')
        # can write this because parent are always visited before
        node.add_features(depth=node.up.depth + 1)
    node_visited = tree._euler_visit([])
    # print tree.get_ascii(show_internal= True, attributes=['name', 'euler_visit', 'depth'])
    # number of element in array
    n = len(node_visited)
    m = int(np.ceil(np.log2(n)))
    rmq_array = np.zeros((n, m), dtype=int)
    for i in xrange(n):
        rmq_array[i, 0] = i
    for j in xrange(1, m):
        i = 0
        while(i + 2**j < n - 1):
            if(node_visited[rmq_array[i, j - 1]].depth < node_visited[rmq_array[(i + 2**(j - 1)), j - 1]].depth):
                rmq_array[i, j] = rmq_array[i, j - 1]
            else:
                rmq_array[i, j] = rmq_array[i + 2**(j - 1), j - 1]
            i += 1

    node_map = ddict()
    name2ind = ddict()
    for i in xrange(n):
        cur_node = node_visited[i]
        # bad practice
        try:
            node_map[cur_node] = min(node_map[cur_node], i)
        except:
            node_map[cur_node] = i
        name2ind[cur_node.name] = node_map[cur_node]

    tree.add_features(lcaprocess=True)
    tree.add_features(rmqmat=rmq_array)
    tree.add_features(ind2node=node_visited)
    tree.add_features(node2ind=node_map)
    tree.add_features(name2ind=name2ind)


def getLca(sptree, species):
    """This should be a faster lcamapping
    species should be a list of node"""
    if not sptree.has_feature('lcaprocess', True):
        lcaPreprocess(sptree)
    A = sptree.ind2node
    M = sptree.rmqmat
    # using the biggest interval should return the lca of all species
    if len(species) > 1:
        s_index = sorted([sptree.node2ind[spec] for spec in species])
    # print "s_index vaut :", s_index , " et taille est : ", len(sptree.node2ind), " et rmq est : ", sptree.rmqmat.shape
    # print sptree.ind2node
        i = s_index[0]
        j = s_index[-1]

    else:
        if isinstance(species[0], str):
            # in this case, we have a leaf
            i = sptree.name2ind[species[0]]
        else:
            # this is an instance of TreeClass
            i = sptree.node2ind[species[0]]
        j = i
    k = int(np.log2(j - i + 1))
    if (A[M[i, k]].depth <= A[M[j - 2**(k) + 1, k]].depth):
        return A[M[i, k]]
    else:
        return A[M[j - 2**(k) + 1, k]]


def lcaMapping(genetree, specietree, multspeciename=True):
    """LCA mapping between a genetree and a specietree
    :argument genetree: your genetree, All leave in the genetree should already have feature 'specie' (set_specie was called)
    :argument specietree: your specietree
    :argument multspeciename: A flag to use in order to accept multi specie name at genetree internal node.
    """

    smap = {}  # a dict that map specie name to specie node in specietree
    mapping = {}
    if not specietree.has_feature('lcaprocess', True):
        lcaPreprocess(specietree)

    for node in genetree.traverse(strategy="postorder"):

        if node.is_leaf():
            mapping[node] = getLca(specietree, [node.species])
        else:
            # ML ADDED THIS
            species = list(set([mapping[n] for n in node.get_children()]))
            mapping[node] = getLca(specietree, species)
            if(multspeciename):
                node.add_features(
                    species=",".join(sorted([x.name for x in species])))
            else:
                node.add_features(species=mapping[node].name)

    genetree.add_features(lcaMap=mapping)
    return mapping


def reconcile(genetree=None, lcaMap=None, lost=False, lost_label_fn=None):
    """Reconcile genetree topology to a specietree, using an adequate mapping obtained with lcaMapping.
    'reconcile' will infer evolutionary events like gene lost, gene speciation and gene duplication with distinction between AD and NAD
    """

    if(lcaMap is None or genetree is None):
        raise Exception("lcaMapping or genetree not found")
    else:
        lost_count = 1
        for node in genetree.traverse("levelorder"):
            node.add_features(type=TreeClass.SPEC)
            node.add_features(dup=False)

            # print node.name , node.species, " and children name ",
            # node.get_children_name()," and children species ",
            # node.get_children_species()
            if(not node.is_leaf() and (lcaMap[node] == lcaMap[node.get_child_at(0)] or lcaMap[node] == lcaMap[node.get_child_at(1)])):
                node.dup = True
                node.type = TreeClass.AD
                # print "\n\nnode = ", node, "\n\nand children : ",
                # node.children
                if not (set(node.get_child_at(0).get_species()).intersection(set(node.get_child_at(1).get_species()))):
                    node.type = TreeClass.NAD

        if (isinstance(lost, basestring) and lost.upper() == "YES") or lost:
            for node in genetree.traverse("postorder"):
                children_list = node.get_children()
                node_is_dup = (
                    node.type == TreeClass.NAD or node.type == TreeClass.AD)
                for child_c in children_list:
                    if((node_is_dup and lcaMap[child_c] != lcaMap[node]) or (not node_is_dup and (lcaMap[child_c].up != lcaMap[node]))):

                        while((lcaMap[child_c].up != lcaMap[node] and node.type == TreeClass.SPEC) or (lcaMap[child_c] != lcaMap[node] and node.type != TreeClass.SPEC)):
                            lostnode = TreeClass()
                            intern_lost = TreeClass()
                            intern_lost.add_features(type=TreeClass.SPEC)
                            intern_lost.add_features(dup=False)

                            if lcaMap[child_c].is_root():
                                intern_lost.species = ",".join(
                                    lcaMap[child_c].get_leaf_names())
                                lcaMap.update({intern_lost: lcaMap[child_c]})

                            else:
                                intern_lost.species = ",".join(
                                    lcaMap[child_c].up.get_leaf_names())
                                lcaMap.update(
                                    {intern_lost: lcaMap[child_c].up})

                            # change here to display a subtree and not a leaf
                            # with a lot of specie
                            lostnode.species = ",".join(
                                set(lcaMap[intern_lost].get_leaf_names()) - set(lcaMap[child_c].get_leaf_names()))
                            splist = lostnode.species.split(',')
                            if(len(splist) > 1):
                                if lost_label_fn:
                                    lostnode.name = lost_label_fn(splist)
                                else:
                                    lostnode.name = "lost_" + \
                                        str(lost_count) + "_" + \
                                        "|".join([s[0:3] for s in splist])

                            else:
                                if lost_label_fn:
                                    lostnode.name = lost_label_fn(
                                        lostnode.species)
                                else:
                                    lostnode.name = "lost_" + lostnode.species

                            lostnode.add_features(type=TreeClass.LOST)
                            lostnode.add_features(dup=False)

                            lost_count += 1
                            child_c.detach()
                            # print "***********************\n\n** node : ", node, "\n\n** child_c: ", child_c, "\n\n** child parent", child_c.up
                            # node.remove_child(child_c)
                            intern_lost.add_child(child=lostnode)
                            intern_lost.add_child(child=child_c)
                            child_c = intern_lost
                        node.add_child(child_c)
                        children_list.append(child_c)

                # Case of polytomie in species tree....
                if not node.is_leaf():
                    specie_list = ",".join(
                        [",".join(lcaMap[child_c].get_leaf_names()) for child_c in node.get_children()])
                    child_specie_set = set(specie_list.split(","))
                    real_specie_list = set(lcaMap[node].get_leaf_names())
                    unadded_specie = real_specie_list - child_specie_set
                    # print unadded_specie, child_specie_set, real_specie_list
                    # print node.species
                    if(unadded_specie):
                        lostnode = TreeClass()
                        lostnode.add_features(type=TreeClass.LOST)
                        lostnode.add_features(dup=False)
                        lostnode.species = ",".join(unadded_specie)

                        if(len(unadded_specie) > 1):
                            lostnode.name = "lost_" + \
                                str(lost_count) + "_" + \
                                "|".join([s[0:3] for s in unadded_specie])

                        else:
                            lostnode.name = "lost_" + lostnode.species

                        lost_count += 1
                        node.add_child(lostnode)
    genetree.add_features(reconciled=True)


def computeDLScore(genetree, lcaMap=None, dupcost=None, losscost=None):
    """
    Compute the reconciliation cost
    """
    if not lcaMap and genetree.has_feature('lcaMap'):
        lcaMap = genetree.lcaMap
    dup_score = 0
    loss_score = 0
    if lcaMap:
        for node in genetree.traverse("levelorder"):
            node_is_dup = 0
            child_map = [lcaMap[child] for child in node.get_children()]
            if (lcaMap[node] in child_map):
                node_is_dup = params.getdup(lcaMap[node])
                dup_score += (dupcost if dupcost else node_is_dup)

            for child in node.get_children():
                if node_is_dup:
                    child_map = [lcaMap[node]]
                else:
                    child_map = lcaMap[node].get_children()
                curr_node = lcaMap[child]
                while(curr_node not in child_map):
                    lost_nodes = set(
                        curr_node.up.get_children()) - set([curr_node])
                    if losscost:
                        loss_score += len(lost_nodes) * losscost
                    else:
                        loss_score += np.sum([params.getloss(l)
                                              for l in lost_nodes])
                    curr_node = curr_node.up

    else:
        raise Exception("LcaMapping not provided !!")
    return dup_score, loss_score


def computeDTLScore(genetree, speciestree, Dc=1, Tc=1, Lc=1, flag=True):
    if not speciestree.has_feature('lcaprocess', True):
        speciestree.label_internal_node()
        lcaPreprocess(speciestree)
    leafMap = {}
    for leaf in genetree:
        if not leaf.has_feature('species'):
            raise ValueError("You should set species before calling")
        leafMap[leaf] = speciestree&leaf.species
    
    cost_table = MatrixRep(genetree, speciestree, np.inf)
    spec_table = MatrixRep(genetree, speciestree, np.inf)
    dup_table = MatrixRep(genetree, speciestree, np.inf)
    trf_table = MatrixRep(genetree, speciestree, np.inf)
    in_table = MatrixRep(genetree, speciestree, np.inf)
    inAlt_table = MatrixRep(genetree, speciestree, np.inf)
    out_table = MatrixRep(genetree, speciestree, np.inf)
    
    for gleaf in genetree:
        glsmap = leafMap[gleaf]
        cost_table[gleaf, glsmap] = 0
        comp_spec = glsmap
        while comp_spec is not None:
            inAlt_table[gleaf, comp_spec] = 0
            in_table[gleaf, comp_spec] = Lc*(-comp_spec.depth + glsmap.depth)
            comp_spec =  comp_spec.up
    for gnode in genetree.iter_internal_node(strategy="postorder", enable_root=True):
        for snode in speciestree.traverse("postorder"):
            gchild1, gchild2 = gnode.get_children()
            if snode.is_leaf():
                spec_table[gnode, snode] = np.inf
                dup_table[gnode, snode] = Dc + cost_table[gchild1, snode] + cost_table[gchild2, snode]
                # because we can't have transfer at root
                if not snode.is_root():
                    # one child is incomprable and the second
                    # is a descendant    
                    trf_table[gnode, snode] = Tc + min(in_table[gchild1, snode]+out_table[gchild2, snode], in_table[gchild2, snode] + out_table[gchild1, snode])
                
                cost_table[gnode, snode] = min(
                    spec_table[gnode, snode], 
                    dup_table[gnode,snode], 
                    trf_table[gnode, snode]
                    )
                
                in_table[gnode, snode] = cost_table[gnode, snode]
                inAlt_table[gnode, snode] = cost_table[gnode, snode]

            else:
                schild1, schild2 = snode.get_children()
                spec_table[gnode, snode] =  min(
                    in_table[gchild1, schild1] + in_table[gchild2, schild2], 
                    in_table[gchild1, schild2] + in_table[gchild2, schild1]
                    )
                dcost_g_s = 0
                if flag :
                    dcost_g_s = min(
                        cost_table[gchild1, snode] + in_table[gchild2, schild1] + Lc, # loss in one child
                        cost_table[gchild1, snode] + in_table[gchild2, schild2] + Lc,
                        cost_table[gchild2, snode] + in_table[gchild1, schild1] + Lc,
                        cost_table[gchild2, snode] + in_table[gchild1, schild2] + Lc,
                        cost_table[gchild2, snode] + cost_table[gchild1, snode], #both map to snode
                        in_table[gchild1, schild1] + in_table[gchild2, schild1] + 2*Lc, #both map to descendant of snode
                        in_table[gchild1, schild1] + in_table[gchild2, schild2] + 2*Lc, #both map to descendant of snode
                        in_table[gchild1, schild2] + in_table[gchild2, schild2] + 2*Lc, #both map to descendant of snode
                        in_table[gchild1, schild2] + in_table[gchild2, schild1] + 2*Lc, #both map to descendant of snode

                    )
                    
                else: 
                    dcost_g_s = in_table[gchild1, snode] + in_table[gchild2, snode]
                dup_table[gnode, snode] = Dc + dcost_g_s
                if not snode.is_root():
                    trf_table[gnode, snode] = Tc + min(
                        in_table[gchild1, snode]+out_table[gchild2, snode], 
                        in_table[gchild2, snode] + out_table[gchild1, snode]
                        )

                cost_table[gnode, snode] = min(
                    spec_table[gnode, snode],
                    dup_table[gnode, snode],
                    trf_table[gnode, snode]
                )
                in_table[gnode, snode] = min(
                    cost_table[gnode, snode],
                    in_table[gnode, schild1] + Lc,
                    in_table[gnode, schild2] + Lc
                )
                inAlt_table[gnode, snode] = min(
                        cost_table[gnode, snode],
                        inAlt_table[gnode, schild1],
                        inAlt_table[gnode, schild2]
                )

        for snode in speciestree.iter_internal_node("preorder",True):
            schild1, schild2 = snode.get_children()
            out_table[gnode, schild1] = min(
                out_table[gnode,snode],
                inAlt_table[gnode, schild2]
            )
            out_table[gnode, schild2] = min(
                out_table[gnode,snode],
                inAlt_table[gnode, schild1]
            )
    #print cost_table.matrix
    return np.min(cost_table[genetree])

def computeDL(genetree, lcaMap=None):
    """
    Compute the number of duplication and the number of losses
    """
    loss = 0
    dup = 0

    if not lcaMap and genetree.has_feature('lcaMap'):
        lcaMap = genetree.lcaMap

    if lcaMap and not genetree.is_reconcilied():
        for node in genetree.traverse("levelorder"):
            if (lcaMap[node] in [lcaMap[child] for child in node.get_children()]):
                dup += 1
            if(node.up):
                parent = node.up
                parent_is_dup = 0
                if(lcaMap[parent] in [lcaMap[child] for child in parent.get_children()]):
                    parent_is_dup = 1
                loss += (lcaMap[node].depth -
                         lcaMap[parent].depth - 1 + parent_is_dup)

    else:
        if(genetree is None or not genetree.is_reconcilied()):
            raise Exception(
                "LcaMapping not found and your Genetree didn't undergo reconciliation yet")

        for node in genetree.traverse():
            if node.has_feature('type'):
                if(node.type == TreeClass.NAD or node.type == TreeClass.AD):
                    dup += 1
                elif node.type == TreeClass.LOST:
                    loss += 1

    return dup, loss


def cleanFeatures(tree=None, features=[]):
    cleaned = False
    if(tree):
        for node in tree.traverse():
            for f in features:
                if(node.has_feature(f)):
                    node.del_feature(f)
                    cleaned = True
    return cleaned


def __is_dist_elligible(tree):
    """Check whether or not a tree has branch length on all its branch"""
    return not (all([n.dist == 1.0 for n in tree.iter_descendants()]) and tree.dist == 0)


def get_distance_from_tree(tree):
    """Return a distance matrix from input tree
    """
    node_order = tree.get_leaf_names()
    if not __is_dist_elligible(tree):
        raise ValueError(
            "Cannot infer distance matrix from tree branch length. All branch are set to default")
    nl = len(node_order)  # number of leaf
    distance_mat = np.zeros((nl, nl), dtype=float)
    for i in range(nl):
        for j in range(i + 1, nl):
            distance_mat[i, j] = distance_mat[
                j, i] = tree.get_distance(node_order[i], node_order[j])
    np.fill_diagonal(distance_mat, 0)
    return distance_mat, node_order


def binaryRecScore(node, lcamap, dupcost=None, losscost=None):
    """Reconcile genetree topology to a specietree, using an adequate mapping obtained with lcaMapping.
    'reconcile' will infer evolutionary events like gene lost, gene speciation and gene duplication with distinction between AD and NAD
    """
    dup = 0
    lost = 0
    if(lcamap is None or node is None):
        raise Exception("lcaMapping or genetree not found")
    else:
        # print node.name , node.species, " and children name ",
        # node.get_children_name()," and children species ",
        # node.get_children_species()
        if(not node.is_leaf() and (lcamap[node].name == lcamap[node.get_child_at(0)].name or lcamap[node].name == lcamap[node.get_child_at(1)].name)):
            if not dupcost:
                dup += params.getdup(lcamap[node].name)
            else:
                dup += dupcost

        children_list = node.get_children()
        supposed_children_species = lcamap[node].get_children_name()
        child_number = 0
        for child in children_list:
            c = lcamap[child]
            child_number += 1
            child_lost = 0

            if(dup == 0):
                while(c is not None and (c.name not in supposed_children_species)):
                    if not losscost:
                        lost += params.getloss(c.name)
                    else:
                        lost += losscost

                    child_lost += 1
                    c = c.up

            if(dup > 0):
                while(c is not None and c.name != node.species):
                    if not losscost:
                        lost += params.getloss(c.name)
                    else:
                        lost += losscost
                    child_lost += 1
                    c = c.up

    return dup + lost, dup, lost


def totalDuplicationConsistency(tree):
    """Compute the total duplication consistency score for a tree"""
    dps = 0
    if not tree.is_reconcilied():
        raise ValueError("Your tree wasn't reconcilied")
    for node in tree.traverse():
        if(node.type == TreeClass.AD or node.type == TreeClass.NAD):
            try:
                dps += node.compute_dup_cons()
            except AssertionError:
                pass
    return dps


def getTreeFromPhyloxml(xml, saveToFile="default.xml", delFile=True):
    """
    Read a phylogeny tree from a phyloxml string and return a TreeClass object
    or a list of TreeClass object
    """
    project = Phyloxml()
    fo = open(saveToFile, "w+")
    fo.write(xml)
    fo.close()
    project.build_from_file(saveToFile)
    treeList = []
    for tree in project.get_phylogeny():
        treeList.append(TreeClass.import_from_PhyloxmlTree(tree))

    if(delFile):
        os.remove(saveToFile)
    if len(treeList) == 1:
        return treeList[0]
    return treeList


def resetNodeName(tree, sep, spec_pos):
    spec_pos *= -1
    for x in tree.traverse():
        x.name = x.name.split(sep)[spec_pos]
    return tree


def makeRandomTree(names=list(string.lowercase), contract_seuil=0, feature_to_contract='support', random_branches=False):
    """Make a random Gene Tree"""
    tree = TreeClass()
    tree.populate(
        len(names), names_library=names, random_branches=random_branches)
    tree.contract_tree(seuil=contract_seuil, feature=feature_to_contract)
    return tree


def getSpecieCount(tree):
    """Species distribution in the genetree"""
    count = ddict(int)
    for node in tree.get_children():
        count[node.species] += 1
    return count


def getReverseMap(lcamap, use_name=False):
    """Get reverse map from specie to gene"""
    reversedmap = ddict(list)
    for (g, s) in lcamap.items():
        if(use_name):
            reversedmap[s.name].append(g)
        else:
            reversedmap[s].append(g)
    return reversedmap


def getImageTreeNode(genetree, specietree, lcamap):
    """ Get the specie image tree node of a genetree"""

    # get pre(s) for each  node in specietree
    reversedmap = getReverseMap(lcamap)
    k = 0
    # Traverse G in df order and set ih to 0 for internal node
    for node in genetree.iter_internal_node("levelorder", enable_root=True):
        node.add_features(i_h=0)
        node.name = 'n%d' % k

    # Arange the children of each node in G according to the position of their images
    # in post-order traversal of S
    for snode in specietree.traverse("postorder"):
        for gnode in reversedmap[snode]:
            p_gnode = gnode.up
            if(p_gnode):
                gnode_ind = [x for x in xrange(len(p_gnode.children)) if p_gnode.children[
                    x] == gnode][0]
                p_gnode.children[gnode_ind], p_gnode.children[
                    p_gnode.i_h] = p_gnode.children[p_gnode.i_h], p_gnode.children[gnode_ind]
                p_gnode.i_h += 1

    # compute B(s) that contains all the gene tree nodes g / s in I(g) for s
    # in S
    B_array = ddict(list)
    for node in genetree.traverse("postorder"):
        childlist = node.get_children()
        for child in childlist:
            B_array[lcamap[child]].append(node)
        for i in xrange(0, len(childlist) - 1):
            B_array[getLca(specietree, [lcamap[childlist[i]],
                                        lcamap[childlist[i + 1]]])].append(node)

    # Store all the specie tree nodes of the compresses child-image subtree I(g)
    # and construct all I(g)
    # At this step, we are actually certain that the euler tour of S was
    # already computed
    image_tree_nodes = ddict(list)
    for s in specietree.ind2node:
        for h in B_array[s]:
            image_tree_nodes[h].append(s)

    # Here we are going to construct the tree
    image_tree = {}
    for node in genetree.traverse("postorder"):
        nodecopied = {}
        if not image_tree_nodes[node]:
            continue
        el1 = image_tree_nodes[node].pop()
        a = el1._copy_node(features=['name', 'depth'])
        nodecopied[el1] = a
        while len(image_tree_nodes[node]) > 0:
            el2 = image_tree_nodes[node].pop()
            b = nodecopied.get(el2, None)
            if not b:
                b = el2._copy_node(features=['name', 'depth'])
                nodecopied[el2] = b
            if (a != b):
                if a.depth < b.depth:
                    if(b not in a.get_children()):
                        a.add_child(b)
                elif a.depth > b.depth:
                    if(a not in b.get_children()):
                        b.add_child(a)
            a = b
        image_tree[node] = a.get_tree_root()
    return image_tree


def getSpecieGeneMap(genetree, specietree):
    """Find the reversed map (map between specietree node and genetree node)"""
    mapGene = {}
    for node in specietree.traverse():
        mapGene[node] = genetree.search_nodes(species=node.name)

    return mapGene


def treeHash(tree, addinfos=''):
    """Hashing the tree based on the sorted node name"""
    newick_str = re.sub(
        "(?<=\()([^()]+?)(?=\))", lambda m: ",".join(sorted(m.group(1).split(","))), tree.write(format=9))
    # print "newick: ", tree.write(format=9), "parsing: ", newick_str
    return hashlib.sha384(newick_str + addinfos).hexdigest()


def newickPreprocessing(newick, gene_sep=None):
    """Newick format pre-processing in order to assure its correctness"""
    DEF_SEP_LIST = [';;', '-', '|', '%', ':', ';', '+', '/']

    if isinstance(newick, basestring):
        if os.path.exists(newick):
            nw = open(newick, 'rU').read()
        else:
            nw = newick
        nw = nw.strip()
        if nw.endswith(';'):
            nw = nw[:-1]

        if gene_sep is None:
            i = 0
            while i < len(DEF_SEP_LIST) and DEF_SEP_LIST[i] not in nw:
                i += 1
            if i < len(DEF_SEP_LIST):
                gene_sep = '%%'
                nw = nw.replace(DEF_SEP_LIST[i], gene_sep)

            elif i >= len(DEF_SEP_LIST) or ';' in nw:
                raise NewickError('Unable to format your newick file, Bad gene-specie separator or too much special chars')
        nw += ';'
        return nw, gene_sep
    else:
        raise NewickError("'newick' argument must be either a filename or a newick string.")


def polySolverPreprocessing(genetree, specietree, distance_mat, capitalize=False, gene_sep=None, specie_pos="postfix", nFlagVal=1e305, nFlag=False, smap=None, errorproof=False):
    """Preprocess genetree for polytomysolver
    """

    # genetree input
    speciemap = None
    if isinstance(genetree, basestring) and not smap:
        genetree, gene_sep = newickPreprocessing(genetree, gene_sep)
        genetree = TreeClass(genetree)

    elif smap:
        if isinstance(smap, dict):
            speciemap = smap
        else:
            genetree = TreeClass(genetree) if isinstance(
                genetree, basestring) else genetree
            regexmap = {}
            speciemap = {}
            with open(smap, 'rU') if isinstance(smap, basestring) else smap as INPUT:
                for line in INPUT:
                    g, s = line.strip().split()
                    if ('*') in g and '.*' not in g:
                        g = g.replace('*', '.*')
                    g_regex = re.compile(g, re.IGNORECASE)
                    regexmap[g_regex] = s

            for leaf in genetree:
                for key, value in regexmap.iteritems():
                    if key.match(leaf.name):
                        speciemap[leaf.name] = value

    genetree.set_species(
        speciesMap=speciemap, sep=gene_sep, capitalize=capitalize, pos=specie_pos)

    # genetree check
    if len(genetree) != len(set(genetree.get_leaf_names())):
        tmp_leaf_name = genetree.get_leaf_names()
        duplicates = set(
            [x for x in tmp_leaf_name if tmp_leaf_name.count(x) > 1])
        raise ValueError(
            "Your polytomy contains the following gene multiple times : %s" % ", ".join(duplicates))

    # specietree input
    if isinstance(specietree, basestring):
        specietree, sep = newickPreprocessing(specietree, '')
        specietree = TreeClass(specietree)
    specietree.label_internal_node()

    # distance matrice input
    if(distance_mat):
        if isinstance(distance_mat, basestring):
            gene_matrix, node_order = clu.distMatProcessor(
                distance_mat, nFlagVal, nFlag)
        else:
            # distance mat is provided as a boolean
            # in that case, just try to get it from the genetree
            gene_matrix, node_order = get_distance_from_tree(genetree)
        # Difference check 1
        # pos = node_order.index('ENSDORP00000008194_dordii')
        # print node_order
        # print gene_matrix[pos, :]
        listerr = set(node_order).symmetric_difference(
            set(genetree.get_leaf_names()))
        if listerr:
            if not errorproof:
                raise ValueError(
                    "Different genes in distance matrix and genetree\n : See symmetric difference : %s\n" % ", ".join(listerr))
            else:
                if gene_sep:
                    resetNodeName(genetree, gene_sep, specie_pos == 'postfix')
                else:
                    exib1 = set(node_order).difference(
                        set(genetree.get_leaf_names()))
                    exib2 = set(genetree.get_leaf_names()
                                ).difference(set(node_order))
                    if exib2:
                        raise Exception(
                            'Genes in trees and not in matrix : %s' % (exib2))
                    elif exib1:
                        print("Genes in matrix and not in tree : %s \nAttempt to correct distance matrix" % (
                            ", ".join(exib1)))
                        for l in exib1:
                            try:
                                lpos = node_order.index(l)
                                gene_matrix = clu.remove_ij(
                                    gene_matrix, lpos, lpos)
                                del node_order[lpos]
                            except:
                                raise IndexError(
                                    "Could not remove gene %s from distance matrix" % l)

    else:
        # This is for debug, will never happen
        raise ValueError(
            "distance matrix not provided and could not be infered from tree")
        # gene_matrix = clu.makeFakeDstMatrice(len(node_order), 0, 1)

    # Find list of species in genetree but not in specietree
    specieGeneList = set(genetree.get_leaf_species())
    specieList = set([x.name for x in specietree.get_leaves()])
    if(specieGeneList - specieList):
        if len(specieGeneList.intersection(specieList)) == 0 and gene_sep:
            raise Exception(
                "*** You probably didn't set the correct species position for you input tree !!")
        raise Exception("Species in genetree but not in specietree : %s" % (
            ", ".join(specieGeneList - specieList)))

    return genetree, specietree, gene_matrix, node_order


def exportToOrthoXML(t, database='customdb', handle=sys.stdout):
    """ This function takes a TreeClass instance and export all
    its speciation and duplication events to the OrthoXML format.

    """

    # Creates an empty orthoXML object
    O = orthoxml.orthoXML()

    # Generate the structure containing sequence information
    leaf2id = {}
    sp2genes = {}
    for genid, leaf in enumerate(t.iter_leaves()):
        spname = leaf.species
        if spname not in sp2genes:
            sp = orthoxml.species(spname)
            db = orthoxml.database(name=database)
            genes = orthoxml.genes()
            sp.add_database(db)
            db.set_genes(genes)
            sp2genes[spname] = genes
            # add info to the orthoXML document
            O.add_species(sp)
        else:
            genes = sp2genes[spname]

        gn = orthoxml.gene(protId=leaf.name, id=genid)
        leaf2id[leaf] = genid
        genes.add_gene(gn)

    # Add an ortho group container to the orthoXML document
    ortho_groups = orthoxml.groups()
    O.set_groups(ortho_groups)

    # OrthoXML does not support duplication events at the root
    # of the tree, so we search for the top most speciation events in
    # the tree and export them as separate ortholog groups
    for speciation_root in t.iter_leaves(is_leaf_fn=(lambda n: getattr(n, 'type', "") == "S" or not n.children)):
        # Creates an orthogroup in which all events will be added
        node2event = {}
        node2event[speciation_root] = orthoxml.group()
        ortho_groups.add_orthologGroup(node2event[speciation_root])

        # if root node is a leaf, just export an orphan sequence within the
        # group
        if speciation_root.is_leaf():
            node2event[speciation_root].add_geneRef(
                orthoxml.geneRef(leaf2id[speciation_root]))

        # otherwise, descend the tree and export orthology structure
        for node in speciation_root.traverse("preorder"):
            if node.is_leaf():
                continue
            parent_event = node2event[node]
            for ch in node.children:
                if ch.is_leaf():
                    parent_event.add_geneRef(orthoxml.geneRef(leaf2id[ch]))
                else:
                    node2event[ch] = orthoxml.group()

                    if not (ch.has_feature('type') or ch.has_feature('dup')):
                        raise AttributeError(
                            "\n\nUnknown evolutionary event. %s" % ch.get_ascii())

                    if(ch.type == TreeClass.SPEC):
                        parent_event.add_orthologGroup(node2event[ch])
                    elif ch.type in TreeClass.DUP:
                        parent_event.add_paralogGroup(node2event[ch])
                    else:
                        raise AttributeError(
                            "\n\Internals nodes labeled by losses are not expected in the orthoXML format")

    O.export(handle, 0, namespace_="")


def generateSmap(specietree, output="smap", relaxed=False, suffix=""):
    """
    Generate a specie map from genetree and specietree
    """
    gene_to_spec_map = []
    specie_names = specietree.get_leaf_names()
    for name in specie_names:
        if(relaxed):
            genes = re.compile(".*" + name + suffix + ".*", re.IGNORECASE)
        else:
            genes = re.compile("^" + name + suffix + ".*", re.IGNORECASE)
        gene_to_spec_map.append([genes.pattern, name])
    with open(output, "w") as f:
        f.writelines('\t'.join(line) + "\n" for line in gene_to_spec_map)


def customTreeCompare(original_t, corrected_t, t):

    # Leaves remaining test and original binary node test
    ct_leaves = []
    t_leaves = []
    t_binary = []
    success = []
    ct_binary = []
    for node in original_t.traverse("levelorder"):
        desc_name = set(node.get_leaf_names())
        ct_parent = corrected_t.get_common_ancestor(desc_name)
        t_parent = t.get_common_ancestor(desc_name)
        ctl = set(ct_parent.get_leaf_names())
        tl = set(t_parent.get_leaf_names())
        ct_leaves.append(ctl.difference(desc_name))
        t_leaves.append(tl.difference(desc_name))
        if(node.is_binary() and not node.has_polytomies()):
            ct_binary.append(ct_parent.robinson_foulds(node)[0:3])
            t_binary.append(t_parent.robinson_foulds(node)[0:3])
        success.append(len(tl.difference(ctl)) < 1)

    ct_success = filter(None, map(lambda x: len(x) < 1, ct_leaves))
    t_success = filter(None, map(lambda x: len(x) < 1, t_leaves))

    print("\nCorrected Tree binary list rf_fould\n")
    print("\n".join(map(lambda x: "\t".join([str(v) for v in x]), ct_binary)))
    print("\nTree binary list rf_fould\n")
    print("\n".join(map(lambda x: "\t".join([str(v) for v in x]), t_binary)))

    if(len(ct_success) == len(ct_leaves)):
        print("**Leave remaining success for corrected tree")
        print("\n".join([str(h) for h in t_success]))

    else:
        print("**Corrected tree doesn't follow patern")
        print("\n".join(map(lambda x: "\t".join(
            [str(v) for v in x]), ct_leaves)))

    if(len(t_success) == len(t_leaves)):
        print("**Leave remaining success for tree")
        # print "\n".join([str(h) for h in t_success])
    else:
        print("**Tree doesn't follow patern")
        # print "\n".join(map(lambda x: "\t".join([str(v) for v in x]),
        # t_leaves))

    print("**Compatibility test between tree: ", all(success))
