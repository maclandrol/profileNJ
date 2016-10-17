# This file is part of profileNJ
#
# Date: 02/2014
# TreeClass is a python class derived from TreeNode class from the ete3 package.
# TreeClass add additional specific function

__author__ = "Emmanuel Noutahi"

from ete3 import TreeNode
from ete3.phylo import EvolEvent
import types
from collections import defaultdict as ddict
from itertools import izip
import copy

try:
    import cPickle as pickle
except:
    import pickle


class TreeClass(TreeNode):

    DEFAULT_SPECIE = "Unknown"
    DEFAULT_NAME = ""
    DEFAULT_GENE = "Unknown"
    AD = 1
    LOST = -1
    SPEC = 0
    NAD = 2
    DUP = [AD, NAD]

    def __init__(self, newick=None, format=0, dist=None, support=None, name=None):
        """	Default init for the TreeClass. This works better than wrapping the entire class"""
        TreeNode.__init__(
            self, newick=newick, format=format, dist=dist, support=support, name=name)

    def __repr__(self):
        return "Tree Class '%s' (%s)" % (self.name, hex(self.__hash__()))

    def get_child_at(self, i=0):
        """Return child at a specific position in the children list of a node"""
        children_list = self.get_children()
        if(i < len(children_list)):
            return children_list[i]
        else:
            raise IndexError(
                "Index out of bound, Can't access the child at index: %i" % i)

    def copy(self, method="simplecopy", nw_features=[], nw_format_root_node=True, binary_correct=False):
        """ .. Override of ete TreeNode original copy

        Returns a copy of the current node.

        :var cpickle method: Protocol used to copy the node
        structure. The following values are accepted:

                - "newick": Tree topology, node names, branch lengths and
                branch support values will be copied by as represented in
                the newick string (copy by newick string serialisation).

                - "newick-extended": Tree topology and all node features
                will be copied based on the extended newick format
                representation. Only node features will be copied, thus
                excluding other node attributes. As this method is also
                based on newick serialisation, features will be converted
                into text strings when making the copy.

                - "cpickle": The whole node structure and its content is
                cloned based on cPickle object serialisation (slower, but
                recommended for full tree copying)

                - "deepcopy": The whole node structure and its content is
                copied based on the standard "copy" Python functionality
                (this is the slowest method but it allows to copy complex
                objects even if attributes point to lambda functions,
                etc.)
                - "simplecopy" : Simple recursive tree topology and feature copy

        """
        # nw_features = ["name", "species"]
        if method == "newick":
            new_node = self.__class__(
                self.write(features=["name"], format_root_node=True))
        elif method == "newick-extended":
            new_node = self.__class__(
                self.write(features=nw_features, format_root_node=nw_format_root_node))
        elif method == "deepcopy":
            parent = self.up
            self.up = None
            new_node = copy.deepcopy(self)
            self.up = parent
        elif method == "cpickle":
            parent = self.up
            self.up = None
            new_node = pickle.loads(pickle.dumps(self, 2))
            self.up = parent

        elif method == 'simplecopy':
            parent = self.up
            self.up = None
            new_node = self._recur_simple_copy(nw_features)
            self.up = parent

        else:
            raise ValueError("Invalid copy method")

        return self._correct_copy(new_node) if binary_correct else new_node

    def _recur_simple_copy(self, features=[]):
        """Simple copy of a node by a recursive call"""
        root = self._copy_node(features)
        for node in self.get_children():
            root.add_child(node._recur_simple_copy(features))
        return root

    def _iter_simple_copy(self, features=[]):
        """Iteratif simple copy, this can be optimized"""
        ori_parents = [self]
        copy_parents = [self._copy_node(features)]

        while ori_parents:
            next_parent = []
            next_copy_parent = []
            for i in xrange(len(ori_parents)):
                parent = ori_parents[i]
                root = copy_parents[i]
                for node in parent.get_children():
                    copy_node = node._copy_node(features)
                    root.add_child(copy_node)
                    next_parent.append(node)
                    next_copy_parent.append(copy_node)
            ori_parents = next_parent
            copy_parents = next_copy_parent
        return root.get_tree_root()

    def _copy_node(self, features=[]):
        "Copy a node and its features to a new node"
        copy = TreeClass()
        if not features:
            features = self.features
        for feature in features:
            if(self.has_feature(feature)):
                copy.add_feature(feature, getattr(self, feature))
        return copy

    def edge_exist(self, node):
        """ Return True if there exist an edge between 2 node"""
        return (self.up == node or node.up == self)

    def insert_node_between(self, node, new_node):
        """ insert a new node between self and node"""
        if not self.edge_exist(node):
            return None
        else:
            if(self.up != node):
                self, node = node, self
            self.detach()
            node.add_child(new_node)
            new_node.add_child(self)
            return node

    def replace_by(self, new_node):
        """Replace self by new node"""
        self.children = new_node.children
        for node in self.children:
            node.up = self
        self.features = set([])
        self.name = new_node.name
        for f in new_node.features:
            self.add_feature(f, getattr(new_node, f))

    def _correct_copy(self, copy):
        """Correct the structure of new node copied using newick method"""
        for node in copy.traverse("postorder"):
            if not node.is_root() and (node.name not in list(self.get_descendant_name())):
                node.detach()
            if node.is_internal() and len(node.get_children()) < 2:
                child = node.get_child_at(0)
                if(child.is_leaf()):
                    ori_node = (self & (child.name)).up
                else:
                    ori_node = self.get_common_ancestor(
                        child.get_leaf_name()).up
                node.up.replace_child(node, ori_node.copy("newick-extended"))
        return copy

    def has_ancestor(self, ancestor):
        """Check if the node in `ancestor` are an ancestor of the current Node"""
        ancestor = self.get_tree_root().translate_nodes(ancestor)
        ancestors = self.get_ancestors()
        return True if len(list(filter(lambda x: x in ancestors, ancestor))) == len(ancestor) else False

    def has_descendant(self, descendant):
        """Check if the nodes is `descendant` are a descendant of the current Node"""
        descendant = self.get_tree_root().translate_nodes(descendant)
        descendants = self.get_descendants()
        return True if len(list(filter(lambda x: x in descendants, descendant))) == len(descendant) else False

    def _euler_visit(self, node_visited=[]):
        """Perform a euler tour and return visited nodes"""
        if self:
            node_visited.append(self)
            self.add_features(euler_visit=True)
            for child in self.get_children():
                if(not child.has_feature('euler_visit')):
                    node_visited = child._euler_visit(node_visited)
        if (self.up):
            node_visited.append(self.up)
        return node_visited

    def translate_nodes(self, *target_nodes):
        """Translate list of node name into Node"""

        if len(target_nodes) == 1 and type(target_nodes[0]) in set([set, tuple, list, frozenset]):
            target_nodes = target_nodes[0]

        try:
            target_nodes = [
                n if isinstance(n, self.__class__) else self & n for n in target_nodes]
            return target_nodes

        except (ValueError, IndexError) as e:
            print("You may have names which cannot be associated with a node")
            raise e

    def insert_child_at(self, i, newNode, replace_if_exist=True):
        """Insertion of a node in a specific position"""
        removing_child = None
        try:
            if i < len(self.get_children()):
                removing_child = self.get_children()[i]
                if(replace_if_exist):
                    self.children[i] = newNode
                    newNode.up = self
            else:
                self.add_child(newNode)
        except ValueError, e:
            raise e
        else:
            return removing_child

    def remove_child_at(self, i):
        """Remove a child at a specific position"""

        try:
            child = self.children[i]
            self.children.remove(child)
        except (IndexError, ValueError), e:
            raise e
        else:
            child.up = None
            return child

    def delete_leaf(self, leafList):
        """Delete a list of leaf"""
        if len(leafList) == 1 and type(leafList[0]) in set([set, tuple, list, frozenset]):
            leafList = leafList[0]

        for leaf in leafList:
            if(leaf.is_leaf()):
                leaf.delete()

    def get_degree(self, as_graph=True):
        """Return the degree of the current Node"""
        child_number = len(self.get_children())
        if(as_graph):
            child_number += (1 if not self.is_root() else 0)
        return child_number

    def set_species(self, speciesMap=None, sep="_", capitalize=False, pos="postfix", use_fn=None, **kwargs):
        """Set species feature for each leaf in the tree.

        :argument speciesMap: Default=None. speciesMap is a Map of species for the geneTree. Each key is a leaf name from the genetree and the value is the corresponding specie name
        :argument sep: Default ="_" , the separator for the default species extraction using the leaf name
        :argument pos: Default="postfix", the species position in the leaf name for the default extraction. Should be used with sep. Can take for value, "prefix", which
        means "specie-sep-gene" or "postfix" for "gene-sep-specie"
        argument fn: Pointer to a parsing python function that receives a node as first argument and returns the species name.

        """
        if speciesMap:
            for node in self.traverse():
                node_specie = speciesMap.get(
                    node.name, TreeClass.DEFAULT_SPECIE)
                node_specie = self.__class__._capitalize(
                    node_specie) if capitalize else node_specie
                node.add_features(species=node_specie)
        else:
            for leaf in self:
                if use_fn is not None:
                    leaf.add_features(species=use_fn(leaf, **kwargs))
                else:
                    leaf.add_features(
                        species=leaf._extract_feature_name(separator=sep, order=pos, cap=capitalize))

    def set_genes(self, genesMap=None, sep="_", capitalize=False, pos="postfix", use_fn=None, **kwargs):
        """Set gene feature for each leaf in the tree.

        :argument genesMap: Default=None. genesMap is a Map of genes for the geneTree. Each key is a leaf name from the genetree and the value is the corresponding genes name
        :argument sep: Default ="_" , the separator for the default genes extraction using the leaf name
        :argument pos: Default="postfix", the gene position in the leaf name for the default extraction. Should be used with sep. Can take for value, "postfix", which
        means "specie-sep-gene" or "prefix" for "gene-sep-specie"
        argument fn: Pointer to a parsing python function that receives a node as first argument and returns the genes name.

        """
        for leaf in self:
            if genesMap:
                node_gene = genesMap.get(leaf.name, TreeClass.DEFAULT_GENE)
                node_gene = self.__class__._capitalize(
                    node_gene) if capitalize else node_gene
                node.add_features(genes=node_gene)

            elif use_fn is not None:
                leaf.add_features(genes=use_fn(leaf, **kwargs))
            else:
                leaf.add_features(
                    genes=leaf._extract_feature_name(separator=sep, order=pos, cap=capitalize))

    def get_species(self, sep=","):
        """Return the list of species for the current node
        (under the current node after reconciliation)"""
        return self.species.split(sep)

    def get_genes(self, sep=","):
        """Return the list of genes for the current node
        (under the current node after reconciliation)"""
        return self.genes.split(sep)

    def _extract_feature_name(self, separator=None, order=None, cap=False):
        """Private function, extract feature name (e.g. genes, species) based on the node name"""
        l = self.name.split(separator)
        if len(l) > 1 and order == "postfix":
            feature = l[-1]
        elif len(l) > 1 and order == "prefix":
            feature = l[0]
        else:
            feature = self.name
        if cap:
            feature = self.__class__._capitalize(feature)
        return feature

    def contract_tree(self, seuil=0, feature='support', break_tree_topo=False):
        """ Contract
         tree based on the dist between node, using a threshold. `contract_tree`
        proceed bottom-up. Any branches with a support less than "seuil" will be removed
        if `break_tree_topo` is set to True, all the branch under this node will be recursively removed
        """
        for node in self.traverse("postorder"):
            if(node.has_feature(feature) and node.is_internal() and getattr(node, feature) < seuil):
                node.to_polytomy(break_tree_topo)

    def restrict_to_species(self, species=[]):
        """Restrict the current genetree to the list of species passed in argument"""
        hits = []
        try:
            for value in species:
                hits.extend(self.get_leaves_by_feature(species=value))
            self.prune(hits)
        except Exception as e:
            print("Check if this tree have species as feature")
            raise e

    @classmethod
    def get_path_to_ancestor(cls, node, ancestor):
        """ Get pathway to ancestor."""
        assert node.has_ancestor(
            ancestor), "ancestor is not an actual ancestor in your tree"
        nodes_between = []
        current = node
        while current != ancestor:
            nodes_between.append(current)
            current = current.up
        nodes_between.append(ancestor)
        return nodes_between

    def get_nodes_between(self, node):
        """Get the list of nodes between two nodes"""
        assert isinstance(node, TreeNode), "Node should be an instance of Tree"
        com_anc = self.get_common_ancestor(node)
        path1 = self.get_path_to_ancestor(self, com_anc)
        path2 = self.get_path_to_ancestor(node, com_anc)[::-1]
        return path1[1:] + path2[1:-1]

    def to_polytomy(self, break_tree_topo=False):
        """Move every leaves to the node by deleting all the internals nodes"""
        if(break_tree_topo):
            for node in self.traverse():
                if(node.is_internal()):
                    node.delete()
        else:
            self.delete()

    def to_star(self):
        """Create a star tree from a tree topology"""
        self.to_polytomy(break_tree_topo=True)

    def get_leaves_by_feature(self, **condition):
        """Return leaves that match the features passed as argument"""
        match = self.search_nodes(**condition)
        return [node for node in match if node.is_leaf()]

    def is_polytomy(self):
        """
        Return True if current node is a polytomy.
        """
        return len(self.children) > 2

    def tree_is_polytomy(self):
        """
        Return True if the all tree is a polytomy.
        """
        root = self.get_tree_root()
        return len(root) == len(root.get_children())

    def is_binary(self):
        """
        Return True if current node is a binary node.
        """
        return self.get_degree(False) == 2

    def is_internal(self):
        """
        Return True if current node is an internal node.
        """
        return not (self.is_root() or self.is_leaf())

    def is_reconcilied(self):
        """
        Return whether or not, this genetree is reconcilied
        """
        return self.get_tree_root().has_feature('reconciled', True)

    def get_internal_node(self, strategy="levelorder", enable_root=False):
        """
        Return the list of all internal nodes under the current node
        """
        internal_nodes = []
        for n in self.traverse(strategy=strategy):
            if n.is_internal() or (enable_root and n.is_root()):
                internal_nodes.append(n)
        return internal_nodes

    def iter_internal_node(self, strategy="postorder", enable_root=False):
        """
        Returns an iterator over the list of internal node under the current node
        """
        for n in self.traverse(strategy=strategy):
            if n.is_internal() or (enable_root and n.is_root()):
                yield n

    def get_all_features(self):
        """Return all the features of all nodes under self in a set"""
        features_list = []
        for node in self.traverse():
            features_list.extend(list(node.features))
        return set(features_list)

    def has_feature(self, feature, name=None):
        """Return weither or not this node has feature in its list of features"""
        return (feature in self.features and (name in [None, self.__getattribute__(feature)]))

    def reroot(self, root_node=True):
        """reroot tree at each node"""
        # self.label_internal_node()
        for node in self.iter_descendants():
            c_tree = self.copy("simplecopy", nw_format_root_node=True)
            # c_node =c_tree&node.name
            c_node = c_tree.get_common_ancestor(
                node.get_leaf_name()) if node.is_internal() else c_tree & node.name
            c_tree.set_outgroup(c_node)
            # case where we root at the node and not at the branch
            if(root_node and not node.is_leaf()):
                root = c_tree.get_tree_root()
                # new_child= [child for child in root.get_children() if child !=node.name][0]
                # rooting_node = [child for child in root.get_children() if child.name ==node.name][0]
                new_child = [child for child in root.get_children() if set(
                    child.get_leaf_name()).symmetric_difference(set(node.get_leaf_name()))][0]
                rooting_node = [
                    child for child in root.get_children() if child != new_child][0]
                c_tree = rooting_node.detach()
                new_child.detach()
                # new_child.label_internal_node()
                c_tree.add_child(new_child)
            yield c_tree

    def iter_edges(self):
        """ Iter all over the edges in this tree"""

        for node in self.traverse("levelorder"):
            for child in node.get_children():
                yield (node, child)

    def get_edges(self):
        """ Return the list of edges in current tree"""
        return [edge for edge in self.iter_edges()]

    def edge_reroot(self, unroot=False):
        """ Reroot a tree around each of its edge"""
        i = 0
        for edge in self.iter_edges():
            parent, child = edge
            c_tree = self.copy("simplecopy", nw_format_root_node=True)
            if(unroot):
                c_tree.unroot()
            c_child = c_tree.get_common_ancestor(
                child.get_leaf_name()) if child.is_internal() else c_tree & child.name

            c_parent = c_child.up
            new_root = TreeClass()
            i += 1
            if c_parent and (unroot or c_parent is not c_tree):
                c_child = c_child.detach()
                path_to_root = [c_parent] if (
                    unroot and c_parent is c_tree) else c_tree.get_path_to_ancestor(c_parent, c_tree)
                sisters = [
                    c_tree] if unroot else path_to_root[-2].get_sisters()

                if(len(path_to_root) > 1):
                    last_node = path_to_root[-1]
                    removed_children = []
                    for n in reversed(path_to_root[:-1]):
                        last_node = last_node.remove_child(n)
                        removed_children.append(last_node)

                    cur_node = removed_children.pop(-1)
                    for node in reversed(removed_children):
                        cur_node = cur_node.add_child(node)

                    # adding sister list at the root
                    for c in sisters:
                        cur_node.add_child(c)

                new_root.add_child(c_child)
                new_root.add_child(c_parent)
                # additional security to avoid single internal node
                new_root.delete_single_child_internal()
                yield new_root

    def get_events(self, include_lost=True):
        """Returns a list of **all** duplication and speciation
        events detected after this node.
        """
        all_events = []
        root = self.copy()
        assert self.is_reconcilied(), "Your tree is not reconciled!"
        for node in self.traverse("levelorder"):
            e = EvolEvent()
            e.node = node

            if not node.is_leaf():
                assert node.is_binary(), node
                species_include = [set([n.name for n in node.get_child_at(0).get_leaves() if n.type != TreeClass.LOST]),
                                   set([n.name for n in node.get_child_at(1).get_leaves() if n.type != TreeClass.LOST])]

                if(len(species_include[0]) > 0 and len(species_include[1]) > 0):
                    if node.type > 0:
                        e.etype = 'D'
                        e.dup_score = node.compute_dup_cons()
                        e.paralogs = species_include
                    else:
                        e.etype = 'S'
                        e.orthologs = species_include

            elif(node.has_feature('type', TreeClass.LOST) and include_lost):
                e.etype = 'L'
                e.losses = node.species

            if(e.etype):
                all_events.append(e)

        return all_events

    def _find_closest_descendant_having_feature(self, name, value, descendant=[]):
        """Find the first closest descendant in right and left child
        having a particular feature """

        for node in self.get_children():
            if(node.has_feature(name, value)):
                descendant.append(node)
            else:
                descendant = node._find_closest_descendant_having_feature(
                    name, value, descendant)

        return descendant

    def are_paralogous(self, genelist):
        """Return true if the gene in genelist descend from are paralogous
        """
        if isinstance(genelist, TreeNode):
            genelist = genelist.get_leaf_name()
            genelist.extend(self.get_leaf_name())
            genelist = set(genelist)
        elif len(genelist) == 1 and self.is_leaf():
            genelist.append(self.name)

        if len(genelist) < 2:
            raise TypeError(
                "GeneList incorrect, should be a list of gene name or a node")

        root = self.get_tree_root()
        com_anc = self.get_common_ancestor(genelist)
        paralogous = com_anc.has_feature('dup', True)
        if paralogous:
            return paralogous, [set(genelist)]
        else:
            duplicated_descendants = root._find_closest_descendant_having_feature(
                'dup', True)
            return paralogous, [set(node.get_leaf_names()) for node in duplicated_descendants]

    def is_monophyletic(self, specieSet):
        """ Returns True if species names under this node are all
        included in a given list or set of species names."""

        if type(specieSet) != set:
            specieSet = set(specieSet)
        return self.get_leaf_species().issubset(specieSet)

    def has_polytomies(self):
        """Return whether or not this tree has polytomies
        """
        return len(self.get_polytomies()) > 0

    def get_children_species(self):
        """ Return the species list of the children under this particular node
        """
        c_species = set([])
        for node in self.get_children():
            c_species.add(node.species)
        return c_species

    def get_children_name(self):
        """ Return the names of the children under this particular node
        """
        c_names = set([])
        for node in self.get_children():
            c_names.add(node.name)
        return c_names

    def get_descendant_species(self):
        """ Return the species list of the descendants under this particular node
        """
        c_species = set([])
        for node in self.get_descendants():
            c_species.add(node.species)
        return c_species

    def get_leaf_species(self, is_leaf_fn=None):
        """ Return the species list of the leave under this node
        """
        return set([leaf.species for leaf in self.iter_leaves(is_leaf_fn=is_leaf_fn)])

    def get_descendant_name(self):
        """ Return the names of the descendants under this particular node
        """
        c_names = set([])
        for node in self.get_descendants():
            c_names.add(node.name)
        return c_names

    def get_ancestor_name(self):
        """Return the names of all the ancestor of this node in a set
        """
        c_names = set([])
        parent = self.up
        while parent is not None:
            c_names.add(parent.name)
            parent = parent.up
        return c_names

    def get_leaf_name(self, is_leaf_fn=None):
        return self.get_leaf_names(is_leaf_fn)

    def delete_single_child_internal(self):
        for node in self.traverse("postorder"):
            if(node.is_internal() and len(node.get_children()) < 2):
                node.delete()

    def has_single_child_internal(self):
        for node in self.traverse("postorder"):
            if(node.is_internal() and len(node.children) < 2):
                return True
        return False

    def iter_polytomies(self, is_polytomy_fn=None, strategy="postorder"):
        """
        Returns an iterator over the polytomies starting from the curent node
        :argument None is_polytomy_fn: See :func:`TreeNode.traverse` for
        documentation.
        """
        for n in self.traverse(strategy=strategy):
            if not is_polytomy_fn:
                if n.is_polytomy():
                    yield n
            else:
                if is_polytomy_fn(n):
                    yield n

    def get_polytomies(self, ind=-1, is_polytomy_fn=None):
        """
        Return a list of polytomies under this node
        """
        polytomies = [
            pol for pol in self.iter_polytomies(is_polytomy_fn=is_polytomy_fn)]
        if(ind > -1 and ind < len(polytomies)):
            return polytomies[ind]
        else:
            return polytomies

    def label_internal_node(self):
        """Label the internal node of a specietree for the polysolver algorithm"""
        count = 1
        for node in self.traverse(strategy='levelorder'):
            if not node.is_leaf() and node.name in [TreeClass.DEFAULT_NAME, TreeClass.DEFAULT_GENE, '']:
                node.name = "n%i" % (count)
                count += 1
        return self

    def get_feature_sum(self, feature):
        """Sum all values of feature for all nodes"""
        cost = 0
        for node in self.traverse("levelorder"):
            if(node.has_feature(feature) and getattr(node, feature) is int):
                cost += getattr(node, feature)
        return cost

    def compute_dup_cons(self):
        """Compute duplication consitency score at the node,
        this function will raise an error if the node is a polytomy or a speciation node
        """
        assert(not self.is_leaf() and self.is_binary() and (
            self.type > 0 or self.has_feature('dup', True)))  # self should be a duplication node
        r_child_spec_set = self.get_child_at(0).get_leaf_species()
        l_child_spec_set = self.get_child_at(1).get_leaf_species()
        inter_set = r_child_spec_set.intersection(l_child_spec_set)
        union_set = r_child_spec_set.union(l_child_spec_set)
        self.add_feature('dupcons', len(inter_set) / len(union_set))
        return self.dupcons

    @classmethod
    def import_from_PhyloxmlTree(cls, phyloxml):
        """import Tree structure and useful Tree features from a _Phyloxml.PhyloxmlTree to a TreeClass
        This is really dirty but it does the job!!.
        ****For each attribut, only the first xml value in the clade is transferred to the TreeClass instance***
        accessible feature: Most of the features are list() or dict(). Because phyloxml format support multi Tree and we can have muliple infos per node!
        -code : the taxon code at the node
        -sc_name: scientific_name for the current node, if it was retrieved from emsembl
        -c_name:common_name for the current node
        -taxon_id: taxon_id from the node
        -'type' : types of the sequence (cdna/protein...
        -'symbol': symbol of the gene
        -'name': name of the gene
        -'accession': accession number of the gene
        -'mol_seq': aa sequence or nuc sequence
        """
        if(phyloxml.__class__.__name__ != "PhyloxmlTree"):
            raise ValueError("Please provide a phyloxml class")

        for node in phyloxml:
            clade = node.phyloxml_clade
            sequence = ddict(list)
            taxa = ddict(list)
            for seq in clade.get_sequence():
                sequence['name'].append(seq.get_name())
                sequence['symbol'].append(seq.get_symbol())
                if(seq.get_mol_seq() is not None):
                    sequence['mol_seq'].append(
                        seq.get_mol_seq().get_valueOf_())
                sequence['type'].append(seq.get_type())
                if(seq.get_accession() is not None):
                    sequence['accession'].append(
                        seq.get_accession().get_valueOf_())

            for taxon in clade.get_taxonomy():
                taxa['common_name'].append(taxon.common_name)
                taxa['sc_name'].append(taxon.scientific_name)
                taxa['code'].append(taxon.code)
                if(taxon.id is not None):
                    taxa['taxon_id'].append(taxon.id.get_valueOf_())

            if(len(taxa['code']) >= 1):
                node.add_features(taxacode=taxa['code'][0])
            if(len(taxa['sc_name']) >= 1):
                node.add_features(specie_scienc=taxa['sc_name'][0])
            if(len(taxa['common_name']) >= 1):
                node.add_features(specie_common=taxa['common_name'][0])
            if(len(taxa['taxon_id']) >= 1):
                node.add_features(tax_id=taxa['taxon_id'][0])

            if(len(sequence['accession']) >= 1):
                node.add_features(alt_name=node.name)
                node.add_features(name=sequence['accession'][0])
            if(len(sequence['mol_seq']) >= 1):
                node.add_features(sequence=sequence['mol_seq'][0])

            if(len(sequence['name']) >= 1):
                node.add_features(seqname=sequence['name'][0])
            if(len(sequence['symbol']) >= 1):
                node.add_features(symbol=sequence['symbol'][0])

        return TreeClass(phyloxml.write(features=[], format_root_node=True))

    def replace_child(self, old_child, new_child):
        """ Replace a child by another node in a tree"""
        if (self is None) or (old_child not in self.get_children()):
            raise ValueError(
                "Node is None or old_child is not a child of current node")
        else:
            self.remove_child(old_child)
            self.add_child(new_child)
            return self

    def write_seq_to_fasta(self, out='seq.fasta', comment=1):
        """Save sequence in tree into a fasta file"""
        if("sequence" in self.get_all_features()):
            with open(out, 'w') as outfile:
                for leaf in self:
                    if("sequence" in leaf.features):
                        id = ">%s" % leaf.name
                        if(comment and "alt_name" in leaf.features and "seqname" in leaf.features):
                            id = id + " %s;%s" % (leaf.alt_name, leaf.seqname)
                            if("sc_name" in leaf.features):
                                id = id + (";%s" % leaf.sc_name)
                        seq = leaf.sequence + "\n"
                        id = id + "\n"
                        outfile.write(id)
                        outfile.write(seq)

    @staticmethod
    def _capitalize(line):
        return "".join([line[0].upper(), line[1:]])
