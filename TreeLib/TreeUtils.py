"""
Author: Emmanuel Noutahi
Date: 02/2014
TreeUtils is a python class that offer function related to phylogeny tree, using TreeClass
"""
import ClusterUtils as clu
import hashlib
import os
import params
import re
import string
import sys
import urllib2
import numpy as np
from TreeClass import TreeClass
from collections import defaultdict as ddict
from ete2 import Phyloxml
from ete2 import orthoxml
from ete2.parser.newick import NewickError


# TreeUtils:

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
        #http = httplib2.Http(".cache")
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
        #resp, content = http.request(server+ext, method="GET", headers={"Content-Type":output})
        if not resp.getcode() == 200:
            print "Invalid response: ", resp.getcode()
            raise ValueError('Failled to process request!')

        if(output.lower() != "text/x-phyloxml"):
            return TreeClass(content)
        else:
            return getTreeFromPhyloxml(content)


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
            print "Invalid response: ", resp.status
            raise ValueError('Failled to process request!')
        if(output.lower() != "text/x-phyloxml"):
            return TreeClass(content)
        else:
            return getTreeFromPhyloxml(content)


def lcapreprocess(tree):
    """Make an euler tour of this tree"""
    #root = tree.get_tree_root()
    tree.add_features(depth=0)
    tree.del_feature('euler_visit')
    for node in tree.iter_descendants("levelorder"):
        node.del_feature('euler_visit')
        # can write this because parent are always visited before
        node.add_features(depth=node.up.depth + 1)
    node_visited = tree._euler_visit([])
    #print tree.get_ascii(show_internal= True, attributes=['name', 'euler_visit', 'depth'])
    # number of element in array
    n = len(node_visited)
    m = int(np.ceil(np.log2(n)))
    rmq_array =  np.zeros((n, m), dtype=int)
    for i in xrange(n):
        rmq_array[i,0] = i
    for j in xrange(1, m):
        i = 0
        while(i + 2**j < n-1 ):
            if(node_visited[rmq_array[i,j-1]].depth < node_visited[rmq_array[(i+2**(j-1)), j-1]].depth):
                rmq_array[i,j] = rmq_array[i,j-1]
            else :
                rmq_array[i,j] = rmq_array[i+2**(j-1), j-1]                 
            i += 1

    node_map = ddict()
    for i in xrange(n):
        cur_node = node_visited[i]
        # bad practice
        try:
            node_map[cur_node] = min(node_map[cur_node], i)
        except :
            node_map[cur_node] = i

    tree.add_features(lcaprocess=True)
    tree.add_features(rmqmat=rmq_array)
    tree.add_features(ind2node=node_visited)
    tree.add_features(node2ind=node_map)


def get_lca(sptree, species):
    """This should be a faster lcamapping
    species should be a list of node"""
    if not sptree.has_feature('lcaprocess', True):
        lcapreprocess(sptree)
    A = sptree.ind2node
    M = sptree.rmqmat
    # using the biggest interval should return the lca of all species
    s_index = sorted([sptree.node2ind[spec] for spec in species])
    #print "s_index vaut :", s_index , " et taille est : ", len(sptree.node2ind), " et rmq est : ", sptree.rmqmat.shape
    #print sptree.ind2node
    i = s_index[0]
    j = s_index[-1]
    k = int(np.log2(j-i+1))
    if (A[M[i, k]].depth<= A[M[j-2**(k) +1, k]].depth):
        return A[M[i,k]]
    else:
        return A[M[j- 2**(k) +1, k]]


def lcaMapping_old(geneTree, specieTree, multspeciename=True):

    smap = {}
    mapping = {}

    for node in geneTree.traverse(strategy="postorder"):
        try:
        # Children are always visited before parent
        # So for an internal node, we see leaves
        # under an internal node before visiting that node
        
            if not node.is_leaf():

                # ML ADDED THIS
                species = set([mapping[n] for n in node.get_children()])

                if(len(species) > 1):
                    mapping[node] = specieTree.get_common_ancestor(species)
                else:
                    mapping[node] = list(species)[0]

                if(multspeciename):
                    node.add_features(
                        species=",".join([x.name for x in species]))
                else:
                    node.add_features(species=mapping[node].name)

            else:
                sname = node.species
                if not sname in smap:
                    s = specieTree & node.species
                    smap[sname] = s
                else:
                    s = smap[sname]
                mapping[node] = s

        except Exception as e:
            print e
            print("Leaves without species")

    return mapping

def lcaMapping(genetree, specietree, multspeciename=True):
    """LCA mapping between a genetree and a specietree
    :argument genetree: your genetree, All leave in the genetree should already have feature 'specie' (set_specie was called)
    :argument specietree: your specietree
    :argument multspeciename: A flag to use in order to accept multi specie name at genetree internal node.
    """

    smap = {} #a dict that map specie name to specie node in specietree
    mapping = {}
    if not specietree.has_feature('lcaprocess', True):
        lcapreprocess(specietree)

    for node in genetree.traverse(strategy="postorder"):
        try:
        # Children are always visited before parent
        # So for an internal node, we see leaves
        # under an internal node before visiting that node
            if not node.is_leaf():
                # ML ADDED THIS
                species = set([mapping[n] for n in node.get_children()])
                if(len(species) > 1):
                    mapping[node] = get_lca(specietree, species)
                else:
                    mapping[node] = list(species)[0]

                if(multspeciename):
                    node.add_features(
                        species=",".join(sorted([x.name for x in species])))
                else:
                    node.add_features(species=mapping[node].name)

            else:
                sname = node.species
                if not sname in smap:
                    s = specietree & node.species
                    smap[sname] = s
                else:
                    s = smap[sname]
                mapping[node] = s

        except Exception as e:
            print e
            print("%s without species" % node.name)
            print node.get_ascii(attributes=['name', 'specie'])

    genetree.add_features(lcaMap=mapping)
    return mapping


def reconcile(genetree=None, lcaMap=None, lost=False):
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

        if (isinstance(lost, basestring) and lost.upper() == "YES") or lost==True:
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
                                set(lcaMap[intern_lost].get_leaf_names()) - set(child_c.species.split(",")))
                            splist = lostnode.species.split(',')
                            if(len(splist) > 1):
                                lostnode.name = "n" + \
                                    str(lost_count) + "_" + \
                                    "".join([s[0] for s in splist])

                            else:
                                lostnode.name = lostnode.species

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
                            lostnode.name = "node_" + \
                                "-".join([s[0] for s in unadded_specie]
                                         ) + "_" + str(lost_count)
                        else:
                            lostnode.name = lostnode.species

                        lost_count += 1
                        node.add_child(lostnode)
    genetree.add_features(reconciled=True)


def ComputeDLScore(genetree, lcaMap=None):
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
                node_is_dup = np.mean([params.getdup(l.name) for l in lcaMap[node].get_leaves()])
                dup_score += node_is_dup
            
            for child in node.get_children():
                if node_is_dup:
                    child_map = [lcaMap[node]]
                else:
                    child_map = lcaMap[node].get_children()
                curr_node = lcaMap[child]
                while(curr_node not in child_map):
                    lost_nodes = set(curr_node.up.get_leaves()) - set(curr_node.get_leaves())
                    loss_score += np.mean([params.getloss(l.name) for l in lost_nodes]) 
                    curr_node = curr_node.up 

    else :
        raise Exception("LcaMapping not provided !!")

    return dup_score, loss_score

def detComputeDLScore(genetree, lcaMap=None):
    """
    Compute the number of duplication and the number of losses
    """
    loss = 0
    dup = 0

    if not lcaMap and genetree.has_feature('lcaMap'):
        lcaMap = genetree.lcaMap

    if lcaMap:  
        for node in genetree.traverse("levelorder"):
            if (lcaMap[node] in [lcaMap[child] for child in node.get_children()]):
                dup += 1
            if(node.up):
                parent = node.up
                parent_is_dup = 0
                if(lcaMap[parent] in [lcaMap[child] for child in parent.get_children()]):
                    parent_is_dup = 1
                loss += (lcaMap[node].depth - lcaMap[parent].depth - 1 + parent_is_dup)
        print loss
    
    else :
        if(genetree is None or 'type' not in genetree.get_all_features()):
            raise Exception("LcaMapping not found and your Genetree didn't undergo reconciliation yet")

        for node in genetree.traverse():
            if node.has_feature('type'):
                if(node.type == TreeClass.NAD or node.type == TreeClass.AD ):
                    dup += 1
                elif node.type == TreeClass.LOST:
                    loss += 1

    return dup, loss


def CleanFeatures(tree=None, features=[]):
    cleaned = False
    if(tree):
        for node in tree.traverse():
            for f in features:
                if(node.has_feature(f)):
                    node.del_feature(f)
                    cleaned = True
    return cleaned


def binary_recon_score(node, lcamap):
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
            dup += 1

        children_list = node.get_children()
        supposed_children_species = lcamap[node].get_children_name()
        child_number = 0
        for child in children_list:
            c = lcamap[child]
            child_number += 1
            child_lost = 0

            if(dup == 0):
                while(c is not None and (c.name not in supposed_children_species)):
                    lost += 1
                    child_lost += 1
                    c = c.up

            if(dup > 0):
                while(c is not None and c.name != node.species):
                    lost += 1
                    child_lost += 1
                    c = c.up

    return dup * params.dupcost + lost * params.losscost, dup, lost


def totalDuplicationConsistency(tree):
    """Compute the total duplication consistency score for a tree"""
    dps = 0
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


def reset_node_name(tree, sep):
    for x in tree.traverse():
        x.name = x.name.split(sep)[0]
    return tree


def make_random_tree(names=list(string.lowercase), contract_seuil=0, feature_to_contract='support', random_branches=False):
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


def newick_preprocessing(newick, gene_sep=None):
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
                raise NewickError, \
                    'Unable to format your newick file, Bad gene-specie separator or too much special chars'
        nw += ';'
        return nw, gene_sep
    else:
        raise NewickError, \
            "'newick' argument must be either a filename or a newick string."


def polySolverPreprocessing(genetree, specietree, distance_file, capitalize=False, gene_sep=None, specie_pos="postfix", nFlagVal=1e305, nFlag=False, smap=None):
    """Preprocess genetree for polytomysolver
    """
    # genetree input
    speciemap = None
    if isinstance(genetree, basestring) and not smap:
        genetree, gene_sep = newick_preprocessing(genetree, gene_sep)
        genetree = TreeClass(genetree)

    elif smap:
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

    # specietree input
    if isinstance(specietree, basestring):
        specietree, sep = newick_preprocessing(specietree, '')
        specietree = TreeClass(specietree)
    specietree.label_internal_node()

    # distance matrice input
    if(distance_file):
        gene_matrix, node_order = clu.distMatProcessor(
            distance_file, nFlagVal, nFlag)
        # Difference check 1
        if set(node_order).difference(set(genetree.get_leaf_names())):
            reset_node_name(genetree, gene_sep)
    else:
        # This is for debug, will never happen
        print "error: dist file not found"
        node_order = genetree.get_leaf_names()
        # Alternative, retrieve aligned sequence and run phyML
        gene_matrix = clu.makeFakeDstMatrice(len(node_order), 0, 1)

    # Find list of species not in genetree
    specieGeneList = set(genetree.get_leaf_species())
    specieList = set([x.name for x in specietree.get_leaves()])
    if(specieGeneList - specieList):
        raise Exception("Species in genetree but not in specietree : %s" % (
            ", ".join(specieGeneList - specieList)))

    return genetree, specietree, gene_matrix, node_order


def export_to_orthoXML(t, database='customdb', handle=sys.stdout):
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

    # OrthoXML does not support duplication events to be at the root
    # of the tree, so we search for the top most speciation events in
    # the tree and export them as separate ortholog groups
    is_speciation = lambda n: getattr(n, 'type', "") == "S" or not n.children
    for speciation_root in t.iter_leaves(is_leaf_fn=is_speciation):
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
                        parent_event.add_paralogGroup(node2event[ch])
                    elif ch.type > 0:
                        parent_event.add_orthologGroup(node2event[ch])
                    else:
                        raise AttributeError(
                            "\n\Internals nodes labeled by losses are not expected in the orthoXML format")

    O.export(handle, 0, namespace_="")


def generate_smap(specietree, output="smap", relaxed=False, suffix=""):
    """
    Generate a specie map from genetree and specietree
    """
    gene_to_spec_map = []
    specie_names = specietree.get_leaf_names()
    for name in specie_names:
        if(relaxed):
            genes = re.compile(".*" + name + suffix + ".*", re.IGNORECASE)
        else :
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

    print "\nCorrected Tree binary list rf_fould\n"
    print "\n".join(map(lambda x: "\t".join([str(v) for v in x]), ct_binary))
    print "\nTree binary list rf_fould\n"
    print "\n".join(map(lambda x: "\t".join([str(v) for v in x]), t_binary))

    if(len(ct_success) == len(ct_leaves)):
        print "**Leave remaining success for corrected tree"
        print "\n".join([str(h) for h in t_success])

    else:
        print "**Corrected tree doesn't follow patern"
        print "\n".join(map(lambda x: "\t".join([str(v) for v in x]), ct_leaves))

    if(len(t_success) == len(t_leaves)):
        print "**Leave remaining success for tree"
        # print "\n".join([str(h) for h in t_success])
    else:
        print "**Tree doesn't follow patern"
        # print "\n".join(map(lambda x: "\t".join([str(v) for v in x]),
        # t_leaves))

    print "**Compatibility test between tree: ", all(success)
