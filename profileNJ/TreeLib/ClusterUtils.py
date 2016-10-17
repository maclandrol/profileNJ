# This file is part of profileNJ
#
# Date: 02/2014
# ClusterUtils contain implementation of nj and upgma clustering algo, and
# required methods

__author__ = "Emmanuel Noutahi"

from TreeClass import TreeClass
import os
import numpy as np
from StringIO import StringIO
import random
try:
    from lxml import etree
    # should work since etree is used by ete
except ImportError:
    try:
        import xml.etree.cElementTree as etree
    except ImporError:
        try:
            import xml.etree.ElementTree as etree
        except:
            pass


np.set_printoptions(precision=3)
numerictypes = np.core.numerictypes.sctype2char
Float = numerictypes(float)


def find_smallest_index(matrice):
    """Return smallest number i,j index in a matrice
    A Tuple (i,j) is returned.
    Warning, the diagonal should have the largest number so it will never be choose
    """

    index = np.tril_indices_from(matrice, -1)
    return np.vstack(index)[:, matrice[index].argmin()]


def condense_matrix(matrice, smallest_index, method='upgma'):
    """Matrice condensation in the next iteration

    Smallest index is returned from find_smallest_index.
    For both leaf at i and j a new distance is calculated from the average of the corresponding
    row or the corresponding columns
    We then replace the first index (row and column) by the average vector obtained
    and the second index by an array with large numbers so that
    it is never chosen again with find_smallest_index.
    Now the new regroupement distance value is at the first position! (on row and column)
    """
    first_index, second_index = smallest_index
    # get the rows and make a new vector by updating distance
    rows = np.take(matrice, smallest_index, 1)
    # default we use upgma
    if(method.lower() == 'nj'):
        new_vector = (
            np.sum(rows, 1) - matrice[first_index, second_index]) * 0.5

    else:
        new_vector = np.average(rows, 1)

    # replace info in the row and column for first index with new_vector
    matrice[second_index] = new_vector
    matrice[:, second_index] = new_vector
    np.fill_diagonal(matrice, 0)
    # replace the info in the row and column for the second index with
    # high numbers so that it is ignored
    return remove_ij(matrice, first_index, first_index)


def remove_ij(x, i, j):
    # Row i and column j divide the array into 4 quadrants
    y = x[:-1, :-1]
    y[:i, j:] = x[:i, j + 1:]
    y[i:, :j] = x[i + 1:, :j]
    y[i:, j:] = x[i + 1:, j + 1:]
    return y


def calculate_Q_ij(matrice, ind, n):
    """Calcutates Q_matrix for two taxa

    """
    return (n - 2) * matrice[ind] - np.sum(matrice[ind[0]]) - np.sum(matrice[ind[1]])


def calculate_Q_matrix(matrice):
    """Calculate Q_matrix for nj algorithm

    """
    n = matrice.shape[0]
    Q_matrix = np.zeros(shape=matrice.shape)
    it = np.nditer(matrice, flags=['multi_index'])
    while not it.finished:
        ind = it.multi_index
        Q_matrix[ind] = calculate_Q_ij(matrice, ind, n)
        it.iternext()

    return Q_matrix


def paired_node_distance(matrice, smallest_index):
    i, j = smallest_index
    # i, j are the index of the recently joined node
    n = matrice.shape[0]

    # http://en.wikipedia.org/wiki/Neighbor_joining#equation_2
    # distance from the pair members to the new node second term
    x = np.sum(matrice[i]) - np.sum(matrice[:, j])

    if(n - 2 > 0):
        dist_i = 0.5 * matrice[i, j] + ((0.5 / (n - 2)) * (x))
        dist_j = matrice[i, j] - dist_i
        return dist_i, dist_j

    else:
        # We have only two node to join (final join)
        # Find the index of the node not already joined

        distance = matrice[i, j]
        # In this case, we split the dist value by two
        return distance / 2.0, distance / 2.0


def condense_node_order(matrice, smallest_index, node_order, method='upgma'):
    """
    condenses two nodes in node_order based on smallest_index info
    This function is used to create a tree while condensing a matrice
    with the condense_matrix function. The smallest_index is retrieved
    with find_smallest_index. The first index is replaced with a node object
    that combines the two nodes corresponding to the indices in node order.
    The second index in smallest_index is replaced with None.
    Also sets the branch length of the nodes to 1/2 of the distance between
    the nodes in the matrice"""
    index1, index2 = smallest_index
    node1 = node_order[index1]
    node2 = node_order[index2]
    # get the distance between the nodes and assign 1/2 the distance to the
    # Length property of each node

    if(method.lower() == 'nj'):
        dist = paired_node_distance(matrice, smallest_index)

    elif(method.lower() == 'upgma'):
        distance = matrice[index1, index2]
        dist = (distance / 2.0, distance / 2.0)

    else:
        dist = (0, 0)

    nodes = [node1, node2]
    pos = [0, 1]

    for ind in pos:
        nodes[ind].add_features(length=dist[ind])
    # combine the two nodes into a new TreeNode object
    new_node = TreeClass()
    new_node.add_child(node1)
    new_node.add_child(node2)
    new_node.add_features(length=sum(dist))
    # replace the object at index1 with the combined node
    node_order[index2] = new_node
    # replace the object at index2 with None
    del node_order[index1]  # distance at i=index2 || j=index2
    return node_order


def NJ_cluster(matrice, node_order, nj_depth=None):
    """
    Node clustering with NJ
    matrice is a np array.
    node_order is a list of PhyloNode objects corresponding to the matrice.

    WARNING: Changes matrice in-place.
    before this function is called.
    """

    # this is for a test, should made it into one function with upgma
    num_entries = len(node_order)
    if not nj_depth or nj_depth > (num_entries - 1):
        nj_depth = num_entries - 1  # default, do all, same as upgma

    tree = None
    smallest_index = []
    for i in range(nj_depth):
        Q_matrix = calculate_Q_matrix(matrice)
        index_1, index_2 = find_smallest_index(Q_matrix)
        smallest_index = (index_1, index_2)
        row_order = condense_node_order(
            matrice, smallest_index, node_order, method='nj')
        matrice = condense_matrix(matrice, smallest_index, method='nj')
        tree = node_order[smallest_index[1]]
    return tree, matrice, smallest_index


def UPGMA_cluster(matrice, node_order, upgma_depth=None):
    """cluster with UPGMA
    matrice is a np array.
    node_order is a list of TreeClass objects corresponding to the matrice.

    WARNING: Changes matrice in-place.
    before this function is called.
    """
    num_entries = len(node_order)
    if not upgma_depth or upgma_depth > (num_entries - 1):
        upgma_depth = num_entries - 1  # default, do all
    tree = None
    smallest_index = []
    for i in range(upgma_depth):
        index_1, index_2 = find_smallest_index(matrice)
        smallest_index = (index_1, index_2)
        assert(index_1 > index_2)
        row_order = condense_node_order(
            matrice, smallest_index, node_order, method='upgma')
        matrice = condense_matrix(matrice, smallest_index, method='upgma')
        tree = node_order[smallest_index[1]]
    return tree, matrice, smallest_index


def RAND_cluster(matrice, node_order, rand_depth=None):
    """Random clustering
    matrice is a np array.
    node_order is a list of PhyloNode objects corresponding to the matrice.s

    WARNING: Changes matrice in-place.
    before this function is called.
    """
    num_entries = len(node_order)
    if not rand_depth or rand_depth > (num_entries - 1):
        rand_depth = num_entries - 1  # default, do all
    tree = None
    smallest_index = []

    for i in range(rand_depth):
        tochoose = [i for i, t in enumerate(node_order) if t is not None]
        index1, index2 = random.sample(tochoose, 2)
        smallest_index = (max(index1, index2), min(index1, index2))
        node_order = condense_node_order(
            matrice, smallest_index, node_order, method='rand')
        tree = node_order[smallest_index[1]]

    return tree, matrice, smallest_index


def treeCluster(matrice, node_order, depth=None, method='upgma'):

    if(len(node_order) == 2):
        smallest_index = (1, 0)
        row_order = condense_node_order(
            matrice, smallest_index, node_order, method='rand')
        tree = node_order[smallest_index[1]]
        return tree, None, smallest_index

    if(method.lower() == 'nj'):
        return NJ_cluster(matrice, node_order, nj_depth=depth)
    elif(method.lower() == 'rand'):
        return RAND_cluster(matrice, node_order, rand_depth=depth)
    else:
        return UPGMA_cluster(matrice, node_order, upgma_depth=depth)


def distMatProcessor(distances, nFlagVal=1e305, nFlag=False, ignoreNodes=[]):
    """Formating distance matrix from a file or string input and node order for
        UPGMA or NJ join
    """

    read_fl = False
    dist_matrix = []
    node_order = []
    matrix = None

    # Read in matrix if file name is given
    if isinstance(distances, basestring) and os.path.exists(distances):
        distances = open(distances, 'rU')

    distances = distances.read()

    distances_lines = distances.splitlines()
    if '<?xml' in distances_lines[0]:
        # this is an xml file
        # parse it differently
        matrix, node_order = parseFastPhyloXml(
            StringIO(distances), nFlagVal, nFlag)
    else:
        x_ind = 0
        for line in distances_lines:
            line = line.strip()
            if(line):
                if not read_fl:
                    read_fl = True
                else:
                    x_ind += 1
                    line_list = [getFloatValue(
                        x.strip(), x_ind, y_ind, nFlagVal, nFlag) for y_ind, x in enumerate(line.split())]
                    dist_matrix.append(line_list[1:])
                    node_order.append(line_list[0])
        matrix = np.array(dist_matrix, dtype=np.float)

    if ignoreNodes:
        for n in ignoreNodes:
            ind = node_order.index(n)
            if ind > -1:
                matrix = remove_ij(matrix, ind, ind)
                node_order.remove(n)

    return matrix, node_order


def makeFakeDstMatrice(n, dmin, dmax):
    """Create a fake distance matrice"""
    b = (dmax - dmin) * np.random.random_sample(size=(n, n)) + dmin
    b_sym = (b + b.T) / 2
    np.fill_diagonal(b_sym, 0)
    return b_sym


def saveMatrix(filename, matrix, node_order):
    # matrix[np.where(matrix==1e305)]=0
    with open(filename, 'w+') as out:
        out.write("\t%i\n" % len(node_order))
        lines = []
        for entry in matrix.tolist():
            line = node_order.pop(0) + "\t" + " ".join(map(str, entry)) + "\n"
            lines.append(line)
        out.writelines(lines)
    return True


def getFloatValue(number, x_ind, y_ind, nFlagVal, nFlag=False):
    """Get a distance matrice validate input from a string"""
    try:
        n = float(number)
        if(n < 0 and nFlag):
            n = nFlagVal
        return 0 if (x_ind == y_ind) else n

    except ValueError:
        return number


def parseFastPhyloXml(infile, nFlagVal, nFlag=False):
    """Parse the fastphylo xml format"""
    xml = etree.parse(infile)
    run = xml.find('//run')
    dimension = int(run.attrib['dim'])
    identities = run.find('identities')
    node_order = [i.attrib['name'] for i in identities.iter('identity')]
    dm = run.find('dms').find('dm')
    distance_mat = np.zeros(shape=(dimension, dimension), dtype=np.float)
    i = 0
    for node in dm.iter('row'):
        j = 0
        for entry in node.iter('entry'):
            val = float(entry.text)
            if(val < 0 and nFlag):
                val = nFlagVal
            distance_mat[i, j] = val
            distance_mat[j, i] = val
            j += 1
        i += 1
    return distance_mat, node_order
