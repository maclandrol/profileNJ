# This file is part of profileNJ
#
# Date: 06/2017
# This file contains several implementation of advanced phylogenetic
# model and function to compute reconciliation likelihood

__author__ = "Emmanuel Noutahi"

import os
import re
import sys
import numpy as np
import random
from TreeClass import TreeClass
from collections import defaultdict as ddict
import itertools


def time_slice(tree, timeframes, leaves_as_extant=True, single_nodes=True):
    """Subdivize a tree into multiple time slice
    use `timeframes` as reference slices.
    `timeframes` can be a list of timestamp or a list of seed nodes.
    if `leaves_as_extant` is True, leaves will always be in the last timeframe.
    (always be the case for ultrametric trees).
    if `single_nodes` is True, nodes with only one child will be added as anchor
    for the slices
    """
    desc = tree.get_descendants()
    if 'brlen' not in tree.get_all_features():
        tree.compute_branches_length()
    tf = set([])
    rootinside = False
    for v in timeframes:
        if v == tree:
            rootinside = True
            tf.add(0.0)
        elif isinstance(v, tree.__class__) and v in desc:
            tf.add(v.brlen)
        elif isinstance(v, float):
            tf.add(v)
    if not rootinside:
        tf.add(0.0)
    tf = sorted(list(tf))
    _subdivize_by_timeframe(tree, tf, leaves_as_extant,
                            single_nodes, desc=desc)


def n_time_slice(tree, nsl, leaves_as_extant=True):
    """Subdivize a tree into `nsl` time slice of equal length
    if `leaves_as_extant` is True, leaves will always be in the last timeframe.
    (always be the case for ultrametric trees).
    """
    if 'brlen' not in tree.get_all_features():
        tree.compute_branches_length()
    nsl = int(nsl)
    assert nsl != 0, "nsl should be a positive integer"
    furthest_from_self = tree.get_farthest_leaf()
    max_dist = furthest_from_self[1]
    interval_t = max_dist / nsl
    tf = [i * interval_t for i in range(nsl)] + [max_dist]
    _subdivize_by_timeframe(tree, tf, leaves_as_extant)


def _subdivize_by_timeframe(tree, timeframes, leaves_as_extant=True, single_nodes=False, desc=None):
    """Subdivize tree into timeframes and add label to nodes
    """
    if not desc:
        desc = tree.get_descendants()
    # edit this to add single node
    cur_time = 0
    max_time = max(timeframes)
    tree.add_features(timestamp=cur_time)
    tree.add_features(time=max_time - timeframes[cur_time])
    sorted_desc = sorted(desc, key=lambda x: x.brlen)
    cur_time = 1
    for node in sorted_desc:
        while cur_time < len(timeframes) and node.brlen > timeframes[cur_time]:
            cur_time += 1
        cur_time = min(len(timeframes) - 1, cur_time)
        node.add_features(timestamp=cur_time)
        node.add_features(time=timeframes[cur_time])
        if node.is_leaf() and leaves_as_extant:
            node.timestamp = len(timeframes) - 1
            node.time = timeframes[len(timeframes) - 1]
        # reverse time to count from leaves here
        node.time = max_time - node.time

    if single_nodes:
        for node in sorted_desc:
            # here we will try to add single node in the tree
            # this is a precaution since TreeClass appear to not
            # be the same type as node
            # don't ask me why...
            ini_tstamp = node.timestamp
            i = -node.up.timestamp + ini_tstamp
            while i > 1:
                i -= 1
                internode = node.__class__()
                node.insert_node_between(node.up, internode)
                internode.add_features(timestamp=(ini_tstamp - i))
                internode.add_features(name=" t%d_  " % internode.timestamp)
                internode.add_features(brlen=timeframes[internode.timestamp])
                node.dist = abs(internode.brlen - node.brlen)
                internode.add_features(
                    dist=(internode.brlen - internode.up.brlen))
                internode.time = max_time - timeframes[internode.timestamp]


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
    rate_var = duration * roeotroe
    if rate_var > 0.0:
        # Kishino, Thorne and Bruno corrected the tendency for the rate to
        #   increase seen in teh TKP, 1998 model
        mu = np.log(starting_rate) - (rate_var / 2.0)
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
    if node.up is not None and node.up.has_feature('rate'):
        starting_rate = node.up.rate
    else:
        starting_rate = kwargs.get('starting_rate', None)
    roeotroe = kwargs.get('roeotroe', None)
    min_rate = kwargs.get('min_rate', None)
    max_rate = kwargs.get('max_rate', None)

    if starting_rate is None or roeotroe is None:
        raise ValueError("starting_rate and roeotroe should not be None")
    if roeotroe * duration <= 0.0:
        if (min_rate and starting_rate < min_rate) or (max_rate and starting_rate > max_rate):
            raise ValueError(
                "Parent rate is out of bounds, but no rate change is possible")

    r = _calc_KTB_rate(starting_rate, duration, roeotroe)
    mr = (starting_rate + r) / 2.0
    if max_rate and r > max_rate:
        assert(starting_rate <= max_rate)
        p_changing = (max_rate - starting_rate) / (r - starting_rate)
        mean_changing = (starting_rate + max_rate) / 2.0
        mr = p_changing * mean_changing + (1.0 - p_changing) * max_rate
        r = max_rate
    elif min_rate and r < min_rate:
        assert(starting_rate >= min_rate)
        p_changing = (starting_rate - min_rate) / (starting_rate - r)
        mean_changing = (starting_rate + min_rate) / 2.0
        mr = p_changing * mean_changing + (1.0 - p_changing) * min_rate
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
    """Transform a tree into an ultrametric one (ultrametric)
    """
    height = 0.0
    cs = tree.get_children()
    if len(cs) == 0:
        return 0
    child_height = []
    for i, child in enumerate(cs):
        child_val = make_clock_like(child) + child.dist
        child_height.append(child_val)
        height += child_val

    height /= len(cs)
    for i, child in enumerate(cs):
        scale_subtree_branches(child, height / child_height[i])

    return height


def scale_subtree_branches(node, factor):
    """Scale the branches lenght to its parent of a node
    by factor
    """
    old_dist = node.dist
    node.dist = old_dist * factor
    for child in node.get_children():
        scale_subtree_branches(child, factor)


def set_height_on_tree(tree):
    """Set height on a tree according to
    http://www.pnas.org/content/suppl/2012/10/04/1202997109.DCSupplemental/sapp.pdf
    """
    N = len(tree)
    for node in tree.traverse():
        if node.has_feature('timestamp'):
            node.add_features(height=1.0 - (node.timestamp + 1.0) / N)


def scale_tree_height(tree, scaledist=True):
    """Scale a tree heigh to a range from 0 to 1
    0 correspond to leaves and 1 to the root
    this change the dist feature"""
    h = tree.time
    for node in tree.traverse("postorder"):
        node.time /= h
        if node.up is not None and scaledist:
            node.dist /= h


def validate_speciation(tree):
    """Check if there are speciation occuring
    at the same time and slightly change the speciation time.
    This should be perfomed before time slicing.
    """
    tmp_d = ddict(int)
    for node in tree.get_descendants():
        if not node.is_leaf():
            # increase by 1e-5 each
            it_blen = tmp_d.get(node.brlen, 0)
            tmp_d[node.brlen] += 1

            node.brlen += it_blen * 1e-5
            node.dist += it_blen * 1e-5
            tmp_d[node.brlen] = 1


def list_time_slice(tree):
    slice_by_rank = ddict(list)
    node_data = []
    i = 0
    for node in tree.traverse("postorder"):
        slice_by_rank[node.timestamp].append(i)
        node.add_features(edge_i=i)
        node_data.append(node)
        i += 1
    return slice_by_rank, node_data


def discretize_timeframe(tree, size=10):
    """Slice a time interval into `size` smaller interval
    by adding single node in the tree"""
    t = tree.copy()
    nodelist = t.get_descendants("postorder")
    for node in nodelist:
        i = 1
        parent_len = node.up.brlen
        while i <= size:
            i += 1
            internode = node.__class__()
            node.insert_node_between(node.up, internode)
            internode.add_features(brlen=node.brlen * (i + 1) - i * parent_len)
            node.dist = abs(internode.brlen - node.brlen)
            internode.add_features(dist=(internode.brlen - internode.up.brlen))
    return t
