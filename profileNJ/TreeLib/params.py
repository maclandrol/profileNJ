# This file is part of profileNJ
#
# params set duplication and losses cost for profileNJ and reconciliation

__author__ = "Emmanuel Noutahi"

import hashlib
import numpy as np

cdup, closs = 1, 1
dupcost, losscost = {}, {}
internal_type = 0


def set(dup, loss, constdlcost=(1, 1), internal_mode='default'):
    global dupcost, losscost
    global cdup, closs
    global internal_type
    dupcost, losscost = dup, loss
    cdup, closs = constdlcost
    internal_type = 1 if internal_mode == 'mean' else 0


def get_hash(splist):
    if not isinstance(splist, basestring):
        splist = ",".join(sorted(splist))
    return hashlib.sha384(splist).hexdigest()


def getdup(specie=None):
    global dupcost, cdup
    slist = specie
    if specie and not isinstance(specie, basestring):
        slist = specie.get_leaf_names()
    if len(slist) > 1:
        return get_internal(specie, getdup, ctype="dup")
    else:
        return dupcost.get(get_hash(slist), cdup)


def getloss(specie=None):
    global losscost, closs
    slist = specie
    if specie and not isinstance(specie, basestring):
        slist = specie.get_leaf_names()
    if len(slist) > 1:
        return get_internal(specie, getloss)
    else:
        return losscost.get(get_hash(slist), closs)


def get_internal(specie, costfun, ctype="loss"):
    global internal_type, closs, cdup
    if not isinstance(specie, basestring) and (specie.is_internal() or specie.is_root()) and internal_type == 1:
        defcost = np.mean([costfun(s) for s in specie.get_leaves()])
    else:
        defcost = closs if ctype == 'loss' else cdup
    return defcost
