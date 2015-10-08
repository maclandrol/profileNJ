# This file is part of profileNJ
#
# params set duplication and losses cost for profileNJ and reconciliation

__author__ = "Emmanuel Noutahi"

import hashlib, numpy as np

dupcost, losscost = 1, 1
internal = 0
def set(dup, loss, internal_mode='default'):
	global dupcost
	global losscost
	global internal
	dupcost, losscost = dup, loss
	internal = 1 if internal_mode=='mean' else 0


def get_hash(splist):
	if not isinstance(splist, basestring):
		splist = ",".join(sorted(splist))
	return hashlib.sha384(splist).hexdigest()


def getdup(specie=None):
	global dupcost
	slist = specie

	if specie and not isinstance(specie, basestring):
		slist = specie.get_leaf_names()

	if isinstance(dupcost, dict):
		return dupcost.get(get_hash(slist), get_internal(specie, getdup))
	else:
		return dupcost

def getloss(specie=None):
	global losscost
	slist = specie
	if specie and not isinstance(specie, basestring):
		slist = specie.get_leaf_names()
	if isinstance(losscost, dict):
		return losscost.get(get_hash(slist), get_internal(specie, getloss))
	else:
		return losscost

def get_internal(specie, costfun=getdup):
	global internal
	defcost = 1
	if not isinstance(specie, basestring) and specie.is_internal() and internal==1:
		defcost = np.mean([costfun(s) for s in specie.get_leaves()])
	return defcost